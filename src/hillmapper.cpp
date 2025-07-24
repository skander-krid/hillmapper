#include <cmath>
#include <functional>
#include <limits>
#include <memory>
#include <mutex>
#include <optional>
#include <random>
#include <thread>
#include <unordered_map>
#include <utility>
#include <vector>

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;

using TableIndex = std::size_t;

using Distance = double;

using RawScanset = std::vector<TableIndex>; // A scanset is a vector of table indices

using Bitmap = std::vector<std::uint64_t>; // Bitmap representation of a scanset

struct Scanset {
    const Bitmap bitmap; // Bitmap representation of a scanset
    const RawScanset tables; // We can assume that all tables are distinct!!!
    size_t n_conflicts = 0; // Number of conflicts in the scanset. We need this when translating scansets.

    Scanset() = default;

    Scanset(const Bitmap&& bm, const RawScanset&& tbl, size_t n_conflicts = 0)
    : bitmap(std::move(bm)), tables(std::move(tbl)), n_conflicts(n_conflicts) {};
};

using Assignment = std::unordered_map<TableIndex, TableIndex>;

Assignment get_random_assignment(TableIndex n, TableIndex k,
                                 std::optional<std::uint64_t> seed = std::nullopt)
{
	static std::mt19937_64 rng; // shared engine
	if (seed.has_value())
		rng.seed(*seed); // deterministic when requested

	std::uniform_int_distribution<TableIndex> dist(1, k);

	Assignment assignment;
	for (TableIndex i = 1; i <= n; ++i)
		assignment[i] = dist(rng);

	return assignment;
}

Bitmap make_bitmap(const RawScanset& X, TableIndex num_tables)
{
    const std::size_t num_words = (num_tables + 64) / 64; // round up
    Bitmap bitmap(num_words, 0); // bitmap must be zero-initialised

    for (TableIndex v : X)
    {
        std::size_t word = v >> 6; // v / 64
        std::size_t bit  = v & 0x3F; // v % 64
        bitmap[word] |= (std::uint64_t(1) << bit);
    }

    return bitmap;
}

Scanset convert_scanset_for_distance(
    const Scanset &scanset,
    const Assignment &assignment,
    size_t num_tables_in_scanset_bitmap)
{
	// For each tpcds table count the number of redset tables that map to it.
	std::unordered_map<TableIndex, size_t> bt_counts;
    size_t n_conflicts = 0;

    // Iterate over the scanset bitmap bit-by-bit
	for (TableIndex r_table : scanset.tables)
	{
		TableIndex b_table = assignment.at(r_table);
		bt_counts[b_table]++;
	}

	// Produce the converted scanset. Iterate over bt_counts.
	RawScanset result;
    for (const auto& [b_table, count] : bt_counts) {
        result.push_back(b_table);
        if (count > 1) {
            n_conflicts += count - 1; // Count conflicts
        }
    }

	return Scanset(std::move(make_bitmap(result, num_tables_in_scanset_bitmap)), std::move(result), n_conflicts);
}

Distance distance(const Scanset &scanset_1, const Scanset &scanset_2)
{
    // Usually only one or very few values.
    // If the whole workload references no more than 64 distinct tables,
    // then any scanset will fit into a single 64-bit integer.
    // So no SIMD needed.

    auto& a = scanset_1.bitmap;
    auto& b = scanset_2.bitmap;

	Distance dist = 0;
    for (size_t i = 0; i < std::min(a.size(), b.size()); i++)
    {
        // XOR the two scansets and count the number of bits set to 1.
        // This is the Hamming distance.
        dist += __builtin_popcountll(a[i] ^ b[i]);
    }
    // If one scanset is longer, add the remaining bits.
    if (a.size() < b.size())
    {
        for (size_t i = a.size(); i < b.size(); i++)
        {
            dist += __builtin_popcountll(b[i]);
        }
    }
    else if (b.size() < a.size())
    {
        for (size_t i = b.size(); i < a.size(); i++)
        {
            dist += __builtin_popcountll(a[i]);
        }
    }
    return dist + scanset_1.n_conflicts + scanset_2.n_conflicts;
}

// BK-tree implementation to speed up nearest-scanset-neighbor search.
// (logarithmic time)
class BKTree
{
private:
	struct Node
	{
		Scanset item;
		std::unordered_map<Distance, std::unique_ptr<Node>> children;

		explicit Node(const Scanset &item_) : item(item_) {}
	};

	std::unique_ptr<Node> root;

public:
	BKTree() = default;

	explicit BKTree(const std::vector<Scanset> &scansets)
	{
		for (const auto &s : scansets)
			add(s);
	}

	void add(const Scanset &item)
	{
		if (!root)
		{
			root = std::make_unique<Node>(item);
			return;
		}

		Node *node = root.get();
		while (true)
		{
			Distance  d  = distance(item, node->item);
			auto it = node->children.find(d);
			if (it != node->children.end())
			{
				node = it->second.get(); // descend
			}
			else
			{
				node->children[d] = std::make_unique<Node>(item); // create new child
				break;
			}
		}
	}

	// Returns {nearest-item, distance}.
	std::pair<const Scanset*, Distance> nearest(const Scanset &query) const
	{
        Distance best_distance = std::numeric_limits<Distance>::max();
        const Scanset* best_item = nullptr;
		if (!root)
			return {best_item, best_distance};

		std::function<void(const Node *)> search = [&](const Node *node)
		{
			Distance d = distance(query, node->item);
			if (d < best_distance) {
                best_distance = d;
                best_item = &node->item;
            }

			for (const auto &[child_dist, child] : node->children)
			{
				if (std::abs(d - best_distance) <= child_dist && child_dist <= d + best_distance)
					search(child.get());
			}
		};

		search(root.get());
		return {best_item, best_distance};
	}
};

// Convenience: vector with indices 0..N-1
inline std::vector<std::size_t> full_rs_index_range(size_t n)
{
	std::vector<std::size_t> idx(n);
	std::iota(idx.begin(), idx.end(), 0);
	return idx;
}

Distance get_assignment_distance(std::vector<size_t> &redset_scanset_counters,
                               std::vector<Scanset> &redset_scansets, const Assignment &assignment,
                               const std::vector<std::size_t> &affected_r_scansets, // Their indexes
                               const BKTree &bk_tree, std::vector<Distance> &per_rs_distance,
                               size_t num_tables_in_scanset_bitmap,
                               bool compute_per_rs_distance, bool update_per_rs_distance = false
                               )
{
	Distance total_distance = 0.0;
	Distance delta          = 0.0;

	for (std::size_t r_scanset_idx : affected_r_scansets)
	{
		Scanset &r_scanset = redset_scansets[r_scanset_idx];
		std::size_t occs = redset_scanset_counters[r_scanset_idx];

		Distance new_distance =
		    bk_tree.nearest(convert_scanset_for_distance(r_scanset, assignment, num_tables_in_scanset_bitmap)).second *
		    static_cast<Distance>(occs);

		if (compute_per_rs_distance)
		{
			total_distance += new_distance;
			per_rs_distance.push_back(new_distance);
		}
		else
		{
			Distance old_distance = per_rs_distance[r_scanset_idx];
			delta += new_distance - old_distance;
			if (update_per_rs_distance)
				per_rs_distance[r_scanset_idx] = new_distance;
		}
	}
	return compute_per_rs_distance ? total_distance : delta;
}

std::pair<Distance, Assignment> thread_work(
    std::size_t n_iterations, std::size_t thread_idx, const BKTree &bk_tree,
    std::vector<size_t>              redset_scanset_counters,
    std::vector<Scanset>             redset_scansets, // TODO: MUST ALL BE DISTINCT!!!!
    std::vector<Scanset>             tpcds_scansets,
    std::vector<std::vector<size_t>> r_table_to_scansets, // Maps redset table to list of
                                                          // scanset indexes it affects
    TableIndex n_redset_tables, TableIndex n_tpcds_tables,
    size_t num_tables_in_scanset_bitmap)
{
	Distance     best_distance = std::numeric_limits<Distance>::infinity();
	Assignment best_assignment{};
	const auto all_rs_indices = full_rs_index_range(redset_scansets.size());

	for (std::size_t itr = 0; itr < n_iterations; itr++)
	{
		Assignment assignment = get_random_assignment(
		    n_redset_tables, n_tpcds_tables, static_cast<std::uint64_t>(thread_idx * 1000 + itr));
		std::vector<Distance> current_per_rs_distance;
		Distance              current_distance =
		    get_assignment_distance(redset_scanset_counters, redset_scansets, assignment,
		                            all_rs_indices, bk_tree, current_per_rs_distance,
                                    num_tables_in_scanset_bitmap,
		                            /*compute_per_rs_distance=*/true); // TODO: Sth missing?
		while (true)
		{
			Distance best_new_distance   = std::numeric_limits<Distance>::infinity();
			std::pair<int, int> best_new_assignment = {-1, -1};
            // TODO: Randomize the order in which I traverse rt and bt to avoid repetitive assignments.
			for (TableIndex rt = 1; rt <= n_redset_tables; rt++)
			{
				TableIndex old_assignment = assignment[rt];
				for (TableIndex bt = 1; bt <= n_tpcds_tables; bt++)
				{
					if (bt == old_assignment)
						continue; // skip the current assignment

					assignment[rt] = bt; // try a new assignment
					Distance new_distance =
					    current_distance + get_assignment_distance(redset_scanset_counters, redset_scansets,
					                                               assignment, r_table_to_scansets[rt],
					                                               bk_tree, current_per_rs_distance,
                                                                   num_tables_in_scanset_bitmap,
					                                               /*compute_per_rs_distance=*/false);

					if (new_distance < best_new_distance)
					{
						best_new_distance   = new_distance;
						best_new_assignment = {rt, bt};
					}
				}
				assignment[rt] = old_assignment; // revert the assignment
			}
            // TODO:!!! When I have the choice over which bt to choose, prefer an unassigned one.
            // To minimize the number of times I suplicate tables later in workload generation.
			if (best_new_distance >= current_distance)
			{
				break; // Local optimum reached
			}

			assignment[best_new_assignment.first] =
			    best_new_assignment.second; // Set the best new global assignment

			// update per-scanset distances in place
			get_assignment_distance(redset_scanset_counters, redset_scansets, assignment, all_rs_indices,
			                        bk_tree, current_per_rs_distance,
                                    num_tables_in_scanset_bitmap,
			                        /*compute_per_rs_distance=*/false,
			                        /*update_per_rs_distance=*/true);

			current_distance = best_new_distance;
		}

		if (current_distance < best_distance)
		{
			best_distance   = current_distance;
			best_assignment = assignment;
		}
	}
	return {best_distance, best_assignment};
}

std::pair<Distance, Assignment> find_optimal_bijection(
    std::size_t n_threads, std::size_t n_iterations_per_thread,
    const std::vector<std::size_t> &redset_scanset_counters,
    const std::vector<RawScanset> &raw_redset_scansets, const std::vector<RawScanset> &raw_tpcds_scansets,
    const std::vector<std::vector<std::size_t>> &r_table_to_scansets, TableIndex n_redset_tables,
    TableIndex n_tpcds_tables)
{
    const TableIndex num_tables_in_scanset_bitmap = std::max(n_redset_tables, n_tpcds_tables);

    // Convert raw scansets to Scanset objects
    std::vector<Scanset> redset_scansets;
    for (const auto &rs : raw_redset_scansets)
    {
        redset_scansets.emplace_back(std::move(make_bitmap(rs, num_tables_in_scanset_bitmap)), std::move(rs));
    }
    std::vector<Scanset> tpcds_scansets;
    for (const auto &ts : raw_tpcds_scansets)
    {
        tpcds_scansets.emplace_back(std::move(make_bitmap(ts, num_tables_in_scanset_bitmap)), std::move(ts));
    }

	// Build the shared BK-tree
	BKTree bk_tree(tpcds_scansets);

	// Start threads
	std::vector<std::thread> workers;
	std::mutex guard; // protects global best
	Distance best_distance = std::numeric_limits<Distance>::infinity();
	Assignment best_assignment;

	for (std::size_t tid = 0; tid < n_threads; ++tid)
	{
		workers.emplace_back(
		    [&, tid]()
		    {
			    auto [dist, assign] =
			        thread_work(n_iterations_per_thread, tid,
			                    bk_tree, // shared
			                    redset_scanset_counters, redset_scansets, tpcds_scansets,
			                    r_table_to_scansets, n_redset_tables, n_tpcds_tables, num_tables_in_scanset_bitmap);

			    std::lock_guard<std::mutex> lk(guard);
			    // TODO: We can do better (no locks) here, but probably not worth the
			    // effort.
			    if (dist < best_distance)
			    {
				    best_distance   = dist;
				    best_assignment = std::move(assign);
			    }
		    });
	}

	// Join threads and return result
	for (auto &t : workers)
		t.join();
	return {best_distance, best_assignment};
}

PYBIND11_MODULE(_hillmapper, m)
{
	m.def("find_optimal_bijection", &find_optimal_bijection, py::arg("n_threads"),
	      py::arg("n_iterations_per_thread"), py::arg("redset_scanset_counters"),
	      py::arg("redset_scansets"), py::arg("tpcds_scansets"), py::arg("r_table_to_scansets"),
	      py::arg("n_redset_tables"), py::arg("n_tpcds_tables"),
	      R"pbdoc(
            Parallel implementation of hill climbing to find the optimal
            assignment of redset tables to tpcds tables that minimizes the
            sum of the scansets manhattan distances.
            return the best (distance, assignment) pair.
        )pbdoc");
}
