#include <cmath>
#include <functional>
#include <limits>
#include <memory>
#include <optional>
#include <random>
#include <unordered_map>
#include <utility>
#include <vector>

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;

using Assignment = std::unordered_map<int, int>;

// TODO: Use different reperesentation: vector where value at index i is the
// number of occurences of table i in the scanset.
// Then for computing the distance, I can avoid any branching + vectorize.
// TODO: Or switch to scanSET semantics instead of scanMULTISETS. Then I can
// represent a scanset as a bitmap, and compute the distance with vectorized
// xor.
using Scanset = std::vector<int>;

Assignment
get_random_assignment(std::size_t n, std::size_t k,
                      std::optional<std::uint64_t> seed = std::nullopt) {
  static std::mt19937_64 rng; // shared engine
  if (seed.has_value())
    rng.seed(*seed); // deterministic when requested

  std::uniform_int_distribution<int> dist(1, static_cast<int>(k));

  Assignment assignment;
  for (int i = 1; i <= static_cast<int>(n); ++i)
    assignment[i] = dist(rng);

  return assignment;
}

Scanset convert_scanset_for_distance(const Scanset &scanset,
                                     const Assignment &assignment) {
  // table id → #occurrences
  std::unordered_map<int, int> ts;
  for (int r_table : scanset)
    ++ts[r_table];

  // Build collision map: b_table -> { all r_tables that map to it }
  std::unordered_map<int, std::unordered_set<int>> collisions;
  std::unordered_set<int> problematic_r_tables;

  for (int r_table : scanset) {
    int b_table = assignment.at(r_table);
    collisions[b_table].insert(r_table);
    problematic_r_tables.insert(r_table);
  }

  // Decide which r_tables survive when collisions occur
  std::unordered_set<int> chosen_r_tables;

  for (const auto &kv : collisions) {
    const auto &r_tables = kv.second;

    if (r_tables.size() == 1) {
      chosen_r_tables.insert(*r_tables.begin());
    } else {
      // Pick the r_table with the highest occurrence count in scanset
      int best_r_table = -1;
      int best_num_occs = 0;

      for (int r_table : r_tables) {
        int occs = ts[r_table];
        if (occs > best_num_occs) {
          best_num_occs = occs;
          best_r_table = r_table;
        }
      }
      chosen_r_tables.insert(best_r_table);
    }
  }

  // Produce the converted scanset
  Scanset result;
  result.reserve(scanset.size());

  for (int table : scanset) {
    if (!problematic_r_tables.count(table) // not a collision
        || chosen_r_tables.count(table))   // or chosen survivor
    {
      result.push_back(assignment.at(table));
    } else {
      result.push_back(-1); // masked out
    }
  }
  return result;
}

int distance(const Scanset &a, const Scanset &b) {
  // TODO: We can do better with sorted vector scansets and double pointer
  // approach.
  std::unordered_map<int, std::ptrdiff_t> delta;

  // Increment counts for the first scanset.
  for (int v : a)
    ++delta[v];

  // Decrement counts for the second scanset.
  for (int v : b)
    --delta[v];

  // Accumulate absolute differences.
  std::size_t dist = 0;
  for (const auto &[value, diff] : delta)
    dist += std::abs(diff);

  return dist;
}

// BK-tree implementation to speed up nearest-scanset-neighbor search.
// (logarithmic time)
class BKTree {
private:
  struct Node {
    Scanset item;
    std::unordered_map<int, std::unique_ptr<Node>> children;

    explicit Node(const Scanset &item_) : item(item_) {}
  };

  std::unique_ptr<Node> root;

public:
  BKTree() = default;

  explicit BKTree(const std::vector<Scanset> &scansets) {
    for (const auto &s : scansets)
      add(s);
  }

  void add(const Scanset &item) {
    if (!root) {
      root = std::make_unique<Node>(item);
      return;
    }

    Node *node = root.get();
    while (true) {
      int d = distance(item, node->item);
      auto it = node->children.find(d);
      if (it != node->children.end()) {
        node = it->second.get(); // descend
      } else {
        node->children[d] = std::make_unique<Node>(item); // create new child
        break;
      }
    }
  }

  // Returns {nearest-item, distance}.
  std::pair<Scanset, int> nearest(const Scanset &query) const {
    std::pair<Scanset, int> best{{}, std::numeric_limits<int>::max()};
    if (!root)
      return best;

    std::function<void(const Node *)> search = [&](const Node *node) {
      int d = distance(query, node->item);
      if (d < best.second)
        best = {node->item, d};

      for (const auto &[child_dist, child] : node->children) {
        if (std::abs(d - best.second) <= child_dist &&
            child_dist <= d + best.second)
          search(child.get());
      }
    };

    search(root.get());
    return best;
  }
};

// Convenience: vector with indices 0..N-1
inline std::vector<std::size_t> full_rs_index_range(size_t n) {
  std::vector<std::size_t> idx(n);
  std::iota(idx.begin(), idx.end(), 0);
  return idx;
}

double get_assignment_distance(
    std::vector<size_t> &redset_scanset_counters,
    std::vector<Scanset> &redset_scansets, const Assignment &assignment,
    const std::vector<std::size_t> &affected_r_scansets, // Their indexes
    const BKTree &bk_tree, std::vector<double> &per_rs_distance,
    bool compute_per_rs_distance, bool update_per_rs_distance = false) {
  double total_distance = 0.0;
  double delta = 0.0;

  for (std::size_t r_scanset_idx : affected_r_scansets) {
    Scanset &r_scanset = redset_scansets[r_scanset_idx];
    std::uint64_t occs = redset_scanset_counters[r_scanset_idx];

    double new_distance =
        bk_tree.nearest(convert_scanset_for_distance(r_scanset, assignment))
            .second *
        static_cast<double>(occs);

    if (compute_per_rs_distance) {
      total_distance += new_distance;
      per_rs_distance.push_back(new_distance);
    } else {
      double old_distance = per_rs_distance[r_scanset_idx];
      delta += new_distance - old_distance;
      if (update_per_rs_distance)
        per_rs_distance[r_scanset_idx] = new_distance;
    }
  }
  return compute_per_rs_distance ? total_distance : delta;
}

std::pair<double, Assignment> thread_work(
    std::size_t n_iterations, std::size_t thread_idx, const BKTree &bk_tree,
    std::vector<size_t> redset_scanset_counters,
    std::vector<Scanset> redset_scansets, // TODO: MUST ALL BE DISTINCT!!!!
    std::vector<Scanset> tpcds_scansets,
    std::vector<std::vector<size_t>>
        r_table_to_scansets, // Maps redset table to list of scanset indexes it
                             // affects
    size_t n_redset_tables, size_t n_tpcds_tables) {
  double best_distance = std::numeric_limits<double>::infinity();
  Assignment best_assignment{};
  const auto all_rs_indices = full_rs_index_range(redset_scansets.size());

  for (std::size_t itr = 0; itr < n_iterations; itr++) {
    Assignment assignment = get_random_assignment(
        n_redset_tables, n_tpcds_tables,
        static_cast<std::uint64_t>(thread_idx * 1000 + itr));
    std::vector<double> current_per_rs_distance;
    double current_distance = get_assignment_distance(
        redset_scanset_counters, redset_scansets, assignment, all_rs_indices,
        bk_tree, current_per_rs_distance,
        /*compute_per_rs_distance=*/true); // TODO: Sth missing?
    while (true) {
      double best_new_distance = std::numeric_limits<double>::infinity();
      std::pair<int, int> best_new_assignment = {-1, -1};
      for (int rt = 1; rt <= static_cast<int>(n_redset_tables); rt++) {
        int old_assignment = assignment[rt];
        for (int bt = 1; bt <= static_cast<int>(n_tpcds_tables); bt++) {
          if (bt == old_assignment)
            continue; // skip the current assignment

          assignment[rt] = bt; // try a new assignment
          double new_distance =
              current_distance +
              get_assignment_distance(redset_scanset_counters, redset_scansets,
                                      assignment, r_table_to_scansets[rt],
                                      bk_tree, current_per_rs_distance,
                                      /*compute_per_rs_distance=*/false);

          if (new_distance < best_new_distance) {
            best_new_distance = new_distance;
            best_new_assignment = {rt, bt};
          }
        }
        assignment[rt] = old_assignment; // revert the assignment
      }
      if (best_new_distance >= current_distance) {
        break; // Local optimum reached
      }

      assignment[best_new_assignment.first] =
          best_new_assignment.second; // Set the best new global assignment

      // update per-scanset distances in place
      get_assignment_distance(redset_scanset_counters, redset_scansets,
                              assignment, all_rs_indices, bk_tree,
                              current_per_rs_distance,
                              /*compute_per_rs_distance=*/false,
                              /*update_per_rs_distance=*/true);

      current_distance = best_new_distance;
    }

    if (current_distance < best_distance) {
      best_distance = current_distance;
      best_assignment = assignment;
    }
  }
  return {best_distance, best_assignment};
}

std::pair<double, Assignment> find_optimal_bijection(
    std::size_t n_threads, std::size_t n_iterations_per_thread,
    const std::vector<std::size_t> &redset_scanset_counters,
    const std::vector<Scanset> &redset_scansets,
    const std::vector<Scanset> &tpcds_scansets,
    const std::vector<std::vector<std::size_t>> &r_table_to_scansets,
    std::size_t n_redset_tables, std::size_t n_tpcds_tables) {
  // Build the shared BK-tree
  BKTree bk_tree(tpcds_scansets);

  // Start threads
  std::vector<std::thread> workers;
  std::mutex guard; // protects global best
  double best_distance = std::numeric_limits<double>::infinity();
  Assignment best_assignment;

  for (std::size_t tid = 0; tid < n_threads; ++tid) {
    workers.emplace_back([&, tid]() {
      auto [dist, assign] =
          thread_work(n_iterations_per_thread, tid,
                      bk_tree, // shared
                      redset_scanset_counters, redset_scansets, tpcds_scansets,
                      r_table_to_scansets, n_redset_tables, n_tpcds_tables);

      std::lock_guard<std::mutex> lk(guard);
      // TODO: We can do better (no locks) here, but probably not worth the
      // effort.
      if (dist < best_distance) {
        best_distance = dist;
        best_assignment = std::move(assign);
      }
    });
  }

  // Join threads and return result
  for (auto &t : workers)
    t.join();
  return {best_distance, best_assignment};
}

PYBIND11_MODULE(bktree, m) {
  py::class_<BKTree>(m, "BKTree")
      .def(py::init<>())                             // empty tree
      .def(py::init<const std::vector<Scanset> &>(), // from data
           py::arg("scansets"))
      .def("add", &BKTree::add, py::arg("item"))
      .def("nearest", &BKTree::nearest, py::arg("query"),
           R"pbdoc(
                 Return a tuple (nearest_scanset, distance).
             )pbdoc");

  m.def("find_optimal_bijection", &find_optimal_bijection, py::arg("n_threads"),
        py::arg("n_iterations_per_thread"), py::arg("redset_scanset_counters"),
        py::arg("redset_scansets"), py::arg("tpcds_scansets"),
        py::arg("r_table_to_scansets"), py::arg("n_redset_tables"),
        py::arg("n_tpcds_tables"),
        R"pbdoc(
            Parallel implementation of hill climbing to find the optimal
            assignment of redset tables to tpcds tables that maximize the
            scanset matches.
            return the best (distance, assignment) pair.
        )pbdoc");
}
