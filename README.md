# HILLMAPPER C++ Extension

```sh
# before
c++ -O3 -Wall -shared -std=c++17 -fPIC $(python3 -m pybind11 --includes) bktree.cpp -o bktree$(python3-config --extension-suffix) $(python3-config --ldflags)
```