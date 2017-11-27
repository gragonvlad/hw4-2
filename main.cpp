#include <iostream>
#include <chrono>
#include <vector>
#include <algorithm>
#include <zconf.h>
#include <fcntl.h>


void drop_caches() {
    int cache_fd = open("/proc/sys/vm/drop_caches", O_WRONLY);
    const auto *d3 = reinterpret_cast<const char *>("3");
    write(cache_fd, d3, sizeof(char));
    fsync(cache_fd);
    close(cache_fd);
}

bool search_sqrt(const std::vector<int>& data, int elem) {
    std::vector<int> pivots;
    std::vector<size_t> idx;
    for(size_t i = 0; i < data.size(); i += sqrt(data.size())) {
        pivots.push_back(data[i]);
        idx.push_back(i);
    }
    idx.push_back(data.size());

    auto first = std::lower_bound(pivots.begin(), pivots.end(), elem);

    for(auto it = data.begin() + *(idx.begin() + (first - pivots.begin()));
        it != data.begin() + *(idx.begin() + (first - pivots.begin() - 1)); it--) {
        if (*it == elem)
            return true;
    }

    return false;
}

bool run_tests(const std::vector<int>& data, int elem, const std::string &type) {
    if (type == "normal") {
        return std::binary_search(data.begin(), data.end(), elem);
    } else if (type == "sqrt") {
        return search_sqrt(data, elem);
    } else if (type == "prefetch") {
        return false;
    } else if (type == "vab") {
        return false;
    } else {
        return false;
    }
}

int main() {
    std::vector<int> data{1, 2, 3, 4, 5, 6, 7, 8};
    std::cout << "Found norm: " << run_tests(data, 3, "normal") << std::endl;
    std::cout << "Found sqrt: " << run_tests(data, 5, "sqrt") << std::endl;
    return 0;
}
