#include <cmath>
#include <iostream>
#include <chrono>
#include <utility>
#include <vector>
#include <algorithm>
#include <zconf.h>
#include <fcntl.h>
#include <functional>
#include <iomanip>


namespace binsearch {
    bool binary_search(std::vector<unsigned long long int>::const_iterator begin, std::vector<unsigned long long int>::const_iterator end, unsigned long key)
    {

        while(begin < end)
        {
            auto mid = begin + (std::distance(begin, end) >> 1);

            if(*mid == key)
            {
                return true;
            }
            else if(*mid > key)
            {
                end = mid;
            }
            else
            {
                begin = mid + 1;
            }
        }

        return false;
    }

    bool prefetch_binary_search(std::vector<unsigned long long int>::const_iterator begin, std::vector<unsigned long long int>::const_iterator end, unsigned long key)
    {
        long dist = 0;
        std::vector<unsigned long long int>::const_iterator mid;
        while(begin < end)
        {
            dist = std::distance(begin, end) >> 1;
            mid = begin + dist;
            __builtin_prefetch(&(*(begin + (dist >> 1))), 0, 3);
            __builtin_prefetch(&(*(mid + 1 + (dist >> 1))), 0, 3);

            if(*mid == key)
            {
                return true;
            }
            else if(*mid > key)
            {
                end = mid;
            }
            else
            {
                begin = mid + 1;
            }
        }

        return false;
    }

    namespace veb {

        inline size_t power_of_two(size_t exponent) { return static_cast<size_t>(1 << exponent); }

        std::vector<std::vector<unsigned long long>> get_precalc(unsigned long long n) {
            auto h = static_cast<unsigned long long>(std::ceil(std::log((double)n + 1) / std::log(2.0)));
            std::vector<std::vector<unsigned long long>> res(h);
            std::vector<unsigned long long> tmp(4);
            for (size_t i = 1; i <= h; ++i) {
                tmp[0] = i >> 1;
                tmp[1] = power_of_two(tmp[0]) - 1;
                tmp[2] = power_of_two(i - tmp[0]) - 1;
                tmp[3] = power_of_two(tmp[0] - 1);
                res[i - 1] = tmp;
                tmp.reserve(4);
            }
            return res;
        }

        // Implicit van Emde Boas binary search
        long long veb_search(const std::vector<unsigned long long>::const_iterator veb_array,
                             size_t length, unsigned long long elt,
                             const std::vector<std::vector<unsigned long long>>& precalc) {
            unsigned long long D, subtree_size, subtree_leaf_count;
            auto h = static_cast<unsigned long long>(std::ceil(std::log((double)length + 1) / std::log(2.0))) - 1;
            D = precalc[h][1];
            subtree_size = precalc[h][2];
            subtree_leaf_count = precalc[h][3];

            if (length > 1) {
                // Recurse on top half of tree
                long long subtree_index = veb_search(veb_array, D, elt, precalc);
                if (subtree_index < 0) return subtree_index;

                unsigned long long offset = subtree_index * subtree_size + D;

                // If not in top half, use subtree index to find place in bottom half
                long long bottom_subtree_index = veb_search(veb_array + offset, subtree_size, elt, precalc);
                return static_cast<long long int>(subtree_leaf_count * subtree_index + bottom_subtree_index);
            } else {
                unsigned long long root = *veb_array;

                if (elt == root) {
                    return -1;
                }

                return (elt < root) ? 0 : 1;
            }
        }

        unsigned long long tree_size(unsigned long long depth) { return power_of_two(depth) - 1; }

        // Given a sorted index in the range [1, 2^{tree_height}-1], return the VEB address
        unsigned long long veb_index(unsigned long long n, unsigned long long tree_height) {
            if (tree_height <= 1) return n;

            // Chop n's bitstring into two halves
            unsigned long long bottom_half = tree_height >> 1;
            unsigned long long top_half = tree_height - bottom_half;

            // Store each half
            unsigned long long top_n = n >> bottom_half;
            unsigned long long bottom_n = n & (power_of_two(bottom_half) - 1);

            // Recurse
            if (bottom_n == 0) return veb_index(top_n, top_half);

            unsigned long long top_address = top_n * tree_size(bottom_half) + tree_size(top_half);
            unsigned long long bot_address = veb_index(bottom_n, bottom_half);

            return top_address + bot_address;
        }


        std::vector<unsigned long long> make_layout(const std::vector<unsigned long long>& data) {

            unsigned long long len = data.size();
            std::vector<std::pair<unsigned long long, unsigned long long> > indices(len);
            std::vector<unsigned long long> idxs(len);

            unsigned long long power = power_of_two(static_cast<unsigned long long>(ceil(std::log((float) len) / std::log(2.0))));
            unsigned long long i;

            for (i = 0; i < len; i++) {
                size_t idx = veb_index(i + 1, power);\
                indices[i] = std::make_pair(idx, data[i]);
            }

            sort(indices.begin(), indices.end());

            for (i = 0; i < len; i++) {
                idxs[i] = indices[i].second;
            }

            return idxs;
        }
    }

}

namespace rands {
    std::mt19937 get_generator() {
        std::random_device rd;
        std::mt19937 gen(rd());
        return gen;
    }

    std::vector<unsigned long long> generate_data(size_t size) {
        auto generator = get_generator();
        std::vector<unsigned long long> v(size);
        std::uniform_int_distribution<unsigned long long> nums(0, static_cast<unsigned long long>(1E12));
        std::function<unsigned long long()> number_generator = std::bind(nums, generator);
        std::generate(v.begin(), v.end(), number_generator);
        std::sort(v.begin(), v.end());
        return v;
    }

    std::vector<unsigned long long> generate_requests(size_t req_count) {
        auto generator = get_generator();
        std::vector<unsigned long long> v(req_count);
        std::uniform_int_distribution<unsigned long long> nums(0, 1E12);
        std::function<unsigned long long()> number_generator = std::bind(nums, generator);
        std::generate(v.begin(), v.end(), number_generator);
        return v;
    }
}

namespace misc {
    std::vector<unsigned long long> read_data(size_t size) {
        std::vector<unsigned long long> v(size);
        for (size_t i = 0; i < size; ++i) {
            std::cin >> v[i];
        }
        return v;
    }

    void drop_caches() {
        int cache_fd = open("/proc/sys/vm/drop_caches", O_WRONLY);
        const auto *d3 = reinterpret_cast<const char *>("3");
        write(cache_fd, d3, sizeof(char));
        fsync(cache_fd);
        close(cache_fd);
    }

    template<class T>
    void vprint(std::vector<T> &d) {
        for (auto &e: d) {
            std::cout << e << " ";
        }
        std::cout << std::endl;
    }
}

namespace tests {

    typedef unsigned long long ull;
    typedef std::vector<ull> vull;
    typedef std::vector<ull>::iterator vullit;
    typedef std::vector<bool> vb;
    typedef std::vector<bool>::iterator vbit;

    typedef std::pair<std::pair<ull, double>, vb> timeres;

    timeres test_basic_self(const vull &data, const vull &requests) {
        std::chrono::time_point<std::chrono::system_clock> start, end;
        vb resullts(requests.size());
        start = std::chrono::system_clock::now();
        for (size_t i = 0; i < requests.size(); ++i) {
            resullts[i] = binsearch::binary_search(data.begin(), data.end(), requests[i]);
        }
        end = std::chrono::system_clock::now();
        ull total = static_cast<ull>(std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count());

        double per_request = total * 1. / requests.size();
        return std::make_pair(std::make_pair(total, per_request), resullts);
    };

    timeres test_basic(const vull &data, const vull &requests) {
        std::chrono::time_point<std::chrono::system_clock> start, end;
        vb resullts(requests.size());
        start = std::chrono::system_clock::now();
        for (size_t i = 0; i < requests.size(); ++i) {
            resullts[i] = std::binary_search(data.begin(), data.end(), requests[i]);
        }
        end = std::chrono::system_clock::now();
        ull total = static_cast<ull>(std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count());

        double per_request = total * 1. / requests.size();
        return std::make_pair(std::make_pair(total, per_request), resullts);
    };

    timeres test_prefetch(const vull &data, const vull &requests) {
        std::chrono::time_point<std::chrono::system_clock> start, end;
        vb resullts(requests.size());
        start = std::chrono::system_clock::now();
        for (size_t i = 0; i < requests.size(); ++i) {
            resullts[i] = binsearch::prefetch_binary_search(data.begin(), data.end(), requests[i]);
        }
        end = std::chrono::system_clock::now();
        ull total = static_cast<ull>(std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count());

        double per_request = total * 1. / requests.size();
        return std::make_pair(std::make_pair(total, per_request), resullts);
    };

    timeres test_sqrt(const vull &data, const vull &requests, ull n = 0) {
        if (!n) {
            n = data.size();
        }

        auto sqrt_n = static_cast<unsigned long long>(ceil(sqrt(n)));
        vull helper(static_cast<unsigned long>(ceil(n / sqrt_n)));

        for (size_t i = 0, j = 0; j < helper.size();
             i = sqrt_n + i >= n
                 ? n - 1
                 : i + sqrt_n,
             j++) {
            helper[j] = data[i];
        }

        misc::drop_caches();

        std::chrono::time_point<std::chrono::system_clock> start, end;
        vb resullts(requests.size());

        start = std::chrono::system_clock::now();
        for (size_t i = 0; i < requests.size(); ++i) {
            auto r = std::lower_bound(helper.begin() + 1, helper.end(), requests[i]);

            auto idx = std::distance(helper.begin(), r);
            unsigned long long le = (idx - 1) * sqrt_n;
            unsigned long long rg = idx * sqrt_n;
            resullts[i] = std::binary_search(data.begin() + le, data.begin() + rg + 1, requests[i]);
        }
        end = std::chrono::system_clock::now();

        ull total = static_cast<ull>(std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count());

        double per_request = total * 1. / requests.size();
        return std::make_pair(std::make_pair(total, per_request), resullts);

    }


    timeres test_veb(const std::vector<unsigned long long>& data,
                const std::vector<unsigned long long>& requests) {
        auto sorted_data = binsearch::veb::make_layout(data);
        auto len = sorted_data.size();

        // Search array
        size_t search_len = requests.size();
        auto prec = binsearch::veb::get_precalc(data.size());

        misc::drop_caches();

        std::chrono::time_point<std::chrono::system_clock> start, end;
        vb resullts(requests.size());

        start = std::chrono::system_clock::now();
        for (size_t i = 0; i < search_len; i++) {
            resullts[i] = (bool)(binsearch::veb::veb_search(sorted_data.begin(), len, requests[i], prec) == -1);
        }
        end = std::chrono::system_clock::now();

        ull total = static_cast<ull>(std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count());

        double per_request = total * 1. / requests.size();
        return std::make_pair(std::make_pair(total, per_request), resullts);
    }

    bool check_results(const std::vector<vb>& results) {
        for (size_t i = 0; i < results.size() - 1; ++i) {
            if (!std::equal(results[i].begin(), results[i].end(), results[i+1].begin())) {
                for (auto j = results[i].begin(), k = results[i+1].begin(); j != results[i].end(), k != results[i + 1].end(); j++, k++) {
                    if (*j != *k) {
                        std::cerr << "Error on: " << i << "\t" << labs(std::distance(j, results[i].begin())) << std::endl;
                    }
                }
                return false;
            }
        }
        return true;
    }

    void trprint(timeres res, const std::string& name) {
        std::cout << "Name: " << name << "\nTotal:\t" << res.first.first << " ns, per request:\t"
                  << std::setprecision(4) << res.first.second << " ns\n" << std::endl;
        std::cout.flush();
    }

    bool run_tests(size_t n, std::vector<unsigned long long> &data,
                   size_t m, std::vector<unsigned long long> &reqs) {
        std::vector<vb> results;

        // default self implemented
        auto basic_self = test_basic_self(data, reqs);
        results.push_back(basic_self.second);
        trprint(basic_self, "Self implemented binsearch");

        // default STL
        auto basic = test_basic(data, reqs);
        results.push_back(basic.second);
        trprint(basic, "STL binsearch");

        // sqrt
        auto bs_sqrt = test_sqrt(data, reqs, n);
        results.push_back(bs_sqrt.second);
        trprint(bs_sqrt, "sqrt-optimized binsearch");

        // prefetch
        auto pref = test_prefetch(data, reqs);
        results.push_back(pref.second);
        trprint(pref, "prefetch-optimized binsearch");

        // veb
        auto veb = test_veb(data, reqs);
        results.push_back(veb.second);
        trprint(veb, "VEB-layout binsearch");

        return check_results(results);

    };
}


int main(int argc, char **argv) {
    unsigned long long n, m;
    std::vector<unsigned long long> data, requests;

    if (!(argc == 2 && std::string(argv[1]) == "--auto")) {
        printf("Starting manual testing mode, provide data size and data, requests count and requests\n\n");

        std::cin >> n;
        data = misc::read_data(n);

        std::cin >> m;
        requests = misc::read_data(m);

        bool res = tests::run_tests(n, data, m, requests);
        if (res) {
            std::cout << "Tests success" << std::endl;
        } else {
            std::cout << "Tests failed" << std::endl;
        }
    } else {
        printf("Starting automated testing mode");

        std::vector<bool> res;
        for (size_t i = 2; i < 9; ++i) {
            n = static_cast<unsigned long long>(pow(10, i));
            m = static_cast<unsigned long long>(pow(10, 6));
            printf("data size: %llu, requests count: %llu\n\n", n, m);

            data = rands::generate_data(n);
            requests = rands::generate_requests(m);

            res.push_back(tests::run_tests(n, data, m, requests));
        }

        if (std::all_of(res.begin(), res.end(), [] (bool r) {return r;})) {
            std::cout << "Tests success" << std::endl;
        } else {
            std::cout << "Tests failed" << std::endl;
        }
    }

    return 0;
}
