#include <cmath>
#include <iostream>
#include <fstream>
#include <chrono>
#include <utility>
#include <vector>
#include <algorithm>
#include <zconf.h>
#include <fcntl.h>
#include <functional>
#include <iomanip>

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

    bool bf_binary_search(std::vector<unsigned long long int>::const_iterator begin, std::vector<unsigned long long int>::const_iterator end, unsigned long key)
    {
        size_t half, dist = static_cast<size_t>(std::abs(std::distance(begin, end)));
        std::vector<unsigned long long int>::const_iterator mid;
        while((half = (dist >> 1))) {
            mid = begin + half;
            begin = (*mid < key) ? mid : begin;
            dist -= half;
        }

        return ((*(begin+dist)) == key);
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


    bool bf_prefetch_binary_search(std::vector<unsigned long long int>::const_iterator begin, std::vector<unsigned long long int>::const_iterator end, unsigned long key)
    {
        long half, dist = std::distance(begin, end);
        std::vector<unsigned long long int>::const_iterator mid;
        while(dist > 1) {
            half = dist >> 1;
            __builtin_prefetch(&(*(begin + (half >> 1))), 0, 3);
            __builtin_prefetch(&(*(begin + half + (half >> 1))), 0, 3);

            mid = begin + half;
            begin = (*mid < key) ? mid : begin;
            dist -= half;
        }

        return ((*mid) == key) | ((*(mid + 1)) == key);
    }
    
    namespace veb {
        const size_t NULLNODEPTR = UINTMAX_MAX;

        class Node {
        public:
            explicit Node(unsigned long long val = 0,
                          size_t left = NULLNODEPTR,
                          size_t right = NULLNODEPTR)
                    : val(val), left(left), right(right) {}

            unsigned long long val = 0;
            size_t left = NULLNODEPTR;
            size_t right = NULLNODEPTR;
        };

        typedef std::vector<Node> Tree;

        namespace implicit_veb {
            class Helpers {
            public:
                Helpers()
                        : T(0)
                        , B(0)
                        , D(0)
                {}

                explicit Helpers(size_t size)
                        : T(size)
                        , B(size)
                        , D(size)
                {}

                std::vector<size_t> T, B, D;
            };

            class VebLayout {
            private:
                Tree tree;
                Helpers helps;
                std::vector<unsigned long long> layout;
                std::vector<size_t> p;

                void setTree(size_t size) {

                    tree = Tree(size);
                    for (size_t i = 0; i < size; ++i) {
                        size_t iLeft = 2 * i + 1, iRight = 2 * i + 2;
                        tree[i].left = iLeft < size ? iLeft : NULLNODEPTR;
                        tree[i].right = iRight < size ? iRight : NULLNODEPTR;
                    }
                }

                void putDataTree(const std::vector<unsigned long long>& data) {
                    setTree(data.size());
                    putData(0, data.begin());
                }

                std::vector<unsigned long long> splitTree(unsigned long long root, size_t height) {

                    if (root == NULLNODEPTR) return std::vector<unsigned long long>(0);

                    if (!height) {
                        std::vector<unsigned long long> subs(0);
                        if (tree[root].left != NULLNODEPTR) {
                            subs.push_back(reinterpret_cast<unsigned long long int &&>(tree[root].left));
                        }
                        tree[root].left = NULLNODEPTR;

                        if (tree[root].right != NULLNODEPTR) {
                            subs.push_back(reinterpret_cast<unsigned long long int &&>(tree[root].right));
                        }
                        tree[root].right = NULLNODEPTR;
                        return subs;
                    }

                    auto left_subtree = splitTree(tree[root].left, height - 1);
                    auto right_subtree = splitTree(tree[root].right, height - 1);
                    left_subtree.insert(left_subtree.end(), right_subtree.begin(), right_subtree.end());

                    return left_subtree;
                }

                size_t recursiveVebLayout(unsigned long long root, size_t rDepth,
                                          bool norm=false) {
                    auto h = getHeight(root);
                    if (!h) {
                        layout.push_back(reinterpret_cast<unsigned long long int &&>(tree[root].val));
                        return 1;
                    }

                    size_t subtree_s = norm ? h - 1 : h >> 1;
                    size_t depth = rDepth + subtree_s;

                    auto subtrees = splitTree(root, subtree_s);
                    auto nodes_added = recursiveVebLayout(root, rDepth);

                    size_t left_subtree_size = !subtrees.empty()
                                               ? recursiveVebLayout(subtrees[0], depth + 1)
                                               : 0;
                    helps.T[depth] = nodes_added;
                    helps.D[depth] = rDepth;
                    helps.B[depth] = left_subtree_size;

                    nodes_added += left_subtree_size;
                    // layout subtrees
                    for (size_t i = 1; i < subtrees.size(); ++i)
                        nodes_added += recursiveVebLayout(subtrees[i], depth + 1);

                    return nodes_added;
                }
            public:

                explicit VebLayout(const std::vector<unsigned long long>& data)
                    : tree(data.size())
                {
                    makeVebLayout(data);
                }

                inline size_t getPos(size_t bfs, size_t d) {
                    return p[helps.D[d] - 1] + helps.T[d] + (bfs & helps.T[d]) * helps.B[d];
                }

                // implicit search
                bool search(const unsigned long long & el) {
                    size_t i = 1, pos = 0;

                    for (size_t st = 1; i <= layout.size(); st++) {
                        auto e = layout[pos];
                        if (e == el) return true;
                        i = e < el ? 2 * i + 1 : 2 * i;
                        pos = getPos(i, st) - 1;
                        p[st] = pos + 1;
                    }
                    return false;
                }

                void makeVebLayout(const std::vector<unsigned long long>& data) {
                    putDataTree(data);
                    size_t h = getHeight(0);
                    helps = Helpers(h + 2);
                    p = std::vector<size_t>(h + 2);
                    p[0] = 1;
                    recursiveVebLayout(0, 1, true);
                }

                std::vector<unsigned long long int>::const_iterator putData(
                        unsigned long long int node,
                        std::vector<unsigned long long int>::const_iterator data) {
                    if (node == NULLNODEPTR) return data; // end of recursion
                    data = putData(tree[node].left, data);
                    tree[node].val = (unsigned long long) *(data++);
                    return putData(tree[node].right, data);
                }

                size_t getHeight(unsigned long long root) {
                    size_t h = 0;
                    while (tree[root].left != NULLNODEPTR) {
                        h++;
                        root = tree[root].left;
                    }
                    return h;
                }
            };
        }

        namespace explicit_veb {

            class VebLayout {
            private:
                Tree tree;

            public:
                explicit VebLayout(const std::vector<unsigned long long> &data)
                        : tree(data.size())
                {
                    makeVebLayout(data);
                }

                void makeVebLayout(const std::vector<unsigned long long>& data) {
                    std::vector<Node> tmp(data.size());
                    for (size_t i = 0; i < data.size(); ++i) {
                        tmp[i].val = data[i];
                    }

                    putData(tmp.begin(), tmp.end(), 0);
                }

                void putData(std::vector<Node>::const_iterator begin, std::vector<Node>::const_iterator end,
                             unsigned long long start_pos) {
                    auto size = std::distance(begin, end);

                    if (size == 0) {
                        return;
                    }

                    if (size == 1) {
                        tree[start_pos] = *begin;
                        return;
                    }

                    auto height = static_cast<size_t>(floor(log2(size)) + 1);
                    size_t depth = height >> 1;
                    size_t top_height = height - depth;
                    auto top_size = static_cast<size_t>((1 << top_height) - 1);
                    auto bottom_size = static_cast<size_t>((1 << depth) - 1);
                    size_t bottom = start_pos + top_size;
                    size_t bottom_end = start_pos + size;

                    std::vector<Node> top_tree;

                    bool left = true;
                    // create subtrees
                    while (bottom != bottom_end) {
                        size_t curr_bottom_size = std::min(bottom_size, bottom_end - bottom);

                        putData(begin, begin + curr_bottom_size, bottom);
                        begin += curr_bottom_size;

                        if (left) {
                            top_tree.push_back(*begin);
                            top_tree.back().left = bottom;
                            begin++;
                        } else {
                            top_tree.back().right = bottom;
                            if (begin != end) {
                                top_tree.push_back(*begin);
                            }
                            begin++;
                        }
                        left = !left;
                        bottom += curr_bottom_size;
                    }

                    if (std::distance(begin, end) > 0) top_tree.insert(top_tree.end(), begin, end);
                    putData(top_tree.begin(), top_tree.end(), start_pos);
                }


                bool search(unsigned long long key) {
                    size_t i = 0;

                    while (i != NULLNODEPTR) {
                        if (tree[i].val == key) {
                            return true;
                        } else if (tree[i].val < key) {
                            i = tree[i].right;
                        } else {
                            i = tree[i].left;
                        }
                    }

                    return false;
                }
            };
        }
    }

    namespace LUT {

        std::vector<std::pair<size_t, size_t>> Lut;

        void InitLut(const std::vector<unsigned long long> &vals, size_t lutBits)
        {
            const size_t lutSize = 1u << lutBits;
            const size_t shiftBits = 64 - lutBits;
            Lut.resize(lutSize);

            size_t thresh = 0, last = 0;

            for (size_t i = 0; i < vals.size() - 1; i++) {
                const size_t nextThresh = vals[i + 1] >> shiftBits;
                Lut[thresh] = {last, i};

                if (nextThresh > thresh) {
                    last = i + 1;
                    for (size_t j = thresh + 1; j <= nextThresh; j++)
                        Lut[j] = {last, i + 1};
                }

                thresh = nextThresh;
            }

            for (size_t i = thresh; i < Lut.size(); i++) Lut[i] = {last, vals.size() - 1};
        }

        bool LutBinarySearch(const std::vector<unsigned long long>& data,
                             const unsigned long long& key,
                             const size_t& lutBits) {
            auto start = Lut[key >> (64 - lutBits)];
            return std::binary_search(data.begin() + start.first, data.begin() + start.second, key);
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
        std::uniform_int_distribution<unsigned long long> nums(0, static_cast<unsigned long long int>(1E12));
        std::function<unsigned long long()> number_generator = std::bind(nums, generator);
        std::generate(v.begin(), v.end(), number_generator);
        return v;
    }
}

namespace tests {

    typedef unsigned long long ull;
    typedef std::vector<ull> vull;
    typedef std::vector<bool> vb;

    typedef std::pair<std::pair<ull, double>, vb> timeres;

    timeres test_basic_self(const vull &data, const vull &requests) {
        std::chrono::time_point<std::chrono::system_clock> start, end;
        vb results(requests.size());
        start = std::chrono::system_clock::now();
        for (size_t i = 0; i < requests.size(); ++i) {
            results[i] = binsearch::binary_search(data.begin(), data.end(), requests[i]);
        }
        end = std::chrono::system_clock::now();
        ull total = static_cast<ull>(std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count());

        double per_request = total * 1. / requests.size();
        return std::make_pair(std::make_pair(total, per_request), results);
    };

    timeres test_bf_self(const vull &data, const vull &requests) {
        std::chrono::time_point<std::chrono::system_clock> start, end;
        vb results(requests.size());
        start = std::chrono::system_clock::now();
        for (size_t i = 0; i < requests.size(); ++i) {
            results[i] = binsearch::bf_binary_search(data.begin(), data.end(), requests[i]);
        }
        end = std::chrono::system_clock::now();
        ull total = static_cast<ull>(std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count());

        double per_request = total * 1. / requests.size();
        return std::make_pair(std::make_pair(total, per_request), results);
    };

    timeres test_basic(const vull &data, const vull &requests) {
        std::chrono::time_point<std::chrono::system_clock> start, end;
        vb results(requests.size());
        start = std::chrono::system_clock::now();
        for (size_t i = 0; i < requests.size(); ++i) {
            results[i] = std::binary_search(data.begin(), data.end(), requests[i]);
        }
        end = std::chrono::system_clock::now();
        ull total = static_cast<ull>(std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count());

        double per_request = total * 1. / requests.size();
        return std::make_pair(std::make_pair(total, per_request), results);
    };

    timeres test_prefetch(const vull &data, const vull &requests) {
        std::chrono::time_point<std::chrono::system_clock> start, end;
        vb results(requests.size());
        start = std::chrono::system_clock::now();
        for (size_t i = 0; i < requests.size(); ++i) {
            results[i] = binsearch::prefetch_binary_search(data.begin(), data.end(), requests[i]);
        }
        end = std::chrono::system_clock::now();
        ull total = static_cast<ull>(std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count());

        double per_request = total * 1. / requests.size();
        return std::make_pair(std::make_pair(total, per_request), results);
    };

    timeres test_bf_prefetch(const vull &data, const vull &requests) {
        std::chrono::time_point<std::chrono::system_clock> start, end;
        vb results(requests.size());
        start = std::chrono::system_clock::now();
        for (size_t i = 0; i < requests.size(); ++i) {
            results[i] = binsearch::bf_prefetch_binary_search(data.begin(), data.end(), requests[i]);
        }
        end = std::chrono::system_clock::now();
        ull total = static_cast<ull>(std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count());

        double per_request = total * 1. / requests.size();
        return std::make_pair(std::make_pair(total, per_request), results);
    };

    timeres test_sqrt(const vull &data, const vull &requests) {
        auto n = data.size();
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
        vb results(requests.size());

        start = std::chrono::system_clock::now();
        for (size_t i = 0; i < requests.size(); ++i) {
            auto r = std::lower_bound(helper.begin() + 1, helper.end(), requests[i]);

            auto idx = std::distance(helper.begin(), r);
            unsigned long long le = std::max((idx - 1) * sqrt_n, 0ull);
            unsigned long long rg = std::min(idx * sqrt_n + 1, (ull)(data.size()));
            results[i] = std::binary_search(data.begin() + le, data.begin() + rg, requests[i]);
        }
        end = std::chrono::system_clock::now();

        ull total = static_cast<ull>(std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count());

        double per_request = total * 1. / requests.size();
        return std::make_pair(std::make_pair(total, per_request), results);
    }

    timeres test_bf_sqrt(const vull &data, const vull &requests) {
        auto n = data.size();
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
        vb results(requests.size());

        start = std::chrono::system_clock::now();
        for (size_t i = 0; i < requests.size(); ++i) {
            auto r = std::lower_bound(helper.begin() + 1, helper.end(), requests[i]);

            auto idx = std::distance(helper.begin(), r);
            unsigned long long le = std::max((idx - 1) * sqrt_n, 0ull);
            unsigned long long rg = std::min(idx * sqrt_n + 1, (ull)(data.size()));
            results[i] = binsearch::bf_prefetch_binary_search(data.begin() + le, data.begin() + rg, requests[i]);
        }
        end = std::chrono::system_clock::now();

        ull total = static_cast<ull>(std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count());

        double per_request = total * 1. / requests.size();
        return std::make_pair(std::make_pair(total, per_request), results);
    }

    timeres test_veb(const std::vector<unsigned long long>& data,
                const std::vector<unsigned long long>& requests) {
        auto veblay = binsearch::veb::implicit_veb::VebLayout(data);
        size_t search_len = requests.size();

        misc::drop_caches();

        std::chrono::time_point<std::chrono::system_clock> start, end;
        vb results(requests.size());

        start = std::chrono::system_clock::now();
        for (size_t i = 0; i < search_len; i++) {
            results[i] = (bool) veblay.search(requests[i]);
        }
        end = std::chrono::system_clock::now();

        ull total = static_cast<ull>(std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count());

        double per_request = total * 1. / requests.size();
        return std::make_pair(std::make_pair(total, per_request), results);
    }

    timeres test_explicit_veb(const std::vector<unsigned long long>& data,
                     const std::vector<unsigned long long>& requests) {
        auto veblay = binsearch::veb::explicit_veb::VebLayout(data);

        size_t search_len = requests.size();

        misc::drop_caches();

        std::chrono::time_point<std::chrono::system_clock> start, end;
        vb results(requests.size());

        start = std::chrono::system_clock::now();
        for (size_t i = 0; i < search_len; i++) {
            results[i] = (bool) veblay.search(requests[i]);
        }
        end = std::chrono::system_clock::now();

        ull total = static_cast<ull>(std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count());

        double per_request = total * 1. / requests.size();
        return std::make_pair(std::make_pair(total, per_request), results);
    }

    timeres test_lut(const vull &data, const vull &requests) {
        std::chrono::time_point<std::chrono::system_clock> start, end;
        vb results(requests.size());

        binsearch::LUT::InitLut(data, binsearch::LUT::LUT_BITS);

        misc::drop_caches();

        start = std::chrono::system_clock::now();
        for (size_t i = 0; i < requests.size(); ++i) {
            results[i] = binsearch::LUT::LutBinarySearch(data, requests[i], binsearch::LUT::LUT_BITS);
        }
        end = std::chrono::system_clock::now();
        ull total = static_cast<ull>(std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count());

        double per_request = total * 1. / requests.size();
        return std::make_pair(std::make_pair(total, per_request), results);
    };

    bool check_results(const std::vector<vb>& results, bool debug = true) {
        for (size_t i = 1; i < results.size(); ++i) {
            if (!std::equal(results[0].begin(), results[0].end(), results[i].begin())) {
                if (debug) {
                    for (auto j = results[0].begin(), k = results[i].begin();
                         j != results[0].end(), k != results[i].end();
                         j++, k++) {
                        if (*j != *k) {
                            std::cerr << "Error on: " << i << "\t" << labs(std::distance(k, results[i].begin()))
                                      << std::endl;
                        }
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

    void save_res(const vb& data, const std::string& postfix) {
        std::ofstream file;
        file.open("result_" + postfix + ".txt");
        for (auto e: data) {
            file << (e ? 1 : 0) << " ";
        }
        file << std::endl;
        file.close();
    }

    bool run_tests(std::vector<unsigned long long> &data, std::vector<unsigned long long> &reqs) {
//        std::vector<vb> results;

        // veb
        auto veb = test_veb(data, reqs);
//        results.push_back(veb.second);
        save_res(veb.second, "veb");
        trprint(veb, "VEB-layout binsearch");

        // explicit veb
        auto explicit_veb = test_explicit_veb(data, reqs);
//        results.push_back(explicit_veb.second);
        save_res(explicit_veb.second, "explicit_veb");
        trprint(explicit_veb, "Explicit VEB-layout binsearch");

        // default STL
        auto basic = test_basic(data, reqs);
//        results.push_back(basic.second);
        save_res(basic.second, "bs");
        trprint(basic, "STL binsearch");

        // default self implemented
        auto basic_self = test_basic_self(data, reqs);
//        results.push_back(basic_self.second);
        save_res(basic_self.second, "bs_self");
        trprint(basic_self, "Self implemented binsearch");

        // branch free self implemented
        auto bf_self = test_bf_self(data, reqs);
//        results.push_back(bf_self.second);
        save_res(bf_self.second, "bf_self");
        trprint(bf_self, "branch-free binsearch");

        // sqrt
        auto bs_sqrt = test_sqrt(data, reqs);
//        results.push_back(bs_sqrt.second);
        save_res(bs_sqrt.second, "sqrt");
        trprint(bs_sqrt, "sqrt-optimized binsearch");

        // sqrt boosted
        auto bf_sqrt = test_bf_sqrt(data, reqs);
//        results.push_back(bf_sqrt.second);
        save_res(bf_sqrt.second, "sqrt");
        trprint(bf_sqrt, "bf prefetch boosted sqrt-optimized binsearch");

        // prefetch
        auto pref = test_prefetch(data, reqs);
//        results.push_back(pref.second);
        save_res(pref.second, "prefetch");
        trprint(pref, "prefetch-optimized binsearch");

        // bf prefetch
        auto bfpref = test_bf_prefetch(data, reqs);
//        results.push_back(bfpref.second);
        save_res(bfpref.second, "bf_prefetch");
        trprint(bfpref, "branch-free prefetch-optimized binsearch");

        // LUT
        auto lut = test_lut(data, reqs);
//        results.push_back(lut.second);
        save_res(lut.second, "lut");
        trprint(lut, "Look-up table based binsearch");

//        return check_results(results);
        return true;
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

        bool res = tests::run_tests(data, requests);
        if (res) {
            std::cout << "Tests success" << std::endl;
        } else {
            std::cout << "Tests failed" << std::endl;
        }
    } else {
        printf("Starting automated testing mode\n\n");

        std::vector<bool> res;
        const size_t TESTS_NUM = 30;
        for (size_t i = 0; i <= TESTS_NUM; ++i) {
            std::cerr << i << "/" << TESTS_NUM << std::endl;
            n = static_cast<unsigned long long int>(std::pow(10, 4) + i * (std::pow(10, 8) - std::pow(10, 4)) / TESTS_NUM);
            m = static_cast<unsigned long long>(pow(10, 5));
            printf("data size: %llu, requests count: %llu\n\n", n, m);

            data = rands::generate_data(n);
            requests = rands::generate_requests(m);

            res.push_back(tests::run_tests(data, requests));
        }

        if (std::all_of(res.begin(), res.end(), [] (bool r) {return r;})) {
            std::cout << "Tests success" << std::endl;
        } else {
            std::cout << "Tests failed" << std::endl;
        }
    }

    return 0;
}
