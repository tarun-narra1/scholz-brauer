#include <cstdint>
#include <vector>
#include <fstream>
#include <iostream>
#include <limits>
#include <string>
#include <algorithm>
#include <atomic>
#include <mutex>
#include <thread>

using u64 = std::uint64_t;

#define main chain_solver_cli_main
#include "chain_solver.cpp"
#undef main

static std::vector<u64> enumerate_compact_kernels() {
    std::vector<u64> kernels;
    kernels.reserve(4096);
    int gaps[6];

    for (int a = 0; a < 4; ++a) {
        gaps[0] = a;
        for (int b = 0; b < 4; ++b) {
            gaps[1] = b;
            for (int c = 0; c < 4; ++c) {
                gaps[2] = c;
                for (int d = 0; d < 4; ++d) {
                    gaps[3] = d;
                    for (int e = 0; e < 4; ++e) {
                        gaps[4] = e;
                        for (int f = 0; f < 4; ++f) {
                            gaps[5] = f;
                            u64 x = 1;
                            for (int i = 0; i < 6; ++i) x = (x << (gaps[i] + 1)) | 1ULL;
                            kernels.push_back(x);
                        }
                    }
                }
            }
        }
    }

    std::sort(kernels.begin(), kernels.end());
    kernels.erase(std::unique(kernels.begin(), kernels.end()), kernels.end());
    return kernels;
}

static int K_search(u64 m) {
    (void)m;
    return 6;
}

static int min_chain_length_wt7(u64 n) {
    if (n == 1) return 0;
    if (n == 2) return 1;
    if ((n & (n - 1)) == 0) return bit_length(n) - 1;
    if (popcount(n) != 7) return min_chain_length(n);
    int lambda = bit_length(n) - 1;
    int lb = lambda;
    int ub = lambda + 6;
    u64 chain[64];
    chain[0] = 1;
    for (int d = lb; d <= ub; ++d) {
        if (dfs_chain(chain, 1, n, d)) return d;
    }
    return ub;
}

static int compute_s_floor_for_kernel(u64 m, int& k_stop_out) {
    (void)K_search;
    int L0 = min_chain_length_wt7(m);
    int s0 = L0 - (bit_length(m) - 1);

    if (s0 <= 3) {
        k_stop_out = 0;
        return s0;
    }

    if (m > (std::numeric_limits<u64>::max() >> 1)) {
        k_stop_out = 0;
        return s0;
    }
    u64 m2 = m << 1;
    int L1 = min_chain_length_wt7(m2);

    if (L1 == L0) {
        k_stop_out = 1;
        return s0 - 1;
    }

    k_stop_out = 0;
    return s0;
}

int main(int argc, char** argv){
    const std::string out_path = (argc >= 2) ? argv[1] : "trailing_floor_cert.csv";
    auto kernels = enumerate_compact_kernels();
    if (kernels.size() != 4096) {
        std::cerr << "Expected 4096 kernels, got " << kernels.size() << "\n";
        return 2;
    }

    struct Row { u64 m; int s_floor; int k_stop; };
    std::vector<Row> rows(kernels.size());

    std::atomic<size_t> next_idx(0);
    std::atomic<size_t> completed(0);
    std::mutex progress_mu;

    unsigned int T = std::thread::hardware_concurrency();
    if (T == 0) T = 8;
    if (T < 8) T = 8;
    if (kernels.size() < T) T = (unsigned int)kernels.size();

    auto worker = [&]() {
        while (true) {
            size_t idx = next_idx.fetch_add(1);
            if (idx >= kernels.size()) return;
            u64 m = kernels[idx];
            int k_stop = -1;
            int s_floor = compute_s_floor_for_kernel(m, k_stop);

            rows[idx] = Row{m, s_floor, k_stop};

            size_t done = completed.fetch_add(1) + 1;
            if ((done % 128) == 0) {
                std::lock_guard<std::mutex> lock(progress_mu);
                std::cerr << "Progress: " << done << " / " << kernels.size() << "\n";
            }
        }
    };

    std::vector<std::thread> threads;
    threads.reserve(T);
    for (unsigned int i = 0; i < T; ++i) threads.emplace_back(worker);
    for (auto& th : threads) th.join();

    std::sort(rows.begin(), rows.end(), [](const Row& a, const Row& b) { return a.m < b.m; });

    std::ofstream out(out_path);
    if (!out) {
        std::cerr << "Could not open output: " << out_path << "\n";
        return 2;
    }
    out << "m,s_floor,k_stop\n";
    for (const auto& r : rows) out << r.m << "," << r.s_floor << "," << r.k_stop << "\n";
    out.flush();
    out.close();
    std::cerr << "Wrote " << out_path << "\n";
    return 0;
}
