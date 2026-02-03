
#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
#include <chrono>
#include <cstdint>
#include <atomic>
#include <mutex>
#include <thread>
#include <bitset>
#include <unordered_set>

#ifdef _OPENMP
#include <omp.h>
#endif

using namespace std;
using u64 = uint64_t;

mutex output_mutex;

inline int bit_length(u64 n) { return n ? 64 - __builtin_clzll(n) : 0; }
inline int popcount(u64 n) { return __builtin_popcountll(n); }





static constexpr int MAX_CHAIN_ELEMS = 256;



static inline u64 next_combination(u64 x) {
    u64 c = x & (~x + 1); 
    u64 r = x + c;
    return (((r ^ x) >> 2) / c) | r;
}

static vector<u64> enumerate_weight7_up_to(u64 N) {
    vector<u64> out;
    
    if (N < 127) return out;
    u64 x = (1ULL << 7) - 1;
    while (x <= N) {
        out.push_back(x);
        u64 nx = next_combination(x);
        if (nx <= x) break; 
        x = nx;
    }
    return out;
}

static inline int trailing_zeros_u64(u64 n) {
    return n ? __builtin_ctzll(n) : 0;
}

static inline int max_internal_gap_u64(u64 n) {
    int prev = -1;
    int best = 0;
    int idx = 0;
    while (n) {
        if (n & 1ULL) {
            if (prev >= 0) best = max(best, idx - prev - 1);
            prev = idx;
        }
        n >>= 1;
        ++idx;
    }
    return best;
}

static inline int max_gap_u64(u64 n) {
    return max(max_internal_gap_u64(n), trailing_zeros_u64(n));
}

static vector<pair<int, int>> compute_addition_certificate(const vector<u64>& chain) {
    vector<pair<int, int>> cert;
    cert.reserve(chain.size());
    if (chain.empty()) return cert;
    cert.push_back({-1, -1}); 

    for (int k = 1; k < (int)chain.size(); ++k) {
        bool ok = false;
        for (int i = k - 1; i >= 0 && !ok; --i) {
            for (int j = i; j >= 0; --j) {
                u64 s = chain[i] + chain[j];
                if (s == chain[k]) {
                    cert.push_back({i, j});
                    ok = true;
                    break;
                }
                if (s < chain[k]) break; 
            }
        }
        if (!ok) {
            cert.push_back({-2, -2}); 
        }
    }
    return cert;
}

static void print_chain_with_certificate(const vector<u64>& chain) {
    auto cert = compute_addition_certificate(chain);
    cout << "len=" << (chain.empty() ? 0 : (int)chain.size() - 1) << " elems=" << chain.size() << "\n";
    for (int k = 0; k < (int)chain.size(); ++k) {
        cout << "  a[" << k << "]=" << chain[k];
        if (k == 0) {
            cout << "\n";
            continue;
        }
        auto [i, j] = cert[k];
        if (i >= 0 && j >= 0) {
            cout << " = a[" << i << "] + a[" << j << "] (" << chain[i] << "+" << chain[j] << ")";
        } else {
            cout << " = [certificate-missing]";
        }
        cout << "\n";
    }
}




bool dfs_star(u64* chain, int chain_len, u64 target, int remaining) {
    u64 last = chain[chain_len - 1];
    
    if (last == target) return true;
    if (remaining <= 0) return false;
    
    
    
    if (last > 0 && remaining > 0) {
        int bl = bit_length(last);
        if (bl + remaining < 64) {
            if ((last << remaining) < target) return false;
        }
        
        
    }
    
    if (last > target) return false;
    
    
    for (int j = chain_len - 1; j >= 0; --j) {
        u64 next = last + chain[j];
        if (next <= last) continue;
        if (next > target) continue;
        
        chain[chain_len] = next;
        if (dfs_star(chain, chain_len + 1, target, remaining - 1)) {
            return true;
        }
    }
    
    return false;
}

int min_star_chain_length(u64 n) {
    if (n == 1) return 0;
    if (n == 2) return 1;
    if ((n & (n - 1)) == 0) return bit_length(n) - 1;  
    
    int lambda = bit_length(n) - 1;
    int weight = popcount(n);
    
    
    
    int lb = lambda + 1;
    
    int ub = lambda + weight - 1;
    
    u64 chain[MAX_CHAIN_ELEMS];
    chain[0] = 1;
    
    for (int d = lb; d <= ub; ++d) {
        if (dfs_star(chain, 1, n, d)) {
            return d;
        }
    }
    
    return ub;
}

int min_star_chain_length_with_chain(u64 n, vector<u64>& out_chain) {
    out_chain.clear();
    if (n == 1) {
        out_chain.push_back(1);
        return 0;
    }
    if (n == 2) {
        out_chain.push_back(1);
        out_chain.push_back(2);
        return 1;
    }
    if ((n & (n - 1)) == 0) {
        int r = bit_length(n) - 1;
        out_chain.reserve(r + 1);
        u64 x = 1;
        out_chain.push_back(x);
        for (int i = 0; i < r; ++i) {
            x *= 2;
            out_chain.push_back(x);
        }
        return r;
    }

    int lambda = bit_length(n) - 1;
    int weight = popcount(n);
    int lb = lambda + 1;
    int ub = lambda + weight - 1;

    u64 chain[MAX_CHAIN_ELEMS];
    chain[0] = 1;

    for (int d = lb; d <= ub; ++d) {
        if (dfs_star(chain, 1, n, d)) {
            out_chain.assign(chain, chain + (d + 1));
            return d;
        }
    }

    out_chain.assign(chain, chain + (ub + 1));
    return ub;
}


bool dfs_chain(u64* chain, int chain_len, u64 target, int remaining) {
    u64 last = chain[chain_len - 1];
    
    if (last == target) return true;
    if (remaining <= 0) return false;
    
    
    
    if (last > 0 && remaining > 0) {
        int bl = bit_length(last);
        if (bl + remaining < 64) {
            if ((last << remaining) < target) return false;
        }
    }
    
    if (last > target) return false;
    
    
    
    vector<u64> cands;
    cands.reserve(chain_len * chain_len);  
    
    for (int i = chain_len - 1; i >= 0; --i) {
        u64 ai = chain[i];
        for (int j = i; j >= 0; --j) {
            u64 s = ai + chain[j];
            if (s <= last) break;  
            if (s <= target) {
                cands.push_back(s);
            }
        }
    }
    
    
    sort(cands.begin(), cands.end(), greater<u64>());
    
    
    cands.erase(unique(cands.begin(), cands.end()), cands.end());
    
    for (u64 val : cands) {
        chain[chain_len] = val;
        if (dfs_chain(chain, chain_len + 1, target, remaining - 1)) {
            return true;
        }
    }
    
    return false;
}

int min_chain_length(u64 n, int upper_bound_hint = -1) {
    if (n == 1) return 0;
    if (n == 2) return 1;
    if ((n & (n - 1)) == 0) return bit_length(n) - 1;
    
    int lambda = bit_length(n) - 1;
    int weight = popcount(n);
    
    
    
    int lb = lambda + 1;
    int ub = (upper_bound_hint > 0) ? upper_bound_hint : (lambda + weight - 1);
    
    u64 chain[MAX_CHAIN_ELEMS];
    chain[0] = 1;
    
    for (int d = lb; d <= ub; ++d) {
        if (dfs_chain(chain, 1, n, d)) {
            return d;
        }
    }
    
    return ub;
}




static bool has_chain_leq(u64 n, int max_len) {
    if (n == 1) return (0 <= max_len);
    if (n == 2) return (1 <= max_len);
    if ((n & (n - 1)) == 0) return (bit_length(n) - 1 <= max_len);
    if (max_len < 0) return false;

    int lambda = bit_length(n) - 1;
    int lb = lambda + 1;  
    if (max_len < lb) return false;

    u64 chain[MAX_CHAIN_ELEMS];
    chain[0] = 1;
    for (int d = lb; d <= max_len; ++d) {
        if (dfs_chain(chain, 1, n, d)) return true;
    }
    return false;
}

int min_chain_length_with_chain(u64 n, vector<u64>& out_chain, int upper_bound_hint = -1) {
    out_chain.clear();
    if (n == 1) {
        out_chain.push_back(1);
        return 0;
    }
    if (n == 2) {
        out_chain.push_back(1);
        out_chain.push_back(2);
        return 1;
    }
    if ((n & (n - 1)) == 0) {
        int r = bit_length(n) - 1;
        out_chain.reserve(r + 1);
        u64 x = 1;
        out_chain.push_back(x);
        for (int i = 0; i < r; ++i) {
            x *= 2;
            out_chain.push_back(x);
        }
        return r;
    }

    int lambda = bit_length(n) - 1;
    int weight = popcount(n);
    int lb = lambda + 1;
    int ub = (upper_bound_hint > 0) ? upper_bound_hint : (lambda + weight - 1);

    u64 chain[MAX_CHAIN_ELEMS];
    chain[0] = 1;

    for (int d = lb; d <= ub; ++d) {
        if (dfs_chain(chain, 1, n, d)) {
            out_chain.assign(chain, chain + (d + 1));
            return d;
        }
    }

    out_chain.assign(chain, chain + (ub + 1));
    return ub;
}


pair<int, int> compute_chain_lengths(u64 n) {
    
    int ls = min_star_chain_length(n);
    
    
    
    int l = min_chain_length(n, ls);
    
    return {l, ls};
}

pair<int, int> compute_chain_lengths_with_witness(u64 n, vector<u64>& chain_general, vector<u64>& chain_star) {
    int ls = min_star_chain_length_with_chain(n, chain_star);
    int l = min_chain_length_with_chain(n, chain_general, ls);
    return {l, ls};
}



static inline u64 odd_part(u64 n) {
    return n ? (n >> __builtin_ctzll(n)) : 0;
}

static inline bool defect_lt_1(u64 n, int l) {
    
    
    __int128 left = (__int128)1 << l;
    __int128 right = (__int128)2 * (__int128)n;
    return left < right;
}

static inline bool defect_lt_d4(u64 n, int l) {
    
    
    __int128 left = (__int128)13 << l;
    __int128 right = (__int128)32 * (__int128)n;
    return left < right;
}

static inline bool defect_lt_d5(u64 n, int l) {
    
    
    __int128 left = (__int128)21 << l;
    __int128 right = (__int128)64 * (__int128)n;
    return left < right;
}

static inline bool is_power_of_two_u64(u64 x) {
    return x && ((x & (x - 1)) == 0);
}

static inline bool in_d4_classification(u64 n) {
    if (popcount(n) <= 2) return true;
    u64 m = odd_part(n);
    return (m == 7 || m == 15 || m == 27);
}

static inline bool in_d5_classification(u64 n) {
    if (popcount(n) <= 2) return true;
    u64 m = odd_part(n);
    if (m == 11 || m == 23 || m == 43 || m == 45 || m == 85) return true;

    
    if (m > 1 && (m - 1) % 3 == 0) {
        u64 t = (m - 1) / 3;
        if (t >= 2 && is_power_of_two_u64(t)) return true;
    }

    
    if (m % 3 == 0) {
        u64 t = m / 3;
        if (t >= 3 && is_power_of_two_u64(t - 1)) return true;
    }

    return false;
}

static inline bool in_lt1_classification(u64 n) {
    return popcount(n) <= 2;
}

static inline string u128_to_string(__uint128_t x) {
    if (x == 0) return "0";
    string s;
    while (x > 0) {
        int digit = (int)(x % 10);
        s.push_back((char)('0' + digit));
        x /= 10;
    }
    reverse(s.begin(), s.end());
    return s;
}

static inline __uint128_t fib_u128(int n) {
    
    if (n <= 0) return 0;
    if (n == 1) return 1;
    __uint128_t a = 0, b = 1;
    for (int i = 2; i <= n; ++i) {
        __uint128_t c = a + b;
        a = b;
        b = c;
    }
    return b;
}




struct Outlier {
    u64 n;
    int gaps[7];
    int l, ls;
    bool is_brauer;
    double time_sec;
};

vector<Outlier> generate_outliers() {
    const u64 N0 = 5784688;
    vector<Outlier> result;
    
    for (int g1 = 0; g1 < 4; g1++)
    for (int g2 = 0; g2 < 4; g2++)
    for (int g3 = 0; g3 < 4; g3++)
    for (int g4 = 0; g4 < 4; g4++)
    for (int g5 = 0; g5 < 4; g5++)
    for (int g6 = 0; g6 < 4; g6++)
    for (int t = 0; t < 4; t++) {
        u64 n = 1;
        n = (n << (g1 + 1)) | 1;
        n = (n << (g2 + 1)) | 1;
        n = (n << (g3 + 1)) | 1;
        n = (n << (g4 + 1)) | 1;
        n = (n << (g5 + 1)) | 1;
        n = (n << (g6 + 1)) | 1;
        n = n << t;
        
        if (n > N0) {
            Outlier o;
            o.n = n;
            o.gaps[0] = g1; o.gaps[1] = g2; o.gaps[2] = g3;
            o.gaps[3] = g4; o.gaps[4] = g5; o.gaps[5] = g6;
            o.gaps[6] = t;
            o.l = o.ls = 0;
            o.is_brauer = false;
            o.time_sec = 0;
            result.push_back(o);
        }
    }
    
    sort(result.begin(), result.end(), [](const Outlier& a, const Outlier& b) {
        return a.n < b.n;
    });
    
    return result;
}

static vector<Outlier> generate_compact_domain() {
    vector<Outlier> result;
    result.reserve(4 * 4 * 4 * 4 * 4 * 4 * 4);

    for (int g1 = 0; g1 < 4; g1++)
    for (int g2 = 0; g2 < 4; g2++)
    for (int g3 = 0; g3 < 4; g3++)
    for (int g4 = 0; g4 < 4; g4++)
    for (int g5 = 0; g5 < 4; g5++)
    for (int g6 = 0; g6 < 4; g6++)
    for (int t = 0; t < 4; t++) {
        u64 n = 1;
        n = (n << (g1 + 1)) | 1;
        n = (n << (g2 + 1)) | 1;
        n = (n << (g3 + 1)) | 1;
        n = (n << (g4 + 1)) | 1;
        n = (n << (g5 + 1)) | 1;
        n = (n << (g6 + 1)) | 1;
        n = n << t;

        Outlier o;
        o.n = n;
        o.gaps[0] = g1; o.gaps[1] = g2; o.gaps[2] = g3;
        o.gaps[3] = g4; o.gaps[4] = g5; o.gaps[5] = g6;
        o.gaps[6] = t;
        o.l = o.ls = 0;
        o.is_brauer = false;
        o.time_sec = 0;
        result.push_back(o);
    }

    sort(result.begin(), result.end(), [](const Outlier& a, const Outlier& b) {
        return a.n < b.n;
    });

    return result;
}

static void run_compact_domain_verification() {
    cout << "================================================================" << endl;
    cout << "SCHOLZ WEIGHT-7: COMPACT DOMAIN VERIFICATION" << endl;
    cout << "Parallel C++ with Star-First Optimization" << endl;
    cout << "Domain: all gap vectors (g1..g6,t) in {0,1,2,3}^7 (|D_C|=4^7=16384)" << endl;
    cout << "================================================================" << endl << endl;

    auto items = generate_compact_domain();

    cout << "Items: " << items.size() << endl;
    cout << "Range: [" << items.front().n << ", " << items.back().n << "]" << endl;

    #ifdef _OPENMP
    int num_threads = omp_get_max_threads();
    cout << "Threads: " << num_threads << endl;
    #else
    cout << "Threads: 1 (OpenMP not enabled)" << endl;
    #endif

    cout << endl << "Starting computation..." << endl << endl;

    atomic<int> completed(0);
    atomic<int> ex0(0), ex1(0), ex2p(0);
    auto total_start = chrono::high_resolution_clock::now();

    #pragma omp parallel for schedule(dynamic)
    for (size_t i = 0; i < items.size(); ++i) {
        auto start = chrono::high_resolution_clock::now();
        auto [l, ls] = compute_chain_lengths(items[i].n);
        auto end = chrono::high_resolution_clock::now();
        double elapsed = chrono::duration<double>(end - start).count();

        items[i].l = l;
        items[i].ls = ls;
        items[i].is_brauer = (l == ls);
        items[i].time_sec = elapsed;

        int ex = ls - l;
        if (ex == 0) ++ex0;
        else if (ex == 1) ++ex1;
        else ++ex2p;

        int done = ++completed;
        if (done % 256 == 0) {
            lock_guard<mutex> lock(output_mutex);
            cout << "Progress: " << done << " / " << items.size() << "\n";
        }
    }

    auto total_end = chrono::high_resolution_clock::now();
    double total_time = chrono::duration<double>(total_end - total_start).count();

    sort(items.begin(), items.end(), [](const Outlier& a, const Outlier& b) {
        return a.n < b.n;
    });

    {
        ofstream f("compact_brauer_results.csv");
        f << "n,l,l_star,ex,is_brauer,g1,g2,g3,g4,g5,g6,t,time_sec" << endl;
        for (const auto& o : items) {
            f << o.n << "," << o.l << "," << o.ls << "," << (o.ls - o.l) << ","
              << (o.is_brauer ? "True" : "False");
            for (int i = 0; i < 7; ++i) f << "," << o.gaps[i];
            f << "," << o.time_sec << endl;
        }
    }

    cout << endl;
    cout << "================================================================" << endl;
    cout << "SUMMARY" << endl;
    cout << "================================================================" << endl;
    cout << "Total time: " << (total_time / 60.0) << " minutes (" << total_time << "s)" << endl;
    cout << "Star-excess counts: ex=0: " << ex0 << "  ex=1: " << ex1 << "  ex>=2: " << ex2p << endl;
    cout << endl;
    cout << "Results saved to: compact_brauer_results.csv" << endl;
}

void run_verification() {
    cout << "================================================================" << endl;
    cout << "SCHOLZ WEIGHT-7: BRAUER VERIFICATION" << endl;
    cout << "Parallel C++ with Star-First Optimization" << endl;
    cout << "================================================================" << endl << endl;
    
    auto outliers = generate_outliers();
    
    cout << "Outliers: " << outliers.size() << endl;
    cout << "Range: [" << outliers.front().n << ", " << outliers.back().n << "]" << endl;
    
    #ifdef _OPENMP
    int num_threads = omp_get_max_threads();
    cout << "Threads: " << num_threads << endl;
    #else
    cout << "Threads: 1 (OpenMP not enabled)" << endl;
    #endif
    
    cout << endl << "Starting computation..." << endl << endl;
    
    atomic<int> completed(0);
    atomic<int> brauer_count(0);
    
    auto total_start = chrono::high_resolution_clock::now();
    
    #pragma omp parallel for schedule(dynamic)
    for (size_t i = 0; i < outliers.size(); ++i) {
        auto start = chrono::high_resolution_clock::now();
        
        auto [l, ls] = compute_chain_lengths(outliers[i].n);
        
        auto end = chrono::high_resolution_clock::now();
        double elapsed = chrono::duration<double>(end - start).count();
        
        outliers[i].l = l;
        outliers[i].ls = ls;
        outliers[i].is_brauer = (l == ls);
        outliers[i].time_sec = elapsed;
        
        if (outliers[i].is_brauer) {
            brauer_count++;
        }
        
        int done = ++completed;
        
        
        {
            lock_guard<mutex> lock(output_mutex);
            printf("[%3d/345] n=%12llu  l=%2d  l*=%2d  %s  (%.2fs)\n",
                   done, (unsigned long long)outliers[i].n,
                   l, ls,
                   outliers[i].is_brauer ? "BRAUER" : "FAIL",
                   elapsed);
            
            if (!outliers[i].is_brauer) {
                cout << "*** COUNTEREXAMPLE FOUND: n=" << outliers[i].n << " ***" << endl;
            }
        }
    }
    
    auto total_end = chrono::high_resolution_clock::now();
    double total_time = chrono::duration<double>(total_end - total_start).count();
    
    
    sort(outliers.begin(), outliers.end(), [](const Outlier& a, const Outlier& b) {
        return a.n < b.n;
    });
    
    {
        ofstream f("brauer_results.csv");
        f << "n,l,l_star,is_brauer,g1,g2,g3,g4,g5,g6,t,time_sec" << endl;
        for (const auto& o : outliers) {
            f << o.n << "," << o.l << "," << o.ls << ","
              << (o.is_brauer ? "True" : "False");
            for (int i = 0; i < 7; ++i) f << "," << o.gaps[i];
            f << "," << o.time_sec << endl;
        }
    }
    
    
    cout << endl;
    cout << "================================================================" << endl;
    cout << "SUMMARY" << endl;
    cout << "================================================================" << endl;
    cout << "Total time: " << (total_time / 60.0) << " minutes (" << total_time << "s)" << endl;
    cout << "Brauer: " << brauer_count << " / " << outliers.size() << endl;
    
    if (brauer_count == (int)outliers.size()) {
        cout << endl;
        cout << "*** ALL 345 COMPACT OUTLIERS ARE BRAUER NUMBERS ***" << endl;
        cout << "Hypothesis 3.4 is now PROVEN." << endl;
        cout << "Scholz holds for the entire compact domain D_C." << endl;
    } else {
        cout << endl;
        cout << "*** SOME OUTLIERS ARE NOT BRAUER ***" << endl;
        cout << "Must investigate Hansen numbers or repair the theory." << endl;
    }
    
    cout << endl;
    cout << "Results saved to: brauer_results.csv" << endl;
}

static void check_fib_envelope(int Rmax) {
    cout << "Checking Fibonacci--doubling envelope via exhaustive pattern enumeration..." << endl;
    cout << "Range: r <= " << Rmax << endl;

    auto eval_pattern = [](unsigned mask, int steps, int& t_out) -> __uint128_t {
        
        __uint128_t prev = 1;
        __uint128_t curr = 2;
        int t = 0;
        for (int i = 0; i < steps; ++i) {
            bool is_S = ((mask >> i) & 1u) != 0;
            if (is_S) {
                __uint128_t next = curr + prev;
                prev = curr;
                curr = next;
                ++t;
            } else {
                prev = curr;
                curr = curr * 2;
            }
        }
        t_out = t;
        return curr;
    };

    for (int r = 1; r <= Rmax; ++r) {
        int steps = r - 1;
        vector<__uint128_t> best(r, 0); 

        unsigned total = (steps >= 31) ? 0u : (1u << steps);
        if (steps >= 31) {
            cerr << "Rmax too large for exhaustive enumeration with 32-bit masks (r=" << r << ")" << endl;
            exit(2);
        }

        for (unsigned mask = 0; mask < total; ++mask) {
            int t = 0;
            __uint128_t terminal = eval_pattern(mask, steps, t);
            if (terminal > best[t]) best[t] = terminal;
        }

        for (int t = 0; t <= r - 1; ++t) {
            __uint128_t expected = fib_u128(t + 3) * ((__uint128_t)1 << (r - t - 1));
            if (best[t] != expected) {
                cout << "MISMATCH at r=" << r << " t=" << t << "\n";
                cout << "best=" << u128_to_string(best[t]) << " expected=" << u128_to_string(expected) << "\n";
                exit(1);
            }
        }
    }

    cout << "OK: envelope verified for all r <= " << Rmax << " (all patterns enumerated)\n";
}

static void check_7_15_27_rays(int Bmax) {
    cout << "Checking ray lengths for 7*2^b, 15*2^b, 27*2^b..." << endl;
    cout << "Range: b <= " << Bmax << endl;

    auto check_one = [&](u64 odd, int base_len) {
        for (int b = 0; b <= Bmax; ++b) {
            if (b >= 63) break;
            
            if (bit_length(odd) + b >= 64) break;
            u64 n = odd << b;
            int l = min_chain_length(n);
            int expected = base_len + b;
            if (l != expected) {
                cout << "MISMATCH for n=" << n << " (odd=" << odd << ", b=" << b << ")\n";
                cout << "l(n)=" << l << " expected=" << expected << "\n";
                exit(1);
            }
        }
    };

    check_one(7, 4);
    check_one(15, 5);
    check_one(27, 6);

    cout << "OK: ray lengths verified for b <= " << Bmax << "\n";
}

static void check_d5_rays(int Bmax) {
    cout << "Checking ray lengths for 11*2^b, 23*2^b, 43*2^b, 45*2^b, 85*2^b..." << endl;
    cout << "Range: b <= " << Bmax << endl;

    auto check_one = [&](u64 odd, int base_len) {
        for (int b = 0; b <= Bmax; ++b) {
            if (b >= 63) break;
            if (bit_length(odd) + b >= 64) break;
            u64 n = odd << b;
            int l = min_chain_length(n);
            int expected = base_len + b;
            if (l != expected) {
                cout << "MISMATCH for n=" << n << " (odd=" << odd << ", b=" << b << ")\n";
                cout << "l(n)=" << l << " expected=" << expected << "\n";
                exit(1);
            }
        }
    };

    check_one(11, 5);
    check_one(23, 6);
    check_one(43, 7);
    check_one(45, 7);
    check_one(85, 8);

    cout << "OK: d5 ray lengths verified for b <= " << Bmax << "\n";
}

static int scan_weight7(u64 N) {
    cout << "Scanning n <= " << N << " over popcount(n)=7 for star-excess > 1..." << endl;
    cout << "Computes exact l(n), l*(n) and reports any counterexample with l*(n)-l(n) >= 2." << endl;

    struct Counterexample {
        u64 n = 0;
        int l = -1;
        int ls = -1;
        int lambda = -1;
    };

    atomic<bool> found(false);
    mutex ce_mutex;
    Counterexample ce;

    auto candidates = enumerate_weight7_up_to(N);
    cout << "Candidates (wt=7): " << candidates.size() << endl;

    
    atomic<u64> progress_total(0);
    atomic<bool> progress_done(false);

    
    mutex stats_mutex;
    u64 total_checked = 0;
    u64 wt7_count = 0;
    u64 ex0 = 0, ex1 = 0, ex2 = 0, ex3p = 0;
    u64 s3 = 0, s4 = 0, s5 = 0, s6p = 0;
    u64 max_ex = 0;
    u64 max_ex_n = 0;

    unsigned int T = std::thread::hardware_concurrency();
    if (T == 0) T = 4;

    thread reporter([&]() {
        using namespace std::chrono_literals;
        while (!progress_done.load(memory_order_relaxed)) {
            this_thread::sleep_for(10s);
            u64 p = progress_total.load(memory_order_relaxed);
            u64 tot = candidates.size();
            if (p > tot) p = tot;
            cerr << "[scan-wt7] progress: ~" << p << "/" << tot << " (" << (tot ? (100.0 * (double)p / (double)tot) : 100.0) << "%)\n";
            cerr.flush();
        }
    });

    constexpr u64 PROGRESS_BATCH = 256;

    auto worker = [&](u64 start, u64 step) {
        u64 local_total = 0;
        u64 local_wt7 = 0;
        u64 local_ex0 = 0, local_ex1 = 0, local_ex2 = 0, local_ex3p = 0;
        u64 local_s3 = 0, local_s4 = 0, local_s5 = 0, local_s6p = 0;
        u64 local_max_ex = 0;
        u64 local_max_ex_n = 0;

        for (u64 idx = start; idx < candidates.size() && !found.load(memory_order_relaxed); idx += step) {
            u64 n = candidates[idx];
            local_total++;
            if ((local_total % PROGRESS_BATCH) == 0) {
                progress_total.fetch_add(PROGRESS_BATCH, memory_order_relaxed);
            }
            local_wt7++;

            auto [l, ls] = compute_chain_lengths(n);
            int lambda = bit_length(n) - 1;
            int ex = ls - l;

            if (local_wt7 == 1) {
                local_max_ex = (u64)ex;
                local_max_ex_n = n;
            }

            if ((u64)ex > local_max_ex) {
                local_max_ex = (u64)ex;
                local_max_ex_n = n;
            }

            if (ex == 0) local_ex0++;
            else if (ex == 1) local_ex1++;
            else if (ex == 2) local_ex2++;
            else if (ex >= 3) local_ex3p++;

            int s = l - lambda;
            if (s == 3) local_s3++;
            else if (s == 4) local_s4++;
            else if (s == 5) local_s5++;
            else if (s >= 6) local_s6p++;

            if (ex >= 2) {
                lock_guard<mutex> lock(ce_mutex);
                if (!found.load(memory_order_relaxed)) {
                    ce = Counterexample{n, l, ls, lambda};
                    found.store(true, memory_order_relaxed);
                }
                return;
            }
        }

        progress_total.fetch_add(local_total % PROGRESS_BATCH, memory_order_relaxed);

        lock_guard<mutex> lock(stats_mutex);
        total_checked += local_total;
        wt7_count += local_wt7;
        ex0 += local_ex0;
        ex1 += local_ex1;
        ex2 += local_ex2;
        ex3p += local_ex3p;
        s3 += local_s3;
        s4 += local_s4;
        s5 += local_s5;
        s6p += local_s6p;
        if (local_max_ex > max_ex || (local_max_ex == max_ex && max_ex_n == 0 && local_max_ex_n != 0)) {
            max_ex = local_max_ex;
            max_ex_n = local_max_ex_n;
        }
    };

    vector<std::thread> threads;
    threads.reserve(T);
    for (u64 i = 0; i < T; ++i) {
        threads.emplace_back(worker, i, (u64)T);
    }
    for (auto& th : threads) th.join();
    progress_done.store(true, memory_order_relaxed);
    reporter.join();

    if (found.load(memory_order_relaxed)) {
        cout << "COUNTEREXAMPLE FOUND at n=" << ce.n << "\n";
        cout << "bits=" << bit_length(ce.n) << " wt=" << popcount(ce.n) << "\n";
        cout << "lambda=" << ce.lambda << " l(n)=" << ce.l << " l*(n)=" << ce.ls << " ex=" << (ce.ls - ce.l) << "\n";
        cout << "bin=" << std::bitset<64>(ce.n).to_string().substr(64 - bit_length(ce.n)) << "\n";
        return 1;
    }

    cout << "OK: no star-excess >= 2 found for popcount(n)=7 with n <= " << N << "\n";
    cout << "Stats: examined " << total_checked << " wt=7 candidates; weight7 count=" << wt7_count << "\n";
    cout << "Star-excess counts: ex=0: " << ex0 << "  ex=1: " << ex1 << "  ex=2: " << ex2 << "  ex>=3: " << ex3p << "\n";
    cout << "Small-steps s=l-lambda counts: s=3: " << s3 << "  s=4: " << s4 << "  s=5: " << s5 << "  s>=6: " << s6p << "\n";
    cout << "Max observed star-excess: " << max_ex << " at n=" << max_ex_n << "\n";
    return 0;
}

static int verify_wt7_gapbound(int gmax, bool require_sparse) {
    cout << "Enumerating all weight-7 gap vectors (g1..g6,t) with entries in [0," << gmax << "] ..." << endl;
    if (require_sparse) {
        cout << "Restricting to sparse cases with max(g1..g6,t) >= 4." << endl;
    }
    cout << "Computes exact l*(n) and checks if any general chain exists of length <= l*(n)-2." << endl;
    cout << "This is sufficient to detect star-excess >= 2 without computing exact l(n) in non-counterexample cases." << endl;

    struct Counterexample {
        u64 n = 0;
        array<int, 7> gaps{};
        int l = -1;
        int ls = -1;
        int lambda = -1;
        int gamma = -1;
    };

    atomic<bool> found(false);
    mutex ce_mutex;
    Counterexample ce;

    mutex stats_mutex;
    u64 total_vectors = 0;
    u64 matched = 0;

    u64 total_possible = 1;
    for (int i = 0; i < 7; ++i) total_possible *= (u64)(gmax + 1);
    cout << "Total gap vectors: " << total_possible << endl;

    atomic<u64> progress_total(0);
    atomic<bool> progress_done(false);
    thread reporter([&]() {
        using namespace std::chrono_literals;
        while (!progress_done.load(memory_order_relaxed)) {
            this_thread::sleep_for(5s);
            u64 p = progress_total.load(memory_order_relaxed);
            if (p > total_possible) p = total_possible;
            cerr << "[verify-wt7-gapbound] progress: ~" << p << "/" << total_possible << " ("
                 << (100.0 * (double)p / (double)total_possible) << "%)\n";
            cerr.flush();
        }
    });

    unsigned int T = std::thread::hardware_concurrency();
    if (T == 0) T = 4;

    constexpr u64 PROGRESS_BATCH = 8;

    const u64 base = (u64)(gmax + 1);

    auto worker = [&](u64 start, u64 step) {
        u64 local_total = 0;
        u64 local_matched = 0;

        for (u64 idx = start; idx < total_possible && !found.load(memory_order_relaxed); idx += step) {
            ++local_total;
            if ((local_total % PROGRESS_BATCH) == 0) {
                progress_total.fetch_add(PROGRESS_BATCH, memory_order_relaxed);
            }

            u64 tmp = idx;
            int g1 = (int)(tmp % base);
            tmp /= base;
            int g2 = (int)(tmp % base);
            tmp /= base;
            int g3 = (int)(tmp % base);
            tmp /= base;
            int g4 = (int)(tmp % base);
            tmp /= base;
            int g5 = (int)(tmp % base);
            tmp /= base;
            int g6 = (int)(tmp % base);
            tmp /= base;
            int t = (int)(tmp % base);

            int gamma = max({g1, g2, g3, g4, g5, g6, t});
            if (require_sparse && gamma < 4) continue;
            ++local_matched;

            
            u64 x = 1;
            x = (x << (g1 + 1)) + 1;
            x = (x << (g2 + 1)) + 1;
            x = (x << (g3 + 1)) + 1;
            x = (x << (g4 + 1)) + 1;
            x = (x << (g5 + 1)) + 1;
            x = (x << (g6 + 1)) + 1;
            u64 n = x << t;

            int lambda = bit_length(n) - 1;
            int ls = min_star_chain_length(n);
            int target = ls - 2;
            
            
            if (target < lambda + 3) continue;
            
            if (target >= 0 && has_chain_leq(n, target)) {
                lock_guard<mutex> lock(ce_mutex);
                if (!found.load(memory_order_relaxed)) {
                    
                    auto [l_exact, ls_exact] = compute_chain_lengths(n);
                    ce = Counterexample{n, {g1, g2, g3, g4, g5, g6, t}, l_exact, ls_exact, lambda, gamma};
                    found.store(true, memory_order_relaxed);
                }
                return;
            }
        }

        progress_total.fetch_add(local_total % PROGRESS_BATCH, memory_order_relaxed);

        lock_guard<mutex> lock(stats_mutex);
        total_vectors += local_total;
        matched += local_matched;
    };

    vector<std::thread> threads;
    threads.reserve(T);
    for (unsigned int i = 0; i < T; ++i) {
        threads.emplace_back(worker, (u64)i, (u64)T);
    }
    for (auto& th : threads) th.join();
    progress_done.store(true, memory_order_relaxed);
    reporter.join();

    if (found.load(memory_order_relaxed)) {
        cout << "COUNTEREXAMPLE FOUND at n=" << ce.n << "\n";
        cout << "gaps(g1..g6,t)=(" << ce.gaps[0] << "," << ce.gaps[1] << "," << ce.gaps[2] << "," << ce.gaps[3] << ","
             << ce.gaps[4] << "," << ce.gaps[5] << "," << ce.gaps[6] << ")\n";
        cout << "gamma=" << ce.gamma << " bits=" << bit_length(ce.n) << " wt=" << popcount(ce.n) << "\n";
        cout << "lambda=" << ce.lambda << " l(n)=" << ce.l << " l*(n)=" << ce.ls << " ex=" << (ce.ls - ce.l) << "\n";
        cout << "bin=" << std::bitset<64>(ce.n).to_string().substr(64 - bit_length(ce.n)) << "\n";
        return 1;
    }

    cout << "OK: no star-excess >= 2 found among enumerated gap vectors with gmax=" << gmax;
    if (require_sparse) cout << " (sparse subset gamma>=4)";
    cout << "\n";
    cout << "Stats: enumerated " << total_vectors << " vectors; matched=" << matched << "\n";
    return 0;
}


static inline int v2_u64(u64 n) { return n ? __builtin_ctzll(n) : 64; }

struct InterfaceKey {
    int r = 0;
    u64 B = 0;
};

static inline bool operator==(const InterfaceKey& a, const InterfaceKey& b) {
    return a.r == b.r && a.B == b.B;
}

struct InterfaceKeyHash {
    size_t operator()(const InterfaceKey& k) const noexcept {
        
        return std::hash<u64>{}((static_cast<u64>(k.r) << 56) ^ k.B);
    }
};




static bool compute_first_gap_interface_shifted(u64 n, int& out_r, u64& out_B) {
    if (n == 0) return false;
    int t = v2_u64(n);
    u64 m = n >> t; 
    if (m == 0) return false;

    
    int e[64];
    int cnt = 0;
    for (int i = 0; i < 64; ++i) {
        if ((m >> i) & 1ULL) e[cnt++] = i;
    }
    if (cnt < 2) return false;
    
    int split_index = -1;
    for (int i = 1; i < cnt; ++i) {
        int gap = e[i] - e[i - 1] - 1;
        if (gap >= 4) {
            split_index = i;
            break;
        }
    }
    if (split_index < 0) return false; 
    int r = e[split_index - 1] + 1; 
    u64 B = m & ((r == 64) ? ~0ULL : ((1ULL << r) - 1ULL));
    if (B == 0) return false;
    out_r = r + t;
    out_B = B << t;
    return true;
}








static int certify_interface_from_gapbound(int gmax, int max_trailing_shift, const string& in_csv, const string& out_csv) {
    
    struct Row {
        int r = 0;
        int R = 0;
        u64 B = 0;
        int certified = 0;
        string method;
        string notes;
    };

    auto split_csv = [](const string& line) {
        vector<string> out;
        string cur;
        bool in_quotes = false;
        for (size_t i = 0; i < line.size(); ++i) {
            char c = line[i];
            if (c == '"') {
                if (in_quotes && i + 1 < line.size() && line[i + 1] == '"') {
                    cur.push_back('"');
                    ++i;
                } else {
                    in_quotes = !in_quotes;
                }
            } else if (c == ',' && !in_quotes) {
                out.push_back(cur);
                cur.clear();
            } else {
                cur.push_back(c);
            }
        }
        out.push_back(cur);
        return out;
    };

    auto csv_escape = [](const string& s) {
        bool need = false;
        for (char c : s) {
            if (c == '"' || c == ',' || c == '\n' || c == '\r') { need = true; break; }
        }
        if (!need) return s;
        string out;
        out.push_back('"');
        for (char c : s) {
            if (c == '"') out += "\"\"";
            else out.push_back(c);
        }
        out.push_back('"');
        return out;
    };

    ifstream fin(in_csv);
    if (!fin) {
        cerr << "Could not open interface CSV: " << in_csv << "\n";
        return 2;
    }
    string header;
    if (!getline(fin, header)) {
        cerr << "Empty interface CSV: " << in_csv << "\n";
        return 2;
    }
    vector<Row> rows;
    rows.reserve(6000);
    string line;
    while (getline(fin, line)) {
        if (line.empty()) continue;
        auto f = split_csv(line);
        if (f.size() < 6) continue;
        Row r;
        r.r = stoi(f[0]);
        r.R = stoi(f[1]);
        r.B = stoull(f[2]);
        r.certified = stoi(f[3]);
        r.method = f[4];
        r.notes = f[5];
        rows.push_back(std::move(r));
    }

    unordered_map<InterfaceKey, size_t, InterfaceKeyHash> row_index;
    row_index.reserve(rows.size() * 2);
    for (size_t i = 0; i < rows.size(); ++i) {
        row_index[InterfaceKey{rows[i].r, rows[i].B}] = i;
    }

    
    vector<int> seen(rows.size(), 0);

    cout << "Interface cert (bounded-gap representatives): gmax=" << gmax
         << " max_trailing_shift=" << max_trailing_shift << "\n";
    cout << "Input: " << in_csv << " (" << rows.size() << " rows)\n";

    
    const u64 base = (u64)(gmax + 1);
    u64 total_possible = 1;
    for (int i = 0; i < 7; ++i) total_possible *= base;

    atomic<bool> found(false);
    mutex ce_mutex;
    struct Counterexample {
        u64 n = 0;
        array<int, 7> gaps{};
        int l = -1;
        int ls = -1;
        int lambda = -1;
        int gamma = -1;
        int r = -1;
        u64 B = 0;
    };
    Counterexample ce;

    atomic<u64> progress_total(0);
    atomic<bool> progress_done(false);
    thread reporter([&]() {
        using namespace std::chrono_literals;
        while (!progress_done.load(memory_order_relaxed)) {
            this_thread::sleep_for(5s);
            u64 p = progress_total.load(memory_order_relaxed);
            if (p > total_possible) p = total_possible;
            cerr << "[interface_cert] progress: ~" << p << "/" << total_possible << " ("
                 << (100.0 * (double)p / (double)total_possible) << "%)\n";
            cerr.flush();
        }
    });

    unsigned int T = std::thread::hardware_concurrency();
    if (T == 0) T = 4;
    constexpr u64 PROGRESS_BATCH = 16;

    auto worker = [&](u64 start, u64 step) {
        u64 local_total = 0;
        for (u64 idx = start; idx < total_possible && !found.load(memory_order_relaxed); idx += step) {
            ++local_total;
            if ((local_total % PROGRESS_BATCH) == 0) progress_total.fetch_add(PROGRESS_BATCH, memory_order_relaxed);

            u64 tmp = idx;
            int g1 = (int)(tmp % base); tmp /= base;
            int g2 = (int)(tmp % base); tmp /= base;
            int g3 = (int)(tmp % base); tmp /= base;
            int g4 = (int)(tmp % base); tmp /= base;
            int g5 = (int)(tmp % base); tmp /= base;
            int g6 = (int)(tmp % base); tmp /= base;
            int t = (int)(tmp % base);

            int gamma = max({g1, g2, g3, g4, g5, g6, t});
            if (gamma < 4) continue; 
            if (t > max_trailing_shift) continue;

            
            u64 x = 1;
            x = (x << (g1 + 1)) + 1;
            x = (x << (g2 + 1)) + 1;
            x = (x << (g3 + 1)) + 1;
            x = (x << (g4 + 1)) + 1;
            x = (x << (g5 + 1)) + 1;
            x = (x << (g6 + 1)) + 1;
            u64 n = x << t;

            int r_if = -1;
            u64 B_if = 0;
            if (!compute_first_gap_interface_shifted(n, r_if, B_if)) {
                continue; 
            }

            auto it = row_index.find(InterfaceKey{r_if, B_if});
            if (it == row_index.end()) {
                
                continue;
            }
            seen[it->second] = 1;

            int lambda = bit_length(n) - 1;
            int ls = min_star_chain_length(n);
            int target = ls - 2;
            
            if (target < lambda + 3) continue;
            if (target >= 0 && has_chain_leq(n, target)) {
                lock_guard<mutex> lock(ce_mutex);
                if (!found.load(memory_order_relaxed)) {
                    auto [l_exact, ls_exact] = compute_chain_lengths(n);
                    ce = Counterexample{n, {g1, g2, g3, g4, g5, g6, t}, l_exact, ls_exact, lambda, gamma, r_if, B_if};
                    found.store(true, memory_order_relaxed);
                }
                return;
            }
        }
        progress_total.fetch_add(local_total % PROGRESS_BATCH, memory_order_relaxed);
    };

    vector<std::thread> threads;
    threads.reserve(T);
    for (unsigned int i = 0; i < T; ++i) threads.emplace_back(worker, (u64)i, (u64)T);
    for (auto& th : threads) th.join();
    progress_done.store(true, memory_order_relaxed);
    reporter.join();

    if (found.load(memory_order_relaxed)) {
        cerr << "COUNTEREXAMPLE FOUND at n=" << ce.n << "\n";
        cerr << "interface (r,B)=(" << ce.r << "," << ce.B << ")\n";
        cerr << "gaps(g1..g6,t)=(" << ce.gaps[0] << "," << ce.gaps[1] << "," << ce.gaps[2] << "," << ce.gaps[3] << ","
             << ce.gaps[4] << "," << ce.gaps[5] << "," << ce.gaps[6] << ")\n";
        cerr << "gamma=" << ce.gamma << " bits=" << bit_length(ce.n) << " wt=" << popcount(ce.n) << "\n";
        cerr << "lambda=" << ce.lambda << " l(n)=" << ce.l << " l*(n)=" << ce.ls << " ex=" << (ce.ls - ce.l) << "\n";
        cerr << "bin=" << std::bitset<64>(ce.n).to_string().substr(64 - bit_length(ce.n)) << "\n";
        
        auto it = row_index.find(InterfaceKey{ce.r, ce.B});
        if (it != row_index.end()) {
            auto& row = rows[it->second];
            row.certified = 0;
            row.method = "bounded-gap-interface-scan(WITNESS_FOUND)";
            std::ostringstream oss;
            oss << row.notes << "; witness_n=" << ce.n << "; witness_ex=" << (ce.ls - ce.l)
                << "; witness_gaps=(" << ce.gaps[0] << "," << ce.gaps[1] << "," << ce.gaps[2] << "," << ce.gaps[3] << ","
                << ce.gaps[4] << "," << ce.gaps[5] << "," << ce.gaps[6] << ")";
            row.notes = oss.str();
        }
        
    } else {
        
        string method = "bounded-gap-interface-scan(gmax=" + std::to_string(gmax) + ",tmax=" + std::to_string(max_trailing_shift) + ")";
        for (size_t i = 0; i < rows.size(); ++i) {
            if (!seen[i]) continue;
            rows[i].certified = 1;
            rows[i].method = method;
            rows[i].notes = rows[i].notes + "; certified_by=chain_solver(--certify-interface-gapbound)";
        }
    }

    
    size_t uncovered = 0;
    size_t certified = 0;
    for (size_t i = 0; i < rows.size(); ++i) {
        if (!seen[i]) uncovered++;
        if (rows[i].certified == 1) certified++;
    }
    cout << "Seen interfaces: " << (rows.size() - uncovered) << "/" << rows.size() << "\n";
    cout << "Certified rows: " << certified << "/" << rows.size() << "\n";

    ofstream fout(out_csv);
    if (!fout) {
        cerr << "Could not open output: " << out_csv << "\n";
        return 2;
    }
    fout << "r,R,B,certified,method,notes\n";
    for (const auto& row : rows) {
        fout << row.r << "," << row.R << "," << row.B << "," << row.certified << ","
             << csv_escape(row.method) << "," << csv_escape(row.notes) << "\n";
    }
    cout << "Wrote " << out_csv << "\n";
    return found.load(memory_order_relaxed) ? 1 : 0;
}

static int list_weight7_ex1(u64 N) {
    cout << "Listing n <= " << N << " over popcount(n)=7 with star-excess = 1 ..." << endl;
    cout << "CSV columns: n,bits,wt,lambda,l,ls,ex,s,t,max_internal_gap,gamma,bin" << endl;

    struct Row {
        u64 n = 0;
        int l = -1;
        int ls = -1;
        int lambda = -1;
        int s = -1;
        int t = -1;
        int g_int = -1;
        int gamma = -1;
    };

    vector<Row> rows;
    mutex rows_mutex;

    auto candidates = enumerate_weight7_up_to(N);
    cout << "Candidates (wt=7): " << candidates.size() << endl;

    unsigned int T = std::thread::hardware_concurrency();
    if (T == 0) T = 4;

    
    atomic<u64> progress_total(0);
    atomic<bool> progress_done(false);
    thread reporter([&]() {
        using namespace std::chrono_literals;
        while (!progress_done.load(memory_order_relaxed)) {
            this_thread::sleep_for(10s);
            u64 p = progress_total.load(memory_order_relaxed);
            u64 tot = candidates.size();
            if (p > tot) p = tot;
            cerr << "[list-wt7-ex1] progress: ~" << p << "/" << tot << " (" << (tot ? (100.0 * (double)p / (double)tot) : 100.0) << "%)\n";
            cerr.flush();
        }
    });

    constexpr u64 PROGRESS_BATCH = 256;

    auto worker = [&](u64 start, u64 step) {
        vector<Row> local;
        u64 local_total = 0;
        for (u64 idx = start; idx < candidates.size(); idx += step) {
            u64 n = candidates[idx];
            local_total++;
            if ((local_total % PROGRESS_BATCH) == 0) {
                progress_total.fetch_add(PROGRESS_BATCH, memory_order_relaxed);
            }
            auto [l, ls] = compute_chain_lengths(n);
            int ex = ls - l;
            if (ex != 1) continue;
            int lambda = bit_length(n) - 1;
            int s = l - lambda;
            int g_int = max_internal_gap_u64(n);
            int t = trailing_zeros_u64(n);
            int gamma = max(g_int, t);
            local.push_back(Row{n, l, ls, lambda, s, t, g_int, gamma});
        }

        progress_total.fetch_add(local_total % PROGRESS_BATCH, memory_order_relaxed);

        lock_guard<mutex> lock(rows_mutex);
        rows.insert(rows.end(), local.begin(), local.end());
    };

    vector<std::thread> threads;
    threads.reserve(T);
    for (u64 i = 0; i < T; ++i) {
        threads.emplace_back(worker, i, (u64)T);
    }
    for (auto& th : threads) th.join();
    progress_done.store(true, memory_order_relaxed);
    reporter.join();

    sort(rows.begin(), rows.end(), [](const Row& a, const Row& b) { return a.n < b.n; });

    for (const auto& row : rows) {
        int bits = bit_length(row.n);
        string bin = std::bitset<64>(row.n).to_string().substr(64 - bits);
        cout << row.n << "," << bits << "," << 7 << "," << row.lambda << "," << row.l << "," << row.ls << ",1," << row.s
             << "," << row.t << "," << row.g_int << "," << row.gamma << "," << bin << "\n";
    }

    cout << "Total ex=1 cases listed: " << rows.size() << "\n";
    return 0;
}

static int scan_weight7_internal_gap(u64 N, int min_gap, int ex_threshold) {
    cout << "Scanning n <= " << N << " over popcount(n)=7 with max_internal_gap >= " << min_gap << " ..." << endl;
    cout << "Reports any n with star-excess >= " << ex_threshold << " in this subset." << endl;

    auto candidates = enumerate_weight7_up_to(N);
    const u64 total_possible = (u64)candidates.size();
    cout << "Weight-7 candidates <= " << N << ": " << total_possible << "\n";

    struct Counterexample {
        u64 n = 0;
        int l = -1;
        int ls = -1;
        int lambda = -1;
        int g_int = -1;
        int t = -1;
    };

    atomic<bool> found(false);
    mutex ce_mutex;
    Counterexample ce;

    
    atomic<u64> progress_total(0);
    atomic<bool> progress_done(false);

    
    mutex stats_mutex;
    u64 total_checked = 0;
    u64 matched = 0;
    u64 ex0 = 0, ex1 = 0, ex2 = 0, ex3p = 0;
    u64 s3 = 0, s4 = 0, s5 = 0, s6p = 0;
    int min_s = 999;
    int max_s = -1;
    u64 min_s_n = 0;
    u64 max_s_n = 0;

    unsigned int T = std::thread::hardware_concurrency();
    if (T == 0) T = 4;

    thread reporter([&]() {
        using namespace std::chrono_literals;
        while (!progress_done.load(memory_order_relaxed)) {
            this_thread::sleep_for(10s);
            u64 p = progress_total.load(memory_order_relaxed);
            if (p > total_possible) p = total_possible;
            cerr << "[scan-wt7-internal] progress: ~" << p << "/" << total_possible << " ("
                 << (100.0 * (double)p / (double)total_possible) << "%)\n";
            cerr.flush();
        }
    });

    constexpr u64 PROGRESS_BATCH = 64;

    auto worker = [&](u64 start, u64 step) {
        u64 local_total = 0;
        u64 local_matched = 0;
        u64 local_ex0 = 0, local_ex1 = 0, local_ex2 = 0, local_ex3p = 0;
        u64 local_s3 = 0, local_s4 = 0, local_s5 = 0, local_s6p = 0;
        int local_min_s = 999;
        int local_max_s = -1;
        u64 local_min_s_n = 0;
        u64 local_max_s_n = 0;

        for (u64 idx = start; idx < total_possible && !found.load(memory_order_relaxed); idx += step) {
            u64 n = candidates[idx];
            local_total++;
            if ((local_total % PROGRESS_BATCH) == 0) progress_total.fetch_add(PROGRESS_BATCH, memory_order_relaxed);
            int g_int = max_internal_gap_u64(n);
            if (g_int < min_gap) continue;
            ++local_matched;

            auto [l, ls] = compute_chain_lengths(n);
            int lambda = bit_length(n) - 1;
            int s = l - lambda;
            int ex = ls - l;

            if (s < local_min_s) {
                local_min_s = s;
                local_min_s_n = n;
            }
            if (s > local_max_s) {
                local_max_s = s;
                local_max_s_n = n;
            }

            if (s == 3) ++local_s3;
            else if (s == 4) ++local_s4;
            else if (s == 5) ++local_s5;
            else if (s >= 6) ++local_s6p;

            if (ex == 0) ++local_ex0;
            else if (ex == 1) ++local_ex1;
            else if (ex == 2) ++local_ex2;
            else if (ex >= 3) ++local_ex3p;

            if (ex >= ex_threshold) {
                lock_guard<mutex> lock(ce_mutex);
                if (!found.load(memory_order_relaxed)) {
                    ce = Counterexample{
                        n,
                        l,
                        ls,
                        lambda,
                        g_int,
                        trailing_zeros_u64(n),
                    };
                    found.store(true, memory_order_relaxed);
                }
                return;
            }
        }

        progress_total.fetch_add(local_total % PROGRESS_BATCH, memory_order_relaxed);

        lock_guard<mutex> lock(stats_mutex);
        total_checked += local_total;
        matched += local_matched;
        ex0 += local_ex0;
        ex1 += local_ex1;
        ex2 += local_ex2;
        ex3p += local_ex3p;
        s3 += local_s3;
        s4 += local_s4;
        s5 += local_s5;
        s6p += local_s6p;
        if (local_min_s < min_s) {
            min_s = local_min_s;
            min_s_n = local_min_s_n;
        }
        if (local_max_s > max_s) {
            max_s = local_max_s;
            max_s_n = local_max_s_n;
        }
    };

    vector<std::thread> threads;
    threads.reserve(T);
    for (u64 i = 0; i < T; ++i) {
        threads.emplace_back(worker, i, (u64)T);
    }
    for (auto& th : threads) th.join();
    progress_done.store(true, memory_order_relaxed);
    reporter.join();

    if (found.load(memory_order_relaxed)) {
        cout << "COUNTEREXAMPLE FOUND at n=" << ce.n << "\n";
        cout << "bits=" << bit_length(ce.n) << " wt=" << popcount(ce.n) << "\n";
        cout << "t=" << ce.t << " max_internal_gap=" << ce.g_int << "\n";
        cout << "lambda=" << ce.lambda << " l(n)=" << ce.l << " l*(n)=" << ce.ls << " ex=" << (ce.ls - ce.l) << "\n";
        cout << "bin=" << std::bitset<64>(ce.n).to_string().substr(64 - bit_length(ce.n)) << "\n";
        return 1;
    }

    cout << "OK: no n found with star-excess >= " << ex_threshold
         << " among popcount(n)=7 candidates <= " << N << " with max_internal_gap >= " << min_gap << "\n";
    cout << "Stats: examined " << total_checked << " candidates; matched=" << matched << "\n";
    cout << "Star-excess counts: ex=0: " << ex0 << "  ex=1: " << ex1 << "  ex=2: " << ex2 << "  ex>=3: " << ex3p << "\n";
    cout << "Small-steps s=l-lambda counts: s=3: " << s3 << "  s=4: " << s4 << "  s=5: " << s5 << "  s>=6: " << s6p << "\n";
    if (matched > 0) {
        cout << "Observed s-range: min_s=" << min_s << " at n=" << min_s_n << " ; max_s=" << max_s << " at n=" << max_s_n << "\n";
    }
    return 0;
}

static int scan_weight7_gamma(u64 N, int gamma_min, int ex_threshold) {
    cout << "Scanning n <= " << N << " over popcount(n)=7 with gamma=max_gap >= " << gamma_min << " ..." << endl;
    cout << "Reports any n with star-excess >= " << ex_threshold << " in this subset." << endl;

    auto candidates = enumerate_weight7_up_to(N);
    const u64 total_possible = (u64)candidates.size();
    cout << "Weight-7 candidates <= " << N << ": " << total_possible << "\n";

    struct Counterexample {
        u64 n = 0;
        int l = -1;
        int ls = -1;
        int lambda = -1;
        int g_int = -1;
        int t = -1;
        int gamma = -1;
    };

    atomic<bool> found(false);
    mutex ce_mutex;
    Counterexample ce;

    
    atomic<u64> progress_total(0);
    atomic<bool> progress_done(false);

    
    mutex stats_mutex;
    u64 total_checked = 0;
    u64 matched = 0;
    u64 ex0 = 0, ex1 = 0, ex2 = 0, ex3p = 0;
    u64 s3 = 0, s4 = 0, s5 = 0, s6p = 0;
    u64 max_ex = 0;
    u64 max_ex_n = 0;

    unsigned int T = std::thread::hardware_concurrency();
    if (T == 0) T = 4;

    thread reporter([&]() {
        using namespace std::chrono_literals;
        while (!progress_done.load(memory_order_relaxed)) {
            this_thread::sleep_for(10s);
            u64 p = progress_total.load(memory_order_relaxed);
            if (p > total_possible) p = total_possible;
            cerr << "[scan-wt7-gamma] progress: ~" << p << "/" << total_possible << " ("
                 << (100.0 * (double)p / (double)total_possible) << "%)\n";
            cerr.flush();
        }
    });

    constexpr u64 PROGRESS_BATCH = 64;

    auto worker = [&](u64 start, u64 step) {
        u64 local_total = 0;
        u64 local_matched = 0;
        u64 local_ex0 = 0, local_ex1 = 0, local_ex2 = 0, local_ex3p = 0;
        u64 local_s3 = 0, local_s4 = 0, local_s5 = 0, local_s6p = 0;
        u64 local_max_ex = 0;
        u64 local_max_ex_n = 0;

        for (u64 idx = start; idx < total_possible && !found.load(memory_order_relaxed); idx += step) {
            u64 n = candidates[idx];
            local_total++;
            if ((local_total % PROGRESS_BATCH) == 0) progress_total.fetch_add(PROGRESS_BATCH, memory_order_relaxed);

            int g_int = max_internal_gap_u64(n);
            int t = trailing_zeros_u64(n);
            int gamma = max(g_int, t);
            if (gamma < gamma_min) continue;
            ++local_matched;

            auto [l, ls] = compute_chain_lengths(n);
            int lambda = bit_length(n) - 1;
            int ex = ls - l;

            if (local_matched == 1) {
                local_max_ex = (u64)ex;
                local_max_ex_n = n;
            }

            if ((u64)ex > local_max_ex) {
                local_max_ex = (u64)ex;
                local_max_ex_n = n;
            }

            if (ex == 0) ++local_ex0;
            else if (ex == 1) ++local_ex1;
            else if (ex == 2) ++local_ex2;
            else if (ex >= 3) ++local_ex3p;

            int s = l - lambda;
            if (s == 3) ++local_s3;
            else if (s == 4) ++local_s4;
            else if (s == 5) ++local_s5;
            else if (s >= 6) ++local_s6p;

            if (ex >= ex_threshold) {
                lock_guard<mutex> lock(ce_mutex);
                if (!found.load(memory_order_relaxed)) {
                    ce = Counterexample{n, l, ls, lambda, g_int, t, gamma};
                    found.store(true, memory_order_relaxed);
                }
                return;
            }
        }

        progress_total.fetch_add(local_total % PROGRESS_BATCH, memory_order_relaxed);

        lock_guard<mutex> lock(stats_mutex);
        total_checked += local_total;
        matched += local_matched;
        ex0 += local_ex0;
        ex1 += local_ex1;
        ex2 += local_ex2;
        ex3p += local_ex3p;
        s3 += local_s3;
        s4 += local_s4;
        s5 += local_s5;
        s6p += local_s6p;
        if (local_max_ex > max_ex || (local_max_ex == max_ex && max_ex_n == 0 && local_max_ex_n != 0)) {
            max_ex = local_max_ex;
            max_ex_n = local_max_ex_n;
        }
    };

    vector<std::thread> threads;
    threads.reserve(T);
    for (u64 i = 0; i < T; ++i) {
        threads.emplace_back(worker, i, (u64)T);
    }
    for (auto& th : threads) th.join();
    progress_done.store(true, memory_order_relaxed);
    reporter.join();

    if (found.load(memory_order_relaxed)) {
        cout << "COUNTEREXAMPLE FOUND at n=" << ce.n << "\n";
        cout << "bits=" << bit_length(ce.n) << " wt=" << popcount(ce.n) << "\n";
        cout << "t=" << ce.t << " max_internal_gap=" << ce.g_int << " gamma=" << ce.gamma << "\n";
        cout << "lambda=" << ce.lambda << " l(n)=" << ce.l << " l*(n)=" << ce.ls << " ex=" << (ce.ls - ce.l) << "\n";
        cout << "bin=" << std::bitset<64>(ce.n).to_string().substr(64 - bit_length(ce.n)) << "\n";
        return 1;
    }

    cout << "OK: no n found with star-excess >= " << ex_threshold
         << " among popcount(n)=7 candidates <= " << N << " with gamma >= " << gamma_min << "\n";
    cout << "Stats: examined " << total_checked << " candidates; matched=" << matched << "\n";
    cout << "Star-excess counts: ex=0: " << ex0 << "  ex=1: " << ex1 << "  ex=2: " << ex2 << "  ex>=3: " << ex3p << "\n";
    cout << "Small-steps s=l-lambda counts: s=3: " << s3 << "  s=4: " << s4 << "  s=5: " << s5 << "  s>=6: " << s6p << "\n";
    cout << "Max observed star-excess: " << max_ex << " at n=" << max_ex_n << "\n";
    return 0;
}







static int scan_weight7_gamma_ex2_fast(u64 N, int gamma_min) {
    cout << "Fast scan: n <= " << N << " over popcount(n)=7 with gamma=max_gap >= " << gamma_min << " ..." << endl;
    cout << "Reports any n with star-excess >= 2 in this subset (may be much faster than exact scans)." << endl;

    auto candidates = enumerate_weight7_up_to(N);
    const u64 total_possible = (u64)candidates.size();
    cout << "Weight-7 candidates <= " << N << ": " << total_possible << "\n";

    struct Counterexample {
        u64 n = 0;
        int l = -1;
        int ls = -1;
        int lambda = -1;
        int g_int = -1;
        int t = -1;
        int gamma = -1;
    };

    atomic<bool> found(false);
    mutex ce_mutex;
    Counterexample ce;

    
    atomic<u64> progress_total(0);
    atomic<bool> progress_done(false);

    unsigned int T = std::thread::hardware_concurrency();
    if (T == 0) T = 4;

    thread reporter([&]() {
        using namespace std::chrono_literals;
        while (!progress_done.load(memory_order_relaxed)) {
            this_thread::sleep_for(10s);
            u64 p = progress_total.load(memory_order_relaxed);
            if (p > total_possible) p = total_possible;
            cerr << "[scan-wt7-gamma-fast] progress: ~" << p << "/" << total_possible << " ("
                 << (100.0 * (double)p / (double)total_possible) << "%)\n";
            cerr.flush();
        }
    });

    constexpr u64 PROGRESS_BATCH = 64;

    auto worker = [&](u64 start, u64 step) {
        u64 local_total = 0;
        for (u64 idx = start; idx < total_possible && !found.load(memory_order_relaxed); idx += step) {
            u64 n = candidates[idx];
            local_total++;
            if ((local_total % PROGRESS_BATCH) == 0) progress_total.fetch_add(PROGRESS_BATCH, memory_order_relaxed);

            int g_int = max_internal_gap_u64(n);
            int t = trailing_zeros_u64(n);
            int gamma = max(g_int, t);
            if (gamma < gamma_min) continue;

            int lambda = bit_length(n) - 1;

            
            int ls = min_star_chain_length(n);

            
            
            if (ls - 2 < lambda + 3) continue;

            
            if (!has_chain_leq(n, ls - 2)) continue;

            
            int l = min_chain_length(n, ls);
            if (ls - l < 2) continue;  

            lock_guard<mutex> lock(ce_mutex);
            if (!found.load(memory_order_relaxed)) {
                ce = Counterexample{n, l, ls, lambda, g_int, t, gamma};
                found.store(true, memory_order_relaxed);
            }
            return;
        }
        progress_total.fetch_add(local_total % PROGRESS_BATCH, memory_order_relaxed);
    };

    vector<std::thread> threads;
    threads.reserve(T);
    for (u64 i = 0; i < T; ++i) {
        threads.emplace_back(worker, i, (u64)T);
    }
    for (auto& th : threads) th.join();
    progress_done.store(true, memory_order_relaxed);
    reporter.join();

    if (found.load(memory_order_relaxed)) {
        cout << "COUNTEREXAMPLE FOUND at n=" << ce.n << "\n";
        cout << "bits=" << bit_length(ce.n) << " wt=" << popcount(ce.n) << "\n";
        cout << "t=" << ce.t << " max_internal_gap=" << ce.g_int << " gamma=" << ce.gamma << "\n";
        cout << "lambda=" << ce.lambda << " l(n)=" << ce.l << " l*(n)=" << ce.ls << " ex=" << (ce.ls - ce.l) << "\n";
        cout << "bin=" << std::bitset<64>(ce.n).to_string().substr(64 - bit_length(ce.n)) << "\n";
        return 1;
    }

    cout << "OK: no n found with star-excess >= 2 among popcount(n)=7 candidates <= " << N
         << " with gamma >= " << gamma_min << "\n";
    return 0;
}




int main(int argc, char** argv) {
    auto print_usage = [&]() {
        cout << "Addition Chain Solver (Parallel + Star-First)" << endl;
        cout << "Usage:" << endl;
        cout << "  " << argv[0] << " <n>        Compute l(n) and l*(n)" << endl;
        cout << "  " << argv[0] << " --witness <n>  Print witness chains (general + star)" << endl;
        cout << "  " << argv[0] << " --help    Show this help" << endl;
        cout << "  " << argv[0] << " --verify   Full 345-outlier verification" << endl;
        cout << "  " << argv[0] << " --verify-compact   Verify all 16384 compact-domain exponents (writes compact_brauer_results.csv)" << endl;
        cout << "  " << argv[0] << " --scan-lt1 [N] Scan n<=N for defect(n)<1 classification (default N=20000)" << endl;
        cout << "  " << argv[0] << " --scan-d4 [N]  Scan n<=N for defect<d4 classification (default N=20000)" << endl;
        cout << "  " << argv[0] << " --scan-d5 [N]  Scan n<=N for defect<d5 classification (default N=20000)" << endl;
        cout << "  " << argv[0] << " --list-d5 [N]  List n<=N with defect<d5 (d5 = 6 - log2 21), as CSV (default N=5000)" << endl;
        cout << "  " << argv[0] << " --stats-d5 [N]  Summary stats for defect<d5 up to N (default N=20000)" << endl;
        cout << "  " << argv[0] << " --check-envelope [R]  Exhaustively verify the Fibonacci envelope for r<=R (default R=22)" << endl;
        cout << "  " << argv[0] << " --check-rays [B]  Verify l(7*2^b), l(15*2^b), l(27*2^b) for b<=B (default B=25)" << endl;
        cout << "  " << argv[0] << " --check-d5-rays [B]  Verify l(m*2^b) for d5 cores m in {11,23,43,45,85} (default B=25)" << endl;
        cout << "  " << argv[0] << " --scan-wt7 [N] Scan n<=N over popcount=7 for star-excess >=2 (default N=65535)" << endl;
        cout << "  " << argv[0] << " --verify-wt7-gapbound [gmax] [sparse_only] Enumerate all wt=7 gap vectors with entries<=gmax and check for star-excess>=2; if sparse_only=1 restrict to gamma>=4 (default gmax=4, sparse_only=1)" << endl;
        cout << "  " << argv[0] << " --certify-interface-gapbound [gmax] [tmax] [in_csv] [out_csv] Interface-focused bounded-gap scan; updates scholz/interface_cert.csv by marking interfaces seen in the scan as certified (default gmax=4, tmax=3, in=scholz/interface_cert.csv, out=scholz/interface_cert.csv)" << endl;
        cout << "  " << argv[0] << " --list-wt7-ex1 [N] List n<=N over popcount=7 with star-excess = 1 (default N=65535)" << endl;
        cout << "  " << argv[0] << " --scan-wt7-internal [N] [min_gap] [ex] Scan n<=N over popcount=7 with max_internal_gap>=min_gap for star-excess>=ex (default N=1000000, min_gap=4, ex=2)" << endl;
        cout << "  " << argv[0] << " --scan-wt7-gamma [N] [gamma] [ex] Scan n<=N over popcount=7 with gamma=max_gap>=gamma for star-excess>=ex (default N=100000, gamma=4, ex=2)" << endl;
        cout << "  " << argv[0] << " --scan-wt7-gamma-fast [N] [gamma] Fast scan n<=N over popcount=7 with gamma>=gamma for star-excess>=2 (default N=1000000, gamma=4)" << endl;
    };

    if (argc == 1) {
        print_usage();
        return 0;
    }
    
    string arg = argv[1];

    if (arg == "--help" || arg == "-h") {
        print_usage();
        return 0;
    }

    if (arg == "--witness") {
        if (argc < 3) {
            cerr << "Missing n for --witness\n\n";
            print_usage();
            return 2;
        }
        u64 n = 0;
        try {
            n = stoull(argv[2]);
        } catch (const std::exception& e) {
            cerr << "Invalid integer argument: " << argv[2] << " (" << e.what() << ")\n\n";
            print_usage();
            return 2;
        }
        cout << "n = " << n << "\n";
        cout << "Bits: " << bit_length(n) << ", Weight: " << popcount(n) << "\n";
        cout << "t=" << trailing_zeros_u64(n) << " max_internal_gap=" << max_internal_gap_u64(n)
             << " gamma=" << max_gap_u64(n) << "\n";

        vector<u64> chain_general;
        vector<u64> chain_star;
        auto t0 = chrono::high_resolution_clock::now();
        auto [l, ls] = compute_chain_lengths_with_witness(n, chain_general, chain_star);
        auto t1 = chrono::high_resolution_clock::now();
        double elapsed = chrono::duration<double>(t1 - t0).count();

        cout << "\n";
        cout << "l(n)  = " << l << "\n";
        cout << "l*(n) = " << ls << "\n";
        cout << "ex    = " << (ls - l) << "\n";
        cout << "Time: " << elapsed << "s\n";

        cout << "\nGENERAL CHAIN\n";
        print_chain_with_certificate(chain_general);
        cout << "\nSTAR CHAIN\n";
        print_chain_with_certificate(chain_star);
        return 0;
    }
    
    if (arg == "--verify") {
        run_verification();
        return 0;
    }

    if (arg == "--verify-compact") {
        run_compact_domain_verification();
        return 0;
    }

    if (arg == "--scan-wt7") {
        u64 N = 65535;
        if (argc >= 3) N = stoull(argv[2]);
        return scan_weight7(N);
    }

    if (arg == "--verify-wt7-gapbound") {
        int gmax = 4;
        int sparse_only = 1;
        if (argc >= 3) gmax = stoi(argv[2]);
        if (argc >= 4) sparse_only = stoi(argv[3]);
        if (gmax < 0) gmax = 0;
        return verify_wt7_gapbound(gmax, sparse_only != 0);
    }

    if (arg == "--certify-interface-gapbound") {
        int gmax = 4;
        int tmax = 3;
        string in_csv = "scholz/interface_cert.csv";
        string out_csv = "scholz/interface_cert.csv";
        if (argc >= 3) gmax = stoi(argv[2]);
        if (argc >= 4) tmax = stoi(argv[3]);
        if (argc >= 5) in_csv = argv[4];
        if (argc >= 6) out_csv = argv[5];
        if (gmax < 0) gmax = 0;
        if (tmax < 0) tmax = 0;
        if (tmax > 3) tmax = 3;
        return certify_interface_from_gapbound(gmax, tmax, in_csv, out_csv);
    }

    if (arg == "--list-wt7-ex1") {
        u64 N = 65535;
        if (argc >= 3) N = stoull(argv[2]);
        return list_weight7_ex1(N);
    }
    
    if (arg == "--scan-wt7-internal") {
        u64 N = 1000000;
        int min_gap = 4;
        int ex_threshold = 2;
        if (argc >= 3) N = stoull(argv[2]);
        if (argc >= 4) min_gap = stoi(argv[3]);
        if (argc >= 5) ex_threshold = stoi(argv[4]);
        return scan_weight7_internal_gap(N, min_gap, ex_threshold);
    }
    
    if (arg == "--scan-wt7-gamma") {
        u64 N = 100000;
        int gamma_min = 4;
        int ex_threshold = 2;
        if (argc >= 3) N = stoull(argv[2]);
        if (argc >= 4) gamma_min = stoi(argv[3]);
        if (argc >= 5) ex_threshold = stoi(argv[4]);
        return scan_weight7_gamma(N, gamma_min, ex_threshold);
    }
    if (arg == "--scan-wt7-gamma-fast") {
        u64 N = (argc >= 3) ? stoull(argv[2]) : 1000000ULL;
        int gamma_min = (argc >= 4) ? stoi(argv[3]) : 4;
        return scan_weight7_gamma_ex2_fast(N, gamma_min);
    }

    if (arg == "--scan-lt1") {
        u64 N = 20000;
        if (argc >= 3) N = stoull(argv[2]);
        cout << "Scanning n <= " << N << " for defect(n) < 1 classification..." << endl;
        cout << "lt1 condition checked via: 2^l(n) < 2 * n" << endl;

        struct Counterexample {
            u64 n = 0;
            int l = -1;
            int ls = -1;
            bool lt = false;
            bool cls = false;
        };

        atomic<bool> found(false);
        mutex ce_mutex;
        Counterexample ce;

        unsigned int T = std::thread::hardware_concurrency();
        if (T == 0) T = 4;

        auto worker = [&](u64 start, u64 step) {
            for (u64 n = start; n <= N && !found.load(memory_order_relaxed); n += step) {
                int l = min_chain_length(n);
                bool lt = defect_lt_1(n, l);
                bool cls = in_lt1_classification(n);
                if (lt != cls) {
                    int ls = min_star_chain_length(n);
                    {
                        lock_guard<mutex> lock(ce_mutex);
                        if (!found.load(memory_order_relaxed)) {
                            ce = Counterexample{n, l, ls, lt, cls};
                            found.store(true, memory_order_relaxed);
                        }
                    }
                    return;
                }
            }
        };

        vector<std::thread> threads;
        threads.reserve(T);
        for (u64 i = 0; i < T; ++i) {
            threads.emplace_back(worker, 1 + i, (u64)T);
        }
        for (auto& th : threads) th.join();

        if (found.load(memory_order_relaxed)) {
            cout << "COUNTEREXAMPLE at n=" << ce.n << "\n";
            cout << "bits=" << bit_length(ce.n) << " wt=" << popcount(ce.n) << "\n";
            cout << "l(n)=" << ce.l << " l*(n)=" << ce.ls << "\n";
            cout << "defect<1? " << (ce.lt ? "yes" : "no") << "  classification? " << (ce.cls ? "yes" : "no") << "\n";
            return 1;
        }

        cout << "OK: no counterexamples for n <= " << N << "\n";
        return 0;
    }

    if (arg == "--scan-d4") {
        u64 N = 20000;
        if (argc >= 3) N = stoull(argv[2]);
        cout << "Scanning n <= " << N << " for defect(n) < d4 classification..." << endl;
        cout << "d4 condition checked via: 13 * 2^l(n) < 32 * n" << endl;

        struct Counterexample {
            u64 n = 0;
            int l = -1;
            int ls = -1;
            bool lt = false;
            bool cls = false;
        };

        atomic<bool> found(false);
        mutex ce_mutex;
        Counterexample ce;

        unsigned int T = std::thread::hardware_concurrency();
        if (T == 0) T = 4;

        auto worker = [&](u64 start, u64 step) {
            for (u64 n = start; n <= N && !found.load(memory_order_relaxed); n += step) {
                int l = min_chain_length(n);
                bool lt = defect_lt_d4(n, l);
                bool cls = in_d4_classification(n);
                if (lt != cls) {
                    int ls = min_star_chain_length(n);
                    {
                        lock_guard<mutex> lock(ce_mutex);
                        if (!found.load(memory_order_relaxed)) {
                            ce = Counterexample{n, l, ls, lt, cls};
                            found.store(true, memory_order_relaxed);
                        }
                    }
                    return;
                }
            }
        };

        vector<std::thread> threads;
        threads.reserve(T);
        for (u64 i = 0; i < T; ++i) {
            threads.emplace_back(worker, 1 + i, (u64)T);
        }
        for (auto& th : threads) th.join();

        if (found.load(memory_order_relaxed)) {
            cout << "COUNTEREXAMPLE at n=" << ce.n << "\n";
            cout << "bits=" << bit_length(ce.n) << " wt=" << popcount(ce.n) << " odd_part=" << odd_part(ce.n) << "\n";
            cout << "l(n)=" << ce.l << " l*(n)=" << ce.ls << "\n";
            cout << "defect<d4? " << (ce.lt ? "yes" : "no") << "  classification? " << (ce.cls ? "yes" : "no") << "\n";
            return 1;
        }

        cout << "OK: no counterexamples for n <= " << N << "\n";
        return 0;
    }

    if (arg == "--scan-d5") {
        u64 N = 20000;
        if (argc >= 3) N = stoull(argv[2]);
        cout << "Scanning n <= " << N << " for defect(n) < d5 classification..." << endl;
        cout << "d5 condition checked via: 21 * 2^l(n) < 64 * n" << endl;

        struct Counterexample {
            u64 n = 0;
            int l = -1;
            int ls = -1;
            bool lt = false;
            bool cls = false;
        };

        atomic<bool> found(false);
        mutex ce_mutex;
        Counterexample ce;

        unsigned int T = std::thread::hardware_concurrency();
        if (T == 0) T = 4;

        auto worker = [&](u64 start, u64 step) {
            for (u64 n = start; n <= N && !found.load(memory_order_relaxed); n += step) {
                int l = min_chain_length(n);
                bool lt = defect_lt_d5(n, l);
                bool cls = in_d5_classification(n);
                if (lt != cls) {
                    int ls = min_star_chain_length(n);
                    {
                        lock_guard<mutex> lock(ce_mutex);
                        if (!found.load(memory_order_relaxed)) {
                            ce = Counterexample{n, l, ls, lt, cls};
                            found.store(true, memory_order_relaxed);
                        }
                    }
                    return;
                }
            }
        };

        vector<std::thread> threads;
        threads.reserve(T);
        for (u64 i = 0; i < T; ++i) {
            threads.emplace_back(worker, 1 + i, (u64)T);
        }
        for (auto& th : threads) th.join();

        if (found.load(memory_order_relaxed)) {
            cout << "COUNTEREXAMPLE at n=" << ce.n << "\n";
            cout << "bits=" << bit_length(ce.n) << " wt=" << popcount(ce.n) << " odd_part=" << odd_part(ce.n) << "\n";
            cout << "l(n)=" << ce.l << " l*(n)=" << ce.ls << "\n";
            cout << "defect<d5? " << (ce.lt ? "yes" : "no") << "  classification? " << (ce.cls ? "yes" : "no") << "\n";
            return 1;
        }

        cout << "OK: no counterexamples for n <= " << N << "\n";
        return 0;
    }

    if (arg == "--list-d5") {
        u64 N = 5000;
        if (argc >= 3) N = stoull(argv[2]);
        struct Row {
            u64 n;
            int l;
            int fl;
            int s;
            int wt;
            u64 op;
        };

        unsigned int T = std::thread::hardware_concurrency();
        if (T == 0) T = 4;

        vector<vector<Row>> buckets(T);

        auto worker = [&](u64 start, u64 step, vector<Row>& out) {
            for (u64 n = start; n <= N; n += step) {
                int l = min_chain_length(n);
                if (!defect_lt_d5(n, l)) continue;
                int fl = bit_length(n) - 1;
                int s = l - fl;
                int wt = popcount(n);
                u64 op = odd_part(n);
                out.push_back(Row{n, l, fl, s, wt, op});
            }
        };

        vector<std::thread> threads;
        threads.reserve(T);
        for (u64 i = 0; i < T; ++i) {
            threads.emplace_back(worker, 1 + i, (u64)T, std::ref(buckets[i]));
        }
        for (auto& th : threads) th.join();

        vector<Row> rows;
        size_t total = 0;
        for (const auto& b : buckets) total += b.size();
        rows.reserve(total);
        for (auto& b : buckets) {
            rows.insert(rows.end(), b.begin(), b.end());
        }
        sort(rows.begin(), rows.end(), [](const Row& a, const Row& b) { return a.n < b.n; });

        cout << "n,l,floorlog,s,wt,odd_part\n";
        for (const auto& r : rows) {
            cout << r.n << "," << r.l << "," << r.fl << "," << r.s << "," << r.wt << "," << r.op << "\n";
        }
        return 0;
    }

    if (arg == "--stats-d5") {
        u64 N = 20000;
        if (argc >= 3) N = stoull(argv[2]);
        cout << "Computing defect<d5 stats for n <= " << N << " (d5 = 6 - log2 21)..." << endl;
        cout << "d5 condition checked via: 21 * 2^l(n) < 64 * n" << endl;

        struct Stats {
            u64 count = 0;
            u64 wt_counts[65] = {0};
            unordered_set<u64> odd_parts_wt_ge_3;
        };

        unsigned int T = std::thread::hardware_concurrency();
        if (T == 0) T = 4;

        vector<Stats> stats(T);

        auto worker = [&](u64 start, u64 step, Stats& st) {
            for (u64 n = start; n <= N; n += step) {
                int l = min_chain_length(n);
                if (!defect_lt_d5(n, l)) continue;
                ++st.count;
                int wt = popcount(n);
                if (wt >= 0 && wt < 65) ++st.wt_counts[wt];
                if (wt >= 3) st.odd_parts_wt_ge_3.insert(odd_part(n));
            }
        };

        vector<std::thread> threads;
        threads.reserve(T);
        for (u64 i = 0; i < T; ++i) {
            threads.emplace_back(worker, 1 + i, (u64)T, std::ref(stats[i]));
        }
        for (auto& th : threads) th.join();

        u64 total_count = 0;
        u64 total_wt[65] = {0};
        unordered_set<u64> odd_parts;
        for (auto& st : stats) {
            total_count += st.count;
            for (int w = 0; w < 65; ++w) total_wt[w] += st.wt_counts[w];
            odd_parts.insert(st.odd_parts_wt_ge_3.begin(), st.odd_parts_wt_ge_3.end());
        }

        cout << "count=" << total_count << "\n";
        cout << "wt_counts:";
        for (int w = 0; w < 65; ++w) {
            if (total_wt[w] == 0) continue;
            cout << " wt" << w << "=" << total_wt[w];
        }
        cout << "\n";

        vector<u64> ops(odd_parts.begin(), odd_parts.end());
        sort(ops.begin(), ops.end());
        cout << "odd_parts_wt>=3(" << ops.size() << "):";
        for (u64 m : ops) cout << " " << m;
        cout << "\n";
        return 0;
    }

    if (arg == "--check-envelope") {
        int R = 22;
        if (argc >= 3) R = stoi(argv[2]);
        check_fib_envelope(R);
        return 0;
    }

    if (arg == "--check-rays") {
        int B = 25;
        if (argc >= 3) B = stoi(argv[2]);
        check_7_15_27_rays(B);
        return 0;
    }

    if (arg == "--check-d5-rays") {
        int B = 25;
        if (argc >= 3) B = stoi(argv[2]);
        check_d5_rays(B);
        return 0;
    }
    
    if (!arg.empty() && arg[0] == '-') {
        cerr << "Unknown option: " << arg << "\n\n";
        print_usage();
        return 2;
    }

    u64 n = 0;
    try {
        n = stoull(arg);
    } catch (const std::exception& e) {
        cerr << "Invalid integer argument: " << arg << " (" << e.what() << ")\n\n";
        print_usage();
        return 2;
    }
    cout << "n = " << n << endl;
    cout << "Bits: " << bit_length(n) << ", Weight: " << popcount(n) << endl;
    
    auto t0 = chrono::high_resolution_clock::now();
    auto [l, ls] = compute_chain_lengths(n);
    auto t1 = chrono::high_resolution_clock::now();
    
    double elapsed = chrono::duration<double>(t1 - t0).count();
    
    cout << endl;
    cout << "l(n)  = " << l << endl;
    cout << "l*(n) = " << ls << endl;
    cout << "Time: " << elapsed << "s" << endl;
    cout << (l == ls ? "=> BRAUER NUMBER" : "=> NOT A BRAUER NUMBER") << endl;
    
    return 0;
}
