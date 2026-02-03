#include <algorithm>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <queue>
#include <set>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>

using u32 = uint32_t;
using u64 = uint64_t;
using ResidueSet = std::vector<u32>;

static inline u64 mask_u64(int R) {
    return (R >= 32) ? 0xFFFFFFFFull : ((1ull << R) - 1ull);
}

static inline int bit_length_u32(u32 x) { return x ? 32 - __builtin_clz(x) : 0; }
static inline int popcount_u32(u32 x) { return __builtin_popcount(x); }
static inline int ctz_u32(u32 x) { return x ? __builtin_ctz(x) : 32; }

static inline u32 odd_part_mod_pow2(u32 x) {
    if (x == 0) return 0;
    return x >> ctz_u32(x);
}

static inline u32 odd_part_mod_2R(u32 x, int R) {
    if (R <= 0) return 0;
    if (x == 0) return 0;
    u32 mask = (R >= 32) ? 0xFFFFFFFFu : ((1u << R) - 1u);
    u32 r = x & mask;
    if (r == 0) return 0;
    r >>= ctz_u32(r);
    return r;
}

static ResidueSet normalize(ResidueSet S) {
    std::sort(S.begin(), S.end());
    S.erase(std::unique(S.begin(), S.end()), S.end());
    return S;
}

static ResidueSet close_under_doubling(ResidueSet S, int R) {
    S = normalize(std::move(S));
    const u64 mask = mask_u64(R);

    std::unordered_set<u32> seen;
    seen.reserve(S.size() * (size_t)(R + 2));
    for (u32 x : S) seen.insert(x);

    for (u32 x : S) {
        u32 y = x;
        for (int k = 0; k < R + 1; ++k) {
            (void)seen.insert(y);
            u32 nxt = (u32)(((u64)y << 1) & mask);
            if (nxt == y) break;
            y = nxt;
            if (y == 0) {
                seen.insert(0);
                break;
            }
        }
    }

    ResidueSet out;
    out.reserve(seen.size());
    for (u32 x : seen) out.push_back(x);
    return normalize(std::move(out));
}

static bool has_two_dirs(const ResidueSet& S, int R) {
    u32 d1 = 0;
    for (u32 x : S) {
        u32 d = odd_part_mod_2R(x, R);
        if (d == 0) continue;
        if (d1 == 0) d1 = d;
        else if (d != d1) return true;
    }
    return false;
}

static bool is_bad(const ResidueSet& S, u32 B, int R) {
    bool has_B = std::binary_search(S.begin(), S.end(), B);
    if (!has_B) return false;
    return has_two_dirs(S, R);
}

struct VecHash {
    size_t operator()(const ResidueSet& v) const noexcept {
        size_t h = 1469598103934665603ull;
        for (u32 x : v) {
            h ^= (size_t)x + 0x9e3779b97f4a7c15ull + (h << 6) + (h >> 2);
        }
        return h;
    }
};

static std::unordered_set<u32> compute_bad_targets(int R, int budget) {
    const u64 mask = mask_u64(R);
    std::unordered_set<u32> bad;
    ResidueSet start = close_under_doubling({0u, 1u}, R);

    std::queue<std::pair<ResidueSet, int>> q;
    q.push({start, 0});
    std::unordered_set<ResidueSet, VecHash> visited;
    visited.reserve(1 << 16);
    visited.insert(start);

    while (!q.empty()) {
        auto [S, d] = std::move(q.front());
        q.pop();

        if (has_two_dirs(S, R)) {
            for (u32 x : S) bad.insert(x);
        }

        if (d == budget) continue;
        for (size_t i = 0; i < S.size(); ++i) {
            for (size_t j = i + 1; j < S.size(); ++j) {
                u32 sum = (u32)(((u64)S[i] + (u64)S[j]) & mask);
                if (std::binary_search(S.begin(), S.end(), sum)) continue;

                ResidueSet T = S;
                T.push_back(sum);
                T = close_under_doubling(std::move(T), R);
                if (visited.insert(T).second) q.push({std::move(T), d + 1});
            }
        }
    }
    return bad;
}

static int ceil_log2_int(int n) {
    if (n <= 1) return 0;
    int p = 0;
    int v = 1;
    while (v < n) {
        v <<= 1;
        ++p;
    }
    return p;
}

static void gen_lower_clusters_rec(int weight, std::vector<u32>& out, int* gaps, int depth) {
    if (weight <= 0) return;
    if (weight == 1) {
        out.push_back(1);
        return;
    }
    if (depth == weight - 1) {
        u32 n = 1;
        for (int i = 0; i < weight - 1; ++i) n = (n << (gaps[i] + 1)) | 1u;
        out.push_back(n);
        return;
    }
    for (int g = 0; g < 4; ++g) {
        gaps[depth] = g;
        gen_lower_clusters_rec(weight, out, gaps, depth + 1);
    }
}

int main() {
    const std::string out_path = "scholz/cert/interface_charge_cert.csv";
    std::ofstream fout(out_path);
    if (!fout) {
        std::cerr << "Could not open " << out_path << " for writing.\n";
        return 2;
    }

    fout << "wtB,B,shift_t,R,wA,cA,required_f,Kres,status\n";

    std::vector<u32> Bs;
    Bs.reserve(2000);
    int gaps[6];
    for (int wtB = 1; wtB <= 6; ++wtB) {
        gen_lower_clusters_rec(wtB, Bs, gaps, 0);
    }
    Bs = normalize(std::move(Bs));

    bool all_ok = true;

    struct Key {
        int R = 0;
        int budget = 0;
        bool operator==(const Key& o) const { return R == o.R && budget == o.budget; }
    };
    struct KeyHash {
        size_t operator()(const Key& k) const noexcept {
            return std::hash<u64>{}(((u64)k.R << 32) ^ (u64)k.budget);
        }
    };
    std::unordered_map<Key, std::unordered_set<u32>, KeyHash> bad_cache;
    bad_cache.reserve(128);

    auto get_bad = [&](int R, int budget) -> const std::unordered_set<u32>& {
        Key k{R, budget};
        auto it = bad_cache.find(k);
        if (it != bad_cache.end()) return it->second;
        std::cerr << "[interface_charge_certify] computing bad-targets for R=" << R << " budget=" << budget << "\n";
        auto bad = compute_bad_targets(R, budget);
        auto [it2, _] = bad_cache.emplace(k, std::move(bad));
        return it2->second;
    };

    for (u32 B0 : Bs) {
        const int wtB_odd = popcount_u32(odd_part_mod_pow2(B0));
        const int wA = 7 - wtB_odd;
        const int cA = ceil_log2_int(wA);
        const int required_f = 2;
        const int Kres = 4 - cA - required_f;

        for (int t = 0; t <= 3; ++t) {
            u32 B = (t >= 32) ? 0u : (B0 << t);
            int r = std::max(1, bit_length_u32(B));
            int R = r + 4;
            if (R < 1) R = 1;
            if (R > 32) R = 32;

            std::string status = "PASS";
            if (Kres >= 0) {
                const auto& bad = get_bad(R, Kres);
                bool ok = (bad.find(B) == bad.end());
                if (!ok) {
                    status = "FAIL";
                    all_ok = false;
                    std::cerr << "FAIL: wtB_odd=" << wtB_odd << " B0=" << B0 << " t=" << t << " B=" << B
                              << " r=" << r << " R=" << R << " wA=" << wA << " cA=" << cA << " Kres=" << Kres << "\n";
                }
            } else {
                status = "PASS(Kres<0)";
            }

            fout << wtB_odd << "," << B << "," << t << "," << R << "," << wA << "," << cA << "," << required_f << "," << Kres << "," << status
                 << "\n";
        }
    }

    fout.flush();
    std::cerr << "Wrote " << out_path << "\n";
    return all_ok ? 0 : 1;
}
