#include <fstream>
#include <iostream>
#include <limits>
#include <sstream>
#include <string>
#include <unordered_map>
#include <vector>
#include <algorithm>

using u64  = std::uint64_t;

#define main chain_solver_cli_main
#include "chain_solver.cpp"
#undef main

static std::string u64_to_dec(u64 x) {
    return std::to_string(x);
}

struct FloorInfo {
    int s_floor;
    int k_stop;
};

static std::unordered_map<u64, FloorInfo> read_floor_csv(const std::string& path) {
    std::unordered_map<u64, FloorInfo> mp;
    std::ifstream in(path);
    if (!in) {
        std::cerr << "ERROR: cannot open floor file: " << path << "\n";
        std::exit(1);
    }
    std::string line;
    bool first = true;
    while (std::getline(in, line)) {
        if (line.empty()) continue;
        if (first) {
            first = false;
            if (line.find("m") != std::string::npos && line.find("s_floor") != std::string::npos) {
                continue;
            }
        }
        std::stringstream ss(line);
        std::string a,b,c;
        if (!std::getline(ss, a, ',')) continue;
        if (!std::getline(ss, b, ',')) continue;
        if (!std::getline(ss, c, ',')) continue;
        u64 m = (u64)std::stoull(a);
        int s_floor = std::stoi(b);
        int k_stop  = std::stoi(c);
        mp[m] = FloorInfo{s_floor, k_stop};
    }
    return mp;
}

static std::vector<u64> enumerate_compact_kernels() {
    std::vector<u64> out;
    out.reserve(4096);
    for (u64 mask = 0; mask < (1ull << 12); ++mask) {
        int gaps[6];
        for (int i=0;i<6;i++) gaps[i] = (mask >> (2*i)) & 3;
        int p[7]; p[0]=0;
        for (int i=0;i<6;i++) p[i+1] = p[i] + gaps[i] + 1;
        u64 m=0;
        for (int i=0;i<7;i++) m |= (1ull << p[i]);
        out.push_back(m);
    }
    std::sort(out.begin(), out.end());
    out.erase(std::unique(out.begin(), out.end()), out.end());
    return out;
}

int main(int argc, char** argv) {
    const std::string floor_path = (argc>=2) ? argv[1] : "trailing_floor_cert.csv";
    const std::string out_path   = (argc>=3) ? argv[2] : "trailing_plateau_cert.csv";

    auto floor = read_floor_csv(floor_path);
    auto kernels = enumerate_compact_kernels();

    std::ofstream out(out_path);
    if (!out) {
        std::cerr << "ERROR: cannot open output file: " << out_path << "\n";
        return 1;
    }
    out << "m,k,n_k,n_k1,s_k,s_k1,plateau,e_k1\n";

    for (u64 m : kernels) {
        auto it = floor.find(m);
        if (it == floor.end()) {
            std::cerr << "ERROR: missing floor info for kernel m=" << m << "\n";
            return 1;
        }
        int s_floor = it->second.s_floor;
        int k_stop  = it->second.k_stop;
        if (k_stop < 0) {
            std::cerr << "ERROR: k_stop<0 for kernel m=" << m << " (floor file invalid)\n";
            return 1;
        }

        if (k_stop == 0) continue;
        if (k_stop != 1) {
            std::cerr << "ERROR: unexpected k_stop=" << k_stop << " for m=" << m << "\n";
            return 1;
        }

        if (m > (std::numeric_limits<u64>::max() >> 1)) {
            std::cerr << "ERROR: u64 overflow at m=" << m << "\n";
            return 1;
        }
        u64 n0 = m;
        u64 n1 = m << 1;

        int lambda0 = bit_length(n0) - 1;
        int s0 = s_floor + 1;
        int s1 = s_floor;
        int L0 = lambda0 + s0;
        int L1 = L0;

        int Lstar1 = min_star_chain_length(n1);
        int e1 = Lstar1 - L1;
        if (e1 > 1) {
            std::cerr << "FAIL: plateau target with e>1 at m=" << m << "\n";
            std::cerr << "  n_1=" << u64_to_dec(n1) << " L=" << L1 << " L*=" << Lstar1 << " e=" << e1 << "\n";
            return 1;
        }

        out << m << ",0," << u64_to_dec(n0) << "," << u64_to_dec(n1)
            << "," << s0 << "," << s1 << ",1," << e1 << "\n";
    }

    out.close();
    std::cerr << "Wrote " << out_path << "\n";
    return 0;
}
