#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <cmath>
#include <algorithm>
#include <iomanip>
#include <set>
#include <map>
#include <bitset>

using namespace std;

struct Config {
    int address_bits;
    int block_size;
    int cache_sets;
    int associativity;
    int offset_bit_count;
    int index_bit_count;
};

struct CacheBlock {
    bool valid;
    unsigned long long tag;
    int ref_bit;
};

struct CacheSet {
    vector<CacheBlock> blocks;
    int clock_hand;
};

vector<unsigned long long> trace_data;
vector<string> trace_strings; 

int log2_int(int n) {
    int res = 0;
    while (n >>= 1) res++;
    return res;
}

unsigned long long binaryStringToInt(const string& s) {
    return stoull(s, nullptr, 2);
}

int get_set_index(unsigned long long addr, const vector<int>& index_bits) {
    int index = 0;
    for (int i = 0; i < index_bits.size(); ++i) {
        int bit = (addr >> index_bits[i]) & 1;
        index |= (bit << i);
    }
    return index;
}

vector<int> calculate_optimized_bits(const Config& cfg, const vector<unsigned long long>& trace) {
    vector<unsigned long long> unique_refs = trace;
    sort(unique_refs.begin(), unique_refs.end());
    unique_refs.erase(unique(unique_refs.begin(), unique_refs.end()), unique_refs.end());
    
    int N = unique_refs.size();
    int M = cfg.address_bits;
    int K = cfg.index_bit_count;
    int Offset = cfg.offset_bit_count;

    vector<int> candidate_bits;
    for (int i = Offset; i < M; ++i) {
        candidate_bits.push_back(i);
    }

    if (candidate_bits.size() <= K) {
        sort(candidate_bits.rbegin(), candidate_bits.rend());
        return candidate_bits;
    }

    map<int, double> Q;
    for (int bit : candidate_bits) {
        int zeros = 0;
        int ones = 0;
        for (unsigned long long addr : unique_refs) {
            if ((addr >> bit) & 1) ones++;
            else zeros++;
        }
        if (zeros == 0 || ones == 0) Q[bit] = 0.0;
        else Q[bit] = (double)min(zeros, ones) / (double)max(zeros, ones);
    }

    map<int, map<int, double>> C;
    for (size_t i = 0; i < candidate_bits.size(); ++i) {
        for (size_t j = i + 1; j < candidate_bits.size(); ++j) {
            int b1 = candidate_bits[i];
            int b2 = candidate_bits[j];
            int same = 0;
            int diff = 0;
            for (unsigned long long addr : unique_refs) {
                int v1 = (addr >> b1) & 1;
                int v2 = (addr >> b2) & 1;
                if (v1 == v2) same++;
                else diff++;
            }
            double val = 0.0;
            if (same == 0 || diff == 0) val = 0.0;
            else val = (double)min(same, diff) / (double)max(same, diff);
            
            C[b1][b2] = val;
            C[b2][b1] = val;
        }
        C[candidate_bits[i]][candidate_bits[i]] = 1.0;
    }

    vector<int> S;
    vector<int> remaining = candidate_bits;

    for (int k = 0; k < K; ++k) {
        int best_bit = -1;
        double max_q = -1.0;

        for (int bit : remaining) {
            if (Q[bit] > max_q) {
                max_q = Q[bit];
                best_bit = bit;
            }
        }

        if (best_bit != -1) {
            S.push_back(best_bit);
            remaining.erase(remove(remaining.begin(), remaining.end(), best_bit), remaining.end());

            for (int bit : remaining) {
                Q[bit] = Q[bit] * C[best_bit][bit];
            }

        } else {
            if (!remaining.empty()) {
                S.push_back(remaining.back());
                remaining.pop_back();
            }
        }
    }

    sort(S.rbegin(), S.rend());
    return S;
}


void simulate(Config cfg, const vector<int>& index_bits, string output_filename) {
    ofstream out(output_filename);
    
    out << "Address bits: " << cfg.address_bits << endl;
    out << "Block size: " << cfg.block_size << endl;
    out << "Cache sets: " << cfg.cache_sets << endl;
    out << "Associativity: " << cfg.associativity << endl;
    out << "Offset bit count: " << cfg.offset_bit_count << endl;
    out << "Indexing bit count: " << cfg.index_bit_count << endl;
    out << "Indexing bits:";
    for (int b : index_bits) {
        out << " " << b;
    }
    out << endl;
    out << ".benchmark testcase1" << endl;

    vector<CacheSet> cache(cfg.cache_sets);
    for (int i = 0; i < cfg.cache_sets; ++i) {
        cache[i].blocks.resize(cfg.associativity);
        cache[i].clock_hand = 0;
        for (int j = 0; j < cfg.associativity; ++j) {
            cache[i].blocks[j].valid = false;
            cache[i].blocks[j].ref_bit = 0;
        }
    }

    int total_miss = 0;

    for (size_t t = 0; t < trace_data.size(); ++t) {
        unsigned long long addr = trace_data[t];
        string addr_str = trace_strings[t];

        int set_idx = get_set_index(addr, index_bits);
        CacheSet& set = cache[set_idx];

        unsigned long long tag = addr >> cfg.offset_bit_count;

        bool hit = false;
        int hit_way = -1;

        for (int i = 0; i < cfg.associativity; ++i) {
            if (set.blocks[i].valid && set.blocks[i].tag == tag) {
                hit = true;
                hit_way = i;
                break;
            }
        }

        if (hit) {
            out << addr_str << " hit" << endl;
            set.blocks[hit_way].ref_bit = 1;
        } else {
            out << addr_str << " miss" << endl;
            total_miss++;

            int ways = cfg.associativity;
            int& hand = set.clock_hand;
            int victim_way = -1;

            while (true) {
                if (!set.blocks[hand].valid) {
                    victim_way = hand;
                    hand = (hand + 1) % ways;
                    break;
                }

                if (set.blocks[hand].ref_bit == 1) {
                    set.blocks[hand].ref_bit = 0;
                    hand = (hand + 1) % ways;
                } else {
                    victim_way = hand;
                    hand = (hand + 1) % ways;
                    break;
                }
            }

            set.blocks[victim_way].valid = true;
            set.blocks[victim_way].tag = tag;
            set.blocks[victim_way].ref_bit = 1;
        }
    }

    out << ".end" << endl;
    out << "Total cache miss count: " << total_miss << endl;
    out.close();
}

int main(int argc, char* argv[]) {
    if (argc != 4) {
        cerr << "Usage: ./project <cache.org> <reference.lst> <index.rpt>" << endl;
        return 1;
    }

    string org_file = argv[1];
    string lst_file = argv[2];
    string rpt_file = argv[3];

    Config cfg;

    ifstream fin_org(org_file);
    if (!fin_org) { cerr << "Error opening " << org_file << endl; return 1; }
    
    string dummy;
    fin_org >> dummy >> dummy >> cfg.address_bits;
    fin_org >> dummy >> cfg.block_size;
    fin_org >> dummy >> cfg.cache_sets;
    fin_org >> dummy >> cfg.associativity;
    fin_org.close();

    cfg.offset_bit_count = log2_int(cfg.block_size);
    cfg.index_bit_count = log2_int(cfg.cache_sets);

    ifstream fin_lst(lst_file);
    if (!fin_lst) { cerr << "Error opening " << lst_file << endl; return 1; }
    
    string line;
    bool reading_benchmark = false;
    while (fin_lst >> line) {
        if (line == ".benchmark") {
            fin_lst >> line;
            reading_benchmark = true;
            continue;
        }
        if (line == ".end") break;
        if (reading_benchmark) {
            trace_strings.push_back(line);
            trace_data.push_back(binaryStringToInt(line));
        }
    }
    fin_lst.close();

    vector<int> index_bits;
    
    //bonus flag
    bool use_optimization = false;

    if (!use_optimization) {
        for (int i = 0; i < cfg.index_bit_count; ++i) {
            index_bits.push_back(cfg.offset_bit_count + i);
        }
        sort(index_bits.rbegin(), index_bits.rend());
    } else {
        index_bits = calculate_optimized_bits(cfg, trace_data);
    }

    simulate(cfg, index_bits, rpt_file);

    return 0;
}