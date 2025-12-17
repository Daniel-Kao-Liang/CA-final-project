#include <iostream>
#include <string>
#include <fstream>
#include <cmath>
#include <map>
#include <algorithm>
#include <vector>
#include <set>
#include <iomanip>

using namespace std;

// --- Helper for Robust Config Reading ---
// 解決 Address_bits 與 Address bits 格式不一致的問題
void read_config_value(ifstream& infile, int& value) {
    string label;
    infile >> label; 
    if (label.back() != ':') { 
        infile >> label; 
    }
    infile >> value;
}

// --- Block Class ---
class Block{
public:
    bool valid;
    int ref_bit; // 1: Recently Used, 0: Replacement Candidate
    string tag;

    Block() {
        valid = false;
        ref_bit = 0; // Init to 0
        tag = "";
    }
};

// --- Cache Class (Fixed Clock Policy) ---
class Cache{
public:
    Cache() {
        associativity = 0;
        clock_hand = 0;
    }
    Cache(int associativity) {
        this->associativity = associativity;
        this->clock_hand = 0;
        if(associativity > 0) blocks.resize(associativity);
    }
    
    // 回傳 true 代表 Hit, false 代表 Miss
    bool access(string tag) {
        // 1. Check Hit
        for(int i = 0; i < associativity; i++) {
            if(blocks[i].valid && blocks[i].tag == tag) {
                blocks[i].ref_bit = 1; // Hit: set ref to 1
                return true; // Hit: No pointer movement
            }
        }

        // 2. Handle Miss (Clock Replacement)
        while(true) {
            // Case A: Empty Slot
            if(!blocks[clock_hand].valid) {
                blocks[clock_hand].valid = true;
                blocks[clock_hand].tag = tag;
                blocks[clock_hand].ref_bit = 1; // Insert: set to 1
                clock_hand = (clock_hand + 1) % associativity; // Move next
                return false;
            }

            // Case B: Second Chance
            if(blocks[clock_hand].ref_bit == 1) {
                blocks[clock_hand].ref_bit = 0; // Reset to 0
                clock_hand = (clock_hand + 1) % associativity; // Move next
            } 
            // Case C: Victim Found
            else {
                // Replace
                blocks[clock_hand].tag = tag;
                blocks[clock_hand].ref_bit = 1; // New block gets 1
                clock_hand = (clock_hand + 1) % associativity; // Move next
                return false;
            }
        }
    }

private:
    vector<Block> blocks;
    int associativity;
    int clock_hand; // [Critical Fix] Clock Policy needs a pointer
};

// --- Simulate Class (Preserved Your Logic) ---
class Simulate{
public:
    Simulate() {
        address_bits = 0; block_size = 0; cache_sets = 0; associativity = 0;
        correlation_matrix = NULL;
    }
    ~Simulate() {
        if(correlation_matrix != NULL) delete correlation_matrix;
    }
    void set_address_bits(int v) { address_bits = v;}
    void set_block_size(int v) {
        block_size = v;
        offset_bits = (int)log2(block_size);
        cache_bits = address_bits - offset_bits;
    }
    void set_cache_sets(int v) {
        cache_sets = v;
        index_bits_num = (int)log2(cache_sets);
    }
    void set_associativity(int v) { associativity = v;}
    
    // Getters
    int get_address_bits() { return address_bits;}
    int get_block_size() { return block_size;}
    int get_cache_sets() { return cache_sets;}
    int get_associativity() { return associativity;}
    int get_offset_bits() { return offset_bits; }
    int get_index_bits_num() { return index_bits_num; }

    void simulation();
    void recursion(int bits_remain, vector<double> quality, set<int> current_idx);
    void initialize();
    void set_corr_matrix(ifstream& infile);
    void output(ofstream& outfile);

private:
    int address_bits;
    int block_size;
    int cache_sets;
    int associativity;
    int offset_bits;
    int index_bits_num;
    int cache_bits; // M - O
    
    vector<string> ref; // Stores reversed binary strings
    set< set < int > > candidates; // Changed name from index_bits to avoid confusion
    string benchmark;
    
    map<string,Cache> simulate_cache; 
    vector< vector< double > >* correlation_matrix;
    int min_miss = 1e9;
    set < int > best_bits;
};

// --- Main Function ---
int main(int argc, char *argv[]) {
    if(argc != 4) {
        cout << "Usage: ./project <cache.org> <reference.lst> <index.rpt>" << endl;
        return 1;
    }

    ifstream infile;
    ofstream outfile;
    Simulate simulate;

    // 1. Read Config (Robust)
    infile.open(argv[1]);
    if (!infile) { cout << "Error opening " << argv[1] << endl; return 1; }
    
    int val;
    read_config_value(infile, val); simulate.set_address_bits(val);
    read_config_value(infile, val); simulate.set_block_size(val);
    read_config_value(infile, val); simulate.set_cache_sets(val);
    read_config_value(infile, val); simulate.set_associativity(val);
    infile.close();

    // 2. Read Trace & Setup
    infile.open(argv[2]);
    if (!infile) { cout << "Error opening " << argv[2] << endl; return 1; }
    simulate.set_corr_matrix(infile); // Reads trace inside
    infile.close();

    // 3. Initialize & Run Optimization
    simulate.initialize(); // Runs recursion
    simulate.simulation(); // Runs cache sim on candidates

    // 4. Output Results
    outfile.open(argv[3]);
    if (!outfile) { cout << "Error opening " << argv[3] << endl; return 1; }
    simulate.output(outfile);
    outfile.close();

    return 0;
}

// --- Implementation of Simulate ---

void Simulate::set_corr_matrix(ifstream& infile) {
    string line;
    // Skip to benchmark name or just find it
    // Handle user's read logic, assuming .benchmark is present
    while(infile >> line) {
        if(line == ".benchmark") {
            infile >> benchmark;
            break;
        }
    }
    
    // Read references
    while(infile >> line) {
        if(line == ".end") break;
        // Reverse so index 0 is LSB (matches your logic)
        reverse(line.begin(), line.end());
        ref.push_back(line);
    }

    // Build Correlation Matrix
    correlation_matrix = new vector< vector<double> >(cache_bits, vector<double>(cache_bits, 0));
    int N = ref.size();
    
    for(int i = 0; i < N; i++) {
        for(int j = 0; j < cache_bits; j++) {
            for(int k = j + 1; k < cache_bits; k++) {
                // j and k are bit positions relative to offset
                if(ref[i][j+offset_bits] == ref[i][k+offset_bits]) {
                    (*correlation_matrix)[j][k]++;
                    (*correlation_matrix)[k][j]++;
                }
            }
        }
    }
    
    // Normalize (Eq 7 in paper)
    for(int i = 0; i < cache_bits; i++) {
        for(int j = 0; j < cache_bits; j++) {
            double E = (*correlation_matrix)[i][j];
            double D = N - E;
            double val = 0;
            if(max(E, D) > 0) val = min(E, D) / max(E, D);
            (*correlation_matrix)[i][j] = val;
        }
    }
}

void Simulate::initialize() {
    // Calc Quality (Eq 5 in paper)
    vector<double> quality(cache_bits, 0);
    int N = ref.size();

    for(int i = 0; i < N; i++) {
        for(int j = 0; j < cache_bits; j++) {
            if(ref[i][j+offset_bits] == '1') quality[j]++;
        }
    }

    for(int i = 0; i < cache_bits; i++) {
        double ones = quality[i];
        double zeros = N - ones;
        if(max(ones, zeros) > 0)
            quality[i] = min(ones, zeros) / max(ones, zeros);
        else 
            quality[i] = 0;
    }

    // Start Recursion to find candidate bits
    recursion(index_bits_num, quality, set<int>());
    
    // Fallback: If recursion didn't produce candidates (e.g. strict filtering), add LSB
    if(candidates.empty()) {
        set<int> lsb_set;
        for(int i=0; i<index_bits_num; i++) lsb_set.insert(i);
        candidates.insert(lsb_set);
    }
}

void Simulate::recursion(int bits_remain, vector<double> quality, set<int> current_idx){
    if(bits_remain == 0) {
        candidates.insert(current_idx);
        return;
    }

    // Find max quality
    double max_q = -1.0;
    for(int j = 0; j < cache_bits; j++) {
        if(quality[j] > max_q) max_q = quality[j];
    }

    // Select bits close to max quality
    bool found = false;
    for(int j = 0; j < cache_bits; j++) {
        // Your logic: if quality is high enough relative to max (or absolute threshold)
        // Adjusted: Pick the absolute max to behave like Greedy (Algorithm 1)
        // Or keep your logic if it's working:
        if(quality[j] >= max_q - 0.001 && quality[j] > 0) { // Simple float comparison
             found = true;
             set<int> next_idx = current_idx;
             next_idx.insert(j);
             
             // Update quality with correlation
             vector<double> next_quality = quality;
             next_quality[j] = -1.0; // Mark as used
             
             for(int k = 0; k < cache_bits; k++) {
                 if(next_quality[k] != -1.0) {
                     next_quality[k] *= (*correlation_matrix)[j][k];
                 }
             }
             recursion(bits_remain - 1, next_quality, next_idx);
             // Break after finding the best to avoid explosion (Greedy approach)
             break; 
        }
    }
    
    // If all qualities are 0 or none found, pick first available
    if(!found) {
         for(int j=0; j<cache_bits; j++) {
             if(current_idx.find(j) == current_idx.end()) {
                 set<int> next_idx = current_idx;
                 next_idx.insert(j);
                 recursion(bits_remain-1, quality, next_idx); // Pass same quality, just pick one
                 break;
             }
         }
    }
}

void Simulate::simulation() {
    min_miss = 1e9;

    for(auto it = candidates.begin(); it != candidates.end(); ++it) {
        int tmp_miss = 0;
        simulate_cache.clear(); // Reset cache for new candidate

        for(int i=0; i<ref.size(); i++){
            // Construct Index and Tag
            // Note: ref[i] is REVERSED. offset_bits is from index 0.
            // j iterates 0 to cache_bits-1. Actual bit index is j + offset_bits.
            
            string index_str = "";
            string tag_str = "";
            
            for(int j = 0; j < cache_bits; j++) {
                char bit = ref[i][j + offset_bits];
                if(it->find(j) != it->end()) {
                    index_str += bit;
                } else {
                    tag_str += bit;
                }
            }
            
            // Access Cache
            if(simulate_cache.find(index_str) == simulate_cache.end()) {
                simulate_cache[index_str] = Cache(associativity);
            }
            
            if(!simulate_cache[index_str].access(tag_str)) {
                tmp_miss++;
            }
        }
        
        if(tmp_miss < min_miss) {
            min_miss = tmp_miss;
            best_bits = *it;
        }
    }
}

void Simulate::output(ofstream& outfile) {
    outfile << "Address bits: " << address_bits << endl;
    outfile << "Block size: " << block_size << endl;
    outfile << "Cache sets: " << cache_sets << endl;
    outfile << "Associativity: " << associativity << endl;
    outfile << "Offset bit count: " << offset_bits << endl;
    outfile << "Indexing bit count: " << index_bits_num << endl;
    
    outfile << "Indexing bits:";
    // Output from MSB to LSB (standard format)
    // best_bits contains indices relative to offset. 
    // Actual bit index = val + offset_bits.
    vector<int> out_bits;
    for(int b : best_bits) out_bits.push_back(b + offset_bits);
    sort(out_bits.rbegin(), out_bits.rend()); // Sort descending
    
    for(int b : out_bits) outfile << " " << b;
    outfile << endl;
    
    outfile << ".benchmark " << benchmark << endl;

    // Run one last time with best bits to print details
    simulate_cache.clear();
    for(int i=0; i<ref.size(); i++){
        string index_str = "";
        string tag_str = "";
        for(int j = 0; j < cache_bits; j++) {
            char bit = ref[i][j + offset_bits];
            if(best_bits.find(j) != best_bits.end()) {
                index_str += bit;
            } else {
                tag_str += bit;
            }
        }
        
        if(simulate_cache.find(index_str) == simulate_cache.end()) {
            simulate_cache[index_str] = Cache(associativity);
        }
        
        bool hit = simulate_cache[index_str].access(tag_str);
        
        // Print original string (need to reverse back)
        string orig = ref[i];
        reverse(orig.begin(), orig.end());
        
        outfile << orig << " " << (hit ? "hit" : "miss") << endl;
    }
    
    outfile << ".end" << endl;
    outfile << "Total cache miss count: " << min_miss << endl;
}
