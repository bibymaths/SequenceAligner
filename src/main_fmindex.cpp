#include <array>
#include <cstring>
#include <iostream>
#include <chrono>
#include <iomanip>
#include <fstream>
#include <vector>
#include <string>
#include <algorithm>
#include <map>
#include <numeric>
#include <stdexcept>
#include <omp.h>
#include <immintrin.h>
#include <mpi.h>
#include <sstream>
#include <unordered_map>
#include <filesystem>
#include <climits>
#include <memory>


#ifndef EDNAFULL_MATRIX_DEFINED
#define EDNAFULL_MATRIX_DEFINED
const int EDNAFULL_SIZE = 15;
int EDNAFULL_matrix[EDNAFULL_SIZE][EDNAFULL_SIZE] = {
// A  C  G  T  R  Y  S  W  K  M  B  D  H  V  N (X)
  {5,-4,-4,-4, 1,-4, 1, 1,-4, 1,-4, 1, 1, 1,-2}, // A
 {-4, 5,-4,-4,-4, 1, 1,-4, 1,-4, 1, 1,-4, 1,-2}, // C
 {-4,-4, 5,-4, 1, 1,-4,-4, 1,-4, 1,-4, 1, 1,-2}, // G
 {-4,-4,-4, 5,-4, 1,-4, 1, 1,-4, 1, 1,-4, 1,-2}, // T (U)
  {1,-4, 1,-4,-1,-4,-2,-2,-2,-2,-3,-2,-2,-2,-1}, // R (A|G)
 {-4, 1, 1, 1,-4,-1,-2,-2,-2,-2,-2,-3,-2,-2,-1}, // Y (C|T)
  {1, 1,-4,-4,-2,-2,-1,-4,-2,-4,-2,-2,-2,-2,-1}, // S (C|G)
  {1,-4,-4, 1,-2,-2,-4,-1,-4,-2,-2,-2,-2,-2,-1}, // W (A|T)
 {-4, 1, 1, 1,-2,-2,-2,-4,-1,-4,-2,-2,-2,-2,-1}, // K (G|T)
  {1,-4, 1,-4,-2,-4,-2,-2,-4,-1,-2,-2,-2,-2,-1}, // M (A|C)
 {-4, 1, 1, 1,-3,-2,-2,-2,-2,-2,-1,-2,-3,-3,-1}, // B (C|G|T)
  {1, 1,-4, 1,-2,-3,-2,-2,-2,-2,-2,-1,-3,-3,-1}, // D (A|G|T)
  {1, 1, 1,-4,-2,-2,-2,-2,-2,-2,-3,-3,-1,-3,-1}, // H (A|C|T)
  {1, 1, 1, 1,-2,-2,-2,-2,-2,-2,-3,-3,-3,-1,-1}, // V (A|C|G)
 {-2,-2,-2,-2,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1}  // N (X) (A|C|G|T)
};
#endif

#ifndef EBLOSUM62_MATRIX_DEFINED
#define EBLOSUM62_MATRIX_DEFINED
const int EBLOSUM62_SIZE = 24; // A..V, B, Z, X, *
int EBLOSUM62_matrix[EBLOSUM62_SIZE][EBLOSUM62_SIZE] = {
// A  R  N  D  C  Q  E  G  H  I  L  K  M  F  P  S  T  W  Y  V  B  Z  X  *
  {4,-1,-2,-2, 0,-1,-1, 0,-2,-1,-1,-1,-1,-2,-1, 1, 0,-3,-2, 0,-2,-1,-0,-4}, // A
 {-1, 5, 0,-2,-3, 1, 0,-2, 0,-3,-2, 2,-1,-3,-2,-1,-1,-3,-2,-3,-1, 0,-1,-4}, // R
 {-2, 0, 6, 1,-3, 0, 0, 0, 1,-3,-3, 0,-2,-3,-2, 1, 0,-4,-2,-3, 3, 0,-1,-4}, // N
 {-2,-2, 1, 6,-3, 0, 2,-1,-1,-3,-4,-1,-3,-3,-1, 0,-1,-4,-3,-3, 4, 1,-1,-4}, // D
  {0,-3,-3,-3, 9,-3,-4,-3,-3,-1,-1,-3,-1,-2,-3,-1,-1,-2,-2,-1,-3,-3,-2,-4}, // C
 {-1, 1, 0, 0,-3, 5, 2,-2, 0,-3,-2, 1, 0,-3,-1, 0,-1,-2,-1,-2, 0, 3,-1,-4}, // Q
 {-1, 0, 0, 2,-4, 2, 5,-2, 0,-3,-3, 1,-2,-3,-1, 0,-1,-3,-2,-2, 1, 4,-1,-4}, // E
  {0,-2, 0,-1,-3,-2,-2, 6,-2,-4,-4,-2,-3,-3,-2, 0,-2,-2,-3,-3,-1,-2,-1,-4}, // G
 {-2, 0, 1,-1,-3, 0, 0,-2, 8,-3,-3,-1,-2,-1,-2,-1,-2,-2, 2,-3, 0, 0,-1,-4}, // H
 {-1,-3,-3,-3,-1,-3,-3,-4,-3, 4, 2,-3, 1, 0,-3,-2,-1,-3,-1, 3,-3,-3,-1,-4}, // I
 {-1,-2,-3,-4,-1,-2,-3,-4,-3, 2, 4,-2, 2, 0,-3,-2,-1,-2,-1, 1,-4,-3,-1,-4}, // L
 {-1, 2, 0,-1,-3, 1, 1,-2,-1,-3,-2, 5,-1,-3,-1, 0,-1,-3,-2,-2, 0, 1,-1,-4}, // K
 {-1,-1,-2,-3,-1, 0,-2,-3,-2, 1, 2,-1, 5, 0,-2,-1,-1,-1,-1, 1,-3,-1,-1,-4}, // M
 {-2,-3,-3,-3,-2,-3,-3,-3,-1, 0, 0,-3, 0, 6,-4,-2,-2, 1, 3,-1,-3,-3,-1,-4}, // F
 {-1,-2,-2,-1,-3,-1,-1,-2,-2,-3,-3,-1,-2,-4, 7,-1,-1,-4,-3,-2,-2,-1,-2,-4}, // P
  {1,-1, 1, 0,-1, 0, 0, 0,-1,-2,-2, 0,-1,-2,-1, 4, 1,-3,-2,-2, 0, 0, 0,-4}, // S
  {0,-1, 0,-1,-1,-1,-1,-2,-2,-1,-1,-1,-1,-2,-1, 1, 5,-2,-2, 0,-1,-1, 0,-4}, // T
 {-3,-3,-4,-4,-2,-2,-3,-2,-2,-3,-2,-3,-1, 1,-4,-3,-2,11, 2,-3,-4,-3,-2,-4}, // W
 {-2,-2,-2,-3,-2,-1,-2,-3, 2,-1,-1,-2,-1, 3,-3,-2,-2, 2, 7,-1,-3,-2,-1,-4}, // Y
  {0,-3,-3,-3,-1,-2,-2,-3,-3, 3, 1,-2, 1,-1,-2,-2, 0,-3,-1, 4,-3,-2,-1,-4}, // V
 {-2,-1, 3, 4,-3, 0, 1,-1, 0,-3,-4, 0,-3,-3,-2, 0,-1,-4,-3,-3, 4, 1,-1,-4}, // B (N|D)
 {-1, 0, 0, 1,-3, 3, 4,-2, 0,-3,-3, 1,-1,-3,-1, 0,-1,-3,-2,-2, 1, 4,-1,-4}, // Z (Q|E)
  {0,-1,-1,-1,-2,-1,-1,-1,-1,-1,-1,-1,-1,-1,-2, 0, 0,-2,-1,-1,-1,-1,-1,-4}, // X
 {-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4, 1}  // *
};
#endif

// Global Variables
bool verbose = false;
bool binary = false;
bool txt = false;
int rank;

enum ScoreMode { MODE_DNA, MODE_PROTEIN };
using ScoreFn = int(*)(char, char);

double GAP_OPEN   = -5.0;
double GAP_EXTEND = -1.0;
static const int LINE_WIDTH = 80;

#define RESET "\033[0m"
#define GREEN "\033[32m"
#define RED   "\033[31m"
#define CYAN  "\033[36m"

// Suffix Array Construction
// This function constructs the suffix array for a given string using a
// O(n log n) algorithm. It sorts the suffixes based on their ranks and
// iteratively refines the ranks until they stabilize.
std::vector<int> suffix_array_construction(const std::string& s) {
    int n = s.length();
    std::vector<int> sa_arr(n);
    if (n == 0) return sa_arr;
    std::iota(sa_arr.begin(), sa_arr.end(), 0);

    std::vector<int> rank_arr(n);
    for (int i = 0; i < n; ++i) {
        rank_arr[i] = static_cast<unsigned char>(s[i]);
    }

    std::vector<int> tmp_rank_arr(n);
    for (int k_iter = 1; ; k_iter <<= 1) {
        std::sort(sa_arr.begin(), sa_arr.end(), [&](int a, int b) {
            if (rank_arr[a] != rank_arr[b]) {
                return rank_arr[a] < rank_arr[b];
            }
            int rank_a_k = (a + k_iter < n) ? rank_arr[a + k_iter] : -1;
            int rank_b_k = (b + k_iter < n) ? rank_arr[b + k_iter] : -1;
            return rank_a_k < rank_b_k;
        });

        tmp_rank_arr[sa_arr[0]] = 0;
        for (int i = 1; i < n; ++i) {
            int prev = sa_arr[i - 1];
            int curr = sa_arr[i];
            bool new_rank_val = (rank_arr[curr] != rank_arr[prev]);
            if (!new_rank_val) {
                 int rank_prev_k = (prev + k_iter < n) ? rank_arr[prev + k_iter] : -1;
                 int rank_curr_k = (curr + k_iter < n) ? rank_arr[curr + k_iter] : -1;
                 if (rank_curr_k != rank_prev_k) {
                     new_rank_val = true;
                 }
            }
            tmp_rank_arr[curr] = tmp_rank_arr[prev] + (new_rank_val ? 1 : 0);
        }
        rank_arr.swap(tmp_rank_arr);

        if (rank_arr[sa_arr[n - 1]] == n - 1) {
            break;
        }
        if (k_iter > n && n > 0) break; 
    }
    return sa_arr;
}

// FM Index Class
// This class implements the FM-index data structure for fast substring search.
class FMIndex {
public:
    std::string text_with_sentinel;
    std::vector<int> sa; 
    std::string bwt;     
    char sentinel_char;
    std::map<char, int> C; 
    std::map<char, std::vector<int>> Occ; 

    FMIndex(const std::string& text, char sentinel = '$')
        : sentinel_char(sentinel) {
        if (text.empty()) { 
             text_with_sentinel = std::string(1, sentinel_char);
        } else {
            text_with_sentinel = text + sentinel_char;
        }
        this->sa = suffix_array_construction(text_with_sentinel);
        this->bwt.reserve(text_with_sentinel.length());
        if (!this->sa.empty()){ 
            for (size_t i = 0; i < this->sa.size(); ++i) {
                if (this->sa[i] == 0) {
                    this->bwt += text_with_sentinel.back();
                } else {
                    this->bwt += text_with_sentinel[this->sa[i] - 1];
                }
            }
        } else if (text_with_sentinel.length() == 1 && text_with_sentinel[0] == sentinel_char) { 
            this->bwt = text_with_sentinel;
        }
        build_c_table();
        build_occ_table();
    }

    FMIndex() : sentinel_char('$') {} 

    void build_c_table() {
        std::map<char, int> counts;
        for (char ch : this->bwt) { counts[ch]++; }
        this->C.clear(); int total = 0;
        for (const auto& pair : counts) { this->C[pair.first] = total; total += pair.second; }
    }

    void build_occ_table() {
        this->Occ.clear();
        for (const auto& pair : this->C) { this->Occ[pair.first].push_back(0); }
        for (char ch_in_bwt : this->bwt) { if (this->Occ.find(ch_in_bwt) == this->Occ.end()) { this->Occ[ch_in_bwt].push_back(0);}}
        for (size_t i = 0; i < this->bwt.length(); ++i) {
            char current_char_in_bwt = this->bwt[i];
            for (auto& pair : this->Occ) {
                char c = pair.first; int prev_count = pair.second.back();
                pair.second.push_back(prev_count + (c == current_char_in_bwt ? 1 : 0));
            }
        }
    }

    std::pair<int, int> backward_search(const std::string& pattern) const {
        if (bwt.empty() || pattern.empty()) return {0,0};
        int l = 0; int r = bwt.length();
        for (int i = pattern.length() - 1; i >= 0; --i) {
            char ch = pattern[i];
            auto c_it = C.find(ch); auto occ_it = Occ.find(ch);
            if (c_it == C.end() || occ_it == Occ.end() || occ_it->second.empty()) return {0, 0};
            if (static_cast<size_t>(l) >= occ_it->second.size() || static_cast<size_t>(r) >= occ_it->second.size()) return {0,0};
            l = c_it->second + occ_it->second[l]; r = c_it->second + occ_it->second[r];
            if (l >= r) return {0, 0};
        } return {l, r};
    }

    std::vector<int> locate(const std::string& pattern) const {
        std::pair<int, int> sa_range = backward_search(pattern); std::vector<int> positions;
        if (sa_range.first < sa_range.second) {
            for (int i = sa_range.first; i < sa_range.second; ++i) {
                 if (static_cast<size_t>(i) < sa.size()) positions.push_back(sa[i]);
            }
        } std::sort(positions.begin(), positions.end()); return positions;
    }

    void save(std::ostream& os) const {
        size_t len = text_with_sentinel.length();
        os.write(reinterpret_cast<const char*>(&len), sizeof(len));
        if (len > 0) os.write(text_with_sentinel.data(), len);
        os.write(reinterpret_cast<const char*>(&sentinel_char), sizeof(sentinel_char));
        len = sa.size(); os.write(reinterpret_cast<const char*>(&len), sizeof(len));
        if (len > 0) os.write(reinterpret_cast<const char*>(sa.data()), len * sizeof(int));
        len = bwt.length(); os.write(reinterpret_cast<const char*>(&len), sizeof(len));
        if (len > 0) os.write(bwt.data(), len);
        len = C.size(); os.write(reinterpret_cast<const char*>(&len), sizeof(len));
        for (const auto& pair : C) { os.write(reinterpret_cast<const char*>(&pair.first), sizeof(char)); os.write(reinterpret_cast<const char*>(&pair.second), sizeof(int));}
        len = Occ.size(); os.write(reinterpret_cast<const char*>(&len), sizeof(len));
        for (const auto& pair : Occ) { os.write(reinterpret_cast<const char*>(&pair.first), sizeof(char)); size_t vec_len = pair.second.size(); os.write(reinterpret_cast<const char*>(&vec_len), sizeof(vec_len)); if (vec_len > 0) os.write(reinterpret_cast<const char*>(pair.second.data()), vec_len * sizeof(int));}
    }

    bool load(std::istream& is) {
        try {
            size_t len;
            is.read(reinterpret_cast<char*>(&len), sizeof(len)); if (!is || len > 2000000000) return false;
            text_with_sentinel.resize(len); if (len > 0) is.read(&text_with_sentinel[0], len);
            is.read(reinterpret_cast<char*>(&sentinel_char), sizeof(sentinel_char));
            is.read(reinterpret_cast<char*>(&len), sizeof(len)); if (!is || len > 2000000000) return false;
            sa.resize(len); if (len > 0) is.read(reinterpret_cast<char*>(sa.data()), len * sizeof(int));
            is.read(reinterpret_cast<char*>(&len), sizeof(len)); if (!is || len > 2000000000) return false;
            bwt.resize(len); if (len > 0) is.read(&bwt[0], len);
            is.read(reinterpret_cast<char*>(&len), sizeof(len)); if (!is) return false; C.clear();
            for (size_t i = 0; i < len; ++i) { char ch_c; int val_c; is.read(reinterpret_cast<char*>(&ch_c), sizeof(char)); is.read(reinterpret_cast<char*>(&val_c), sizeof(int)); if (!is) return false; C[ch_c] = val_c;}
            is.read(reinterpret_cast<char*>(&len), sizeof(len)); if (!is) return false; Occ.clear();
            for (size_t i = 0; i < len; ++i) { char ch_o; size_t vec_len; is.read(reinterpret_cast<char*>(&ch_o), sizeof(char)); is.read(reinterpret_cast<char*>(&vec_len), sizeof(vec_len)); if (!is || vec_len > 2000000000) return false; std::vector<int> occ_row(vec_len); if (vec_len > 0) is.read(reinterpret_cast<char*>(occ_row.data()), vec_len * sizeof(int)); if (!is && vec_len > 0) return false; Occ[ch_o] = occ_row;}
        } catch (const std::ios_base::failure& e) { if (verbose && ::rank == 0) std::cerr << "FMIndex load I/O exception: " << e.what() << std::endl; return false;}
          catch (const std::bad_alloc& e) { if (verbose && ::rank == 0) std::cerr << "FMIndex load memory exception: " << e.what() << std::endl; return false;}
        return is.good() && !is.eof();
    }
};

// -------- Scoring Lookups & Functions --------
static const std::array<uint8_t,256> char2idx = [](){
    std::array<uint8_t,256> m{}; m.fill(255);
    m[static_cast<unsigned char>('A')] =  0;  m[static_cast<unsigned char>('C')] =  1;  m[static_cast<unsigned char>('G')] =  2;  m[static_cast<unsigned char>('T')] =  3;  m[static_cast<unsigned char>('U')] =  3;
    m[static_cast<unsigned char>('R')] =  4;  m[static_cast<unsigned char>('Y')] =  5;  m[static_cast<unsigned char>('S')] =  6;  m[static_cast<unsigned char>('W')] =  7;  m[static_cast<unsigned char>('K')] =  8;
    m[static_cast<unsigned char>('M')] =  9;  m[static_cast<unsigned char>('B')] = 10;  m[static_cast<unsigned char>('D')] = 11;  m[static_cast<unsigned char>('H')] = 12;  m[static_cast<unsigned char>('V')] = 13;
    m[static_cast<unsigned char>('N')] = 14;  m[static_cast<unsigned char>('X')] = 14; return m;
}();
static const std::array<uint8_t,256> prot_idx = [](){
    std::array<uint8_t,256> m{}; m.fill(255);
    m[static_cast<unsigned char>('A')]=0; m[static_cast<unsigned char>('R')]=1; m[static_cast<unsigned char>('N')]=2; m[static_cast<unsigned char>('D')]=3; m[static_cast<unsigned char>('C')]=4; m[static_cast<unsigned char>('Q')]=5; m[static_cast<unsigned char>('E')]=6;
    m[static_cast<unsigned char>('G')]=7; m[static_cast<unsigned char>('H')]=8; m[static_cast<unsigned char>('I')]=9; m[static_cast<unsigned char>('L')]=10; m[static_cast<unsigned char>('K')]=11; m[static_cast<unsigned char>('M')]=12;
    m[static_cast<unsigned char>('F')]=13; m[static_cast<unsigned char>('P')]=14; m[static_cast<unsigned char>('S')]=15; m[static_cast<unsigned char>('T')]=16; m[static_cast<unsigned char>('W')]=17; m[static_cast<unsigned char>('Y')]=18;
    m[static_cast<unsigned char>('V')]=19; m[static_cast<unsigned char>('B')]=20; m[static_cast<unsigned char>('Z')]=21; m[static_cast<unsigned char>('X')]=22; m[static_cast<unsigned char>('*')]=23; return m;
}();

inline int edna_score(char x, char y) {
    uint8_t ix = char2idx[static_cast<uint8_t>(x)];
    uint8_t iy = char2idx[static_cast<uint8_t>(y)];
    if (ix == 255 || iy == 255) throw std::runtime_error(std::string("Invalid DNA code in edna_score: '") + x + "','" + y + "'");
    return EDNAFULL_matrix[ix][iy];
}
inline int blosum62_score(char x, char y) {
    uint8_t ix = prot_idx[static_cast<uint8_t>(x)];
    uint8_t iy = prot_idx[static_cast<uint8_t>(y)];
    if (ix == 255 || iy == 255) throw std::runtime_error(std::string("Invalid protein code in blosum62_score: '") + x + "','" + y + "'");
    return EBLOSUM62_matrix[ix][iy];
}
inline int score(char x, char y, ScoreMode mode) {
    if (mode == MODE_DNA) return edna_score(x, y);
    return blosum62_score(x, y);
}

// -------- Utilities (showProgressBar, FASTA I/O, Alignment Savers, Matrix Writers) --------
void showProgressBar(int progress, int total) {
    using namespace std::chrono;
    static auto start_time = steady_clock::now(); // static per-call, not per-rank here.
                                                 // If called by multiple ranks, this is an issue.
                                                 // Should be managed if progress is collective.
    if (::rank != 0) return; // Only rank 0 prints progress

    auto now = steady_clock::now();
    auto elapsed_s = duration_cast<seconds>(now - start_time).count();
    long eta_s = 0;
    if (progress > 0 && progress < total) {
        eta_s = elapsed_s * (total - progress) / progress;
    }

    auto format_hms_func = [](long secs_val) {
        std::ostringstream os_formatter;
        long h_val = secs_val / 3600;
        long m_val = (secs_val % 3600) / 60;
        long s_val = secs_val % 60;
        if (h_val > 0) os_formatter << h_val << "h";
        if (m_val > 0 || h_val > 0) os_formatter << std::setw(h_val > 0 ? 2:1) << std::setfill(h_val > 0 ? '0':' ') << m_val << "m";
        os_formatter << std::setw(2) << std::setfill('0') << s_val << "s";
        return os_formatter.str();
    };

    constexpr int bar_width = 50; // Adjusted for typical console
    float current_ratio = total > 0 ? static_cast<float>(progress) / total : 0.0f;
    int filled_pos = static_cast<int>(bar_width * current_ratio);

    std::cout << "\r[";
    for (int i_bar = 0; i_bar < bar_width; ++i_bar) {
        if (i_bar < filled_pos) std::cout << "=";
        else if (i_bar == filled_pos) std::cout << ">";
        else std::cout << " ";
    }
    std::cout << "] " << std::setw(3) << static_cast<int>(current_ratio * 100.0f) << "% "
              << progress << "/" << total
              << " | Elapsed: " << format_hms_func(elapsed_s)
              << " | ETA: " << format_hms_func(eta_s) << "   "; // Added spaces to clear previous longer lines
    std::cout << std::flush;
    if (progress == total) {
        std::cout << std::endl;
        start_time = steady_clock::now(); // Reset for next potential progress bar
    }
}

std::string getAccession(const std::string& header, ScoreMode mode) {
    if (mode == MODE_PROTEIN) { // Specific parsing for UniProt-like headers
        size_t firstPipe = header.find('|');
        if (firstPipe != std::string::npos) {
            size_t secondPipe = header.find('|', firstPipe + 1);
            if (secondPipe != std::string::npos) {
                return header.substr(firstPipe + 1, secondPipe - firstPipe - 1);
            }
        }
    }
    // Default: first word for DNA or if protein format not matched
    std::istringstream iss(header);
    std::string accession_val;
    iss >> accession_val; // Takes the first whitespace-separated token
    return accession_val;
}

std::string getGeneSymbol(const std::string& header, ScoreMode mode) {
    if (mode == MODE_DNA) {
        size_t open_p = header.find('(');
        size_t close_p = (open_p != std::string::npos) ? header.find(')', open_p + 1) : std::string::npos;
        if (open_p != std::string::npos && close_p != std::string::npos && close_p > open_p + 1) {
            return header.substr(open_p + 1, close_p - open_p - 1);
        }
    } else if (mode == MODE_PROTEIN) { // Try GN= field first for proteins
        size_t gn_pos = header.find("GN=");
        if (gn_pos != std::string::npos) {
            size_t start_gn = gn_pos + 3;
            size_t end_gn = header.find_first_of(" \t", start_gn); // Gene name ends at space or tab
            if (end_gn == std::string::npos) end_gn = header.length(); // Or end of header
            if (end_gn > start_gn) return header.substr(start_gn, end_gn - start_gn);
        }
        // Fallback for protein if GN= not found (e.g. after second pipe, before underscore)
        size_t firstPipe = header.find('|');
        size_t secondPipe = (firstPipe != std::string::npos) ? header.find('|', firstPipe + 1) : std::string::npos;
        if (secondPipe != std::string::npos) {
            size_t protein_name_start = secondPipe + 1;
            size_t underscore_pos = header.find('_', protein_name_start);
            if (underscore_pos != std::string::npos && underscore_pos > protein_name_start) { // Often GeneSymbol_SPECIES
                return header.substr(protein_name_start, underscore_pos - protein_name_start);
            }
            // If no underscore, might be the whole protein name part
            size_t first_space_after_pipes = header.find(' ', protein_name_start);
             if (first_space_after_pipes != std::string::npos) {
                 return header.substr(protein_name_start, first_space_after_pipes - protein_name_start);
             } else {
                 return header.substr(protein_name_start); // rest of string
             }
        }
    }
    return ""; // No gene symbol found
}

void processFasta(const std::string &filename, std::string &header_out, std::string &sequence_out) {
    std::ifstream file_in(filename);
    if (!file_in) {
        throw std::runtime_error("Error: Unable to open FASTA file " + filename);
    }
    std::string line_buffer;
    header_out = "";
    sequence_out = "";
    bool first_header_found = false;

    while (getline(file_in, line_buffer)) {
        if (line_buffer.empty()) continue;
        if (line_buffer[0] == '>') {
            if (!first_header_found) { // Only process the first sequence in the file
                header_out = line_buffer.substr(1);
                // Remove potential trailing carriage return for Windows/Mac line endings
                if (!header_out.empty() && header_out.back() == '\r') {
                    header_out.pop_back();
                }
                first_header_found = true;
            } else {
                break; // Stop after the first sequence
            }
        } else if (first_header_found) {
            // Remove potential trailing carriage return before appending
            if (!line_buffer.empty() && line_buffer.back() == '\r') {
                line_buffer.pop_back();
            }
            sequence_out += line_buffer;
        }
    }
    if (!first_header_found && filename != "-"){ // Check if file was empty or malformed
         // throw std::runtime_error("Error: No FASTA sequence found in " + filename);
         // Or just allow empty header/sequence if that's acceptable.
    }
}

void savePlainAlignment(const std::string &h1, const std::string &h2, const std::string &a1, const std::string &a2, std::ostream &os) {
    os << '>' << h1 << '\n';
    for (size_t i = 0; i < a1.length(); i += LINE_WIDTH) { os << a1.substr(i, LINE_WIDTH) << '\n'; }
    os << '>' << h2 << '\n';
    for (size_t i = 0; i < a2.length(); i += LINE_WIDTH) { os << a2.substr(i, LINE_WIDTH) << '\n'; }
}

void saveLCS(const std::string &id, const std::string &lcs_str_val, std::ostream &os) {
    os << '>' << id << "_LCS_len=" << lcs_str_val.size() << "\n";
    for (size_t i = 0; i < lcs_str_val.size(); i += LINE_WIDTH) {
        os << lcs_str_val.substr(i, LINE_WIDTH) << "\n";
    }
}

void printColoredAlignment(const std::string &seq1_aln, const std::string &seq2_aln, std::ostream &os = std::cout) {
    size_t aln_len = seq1_aln.length();
    if (aln_len == 0) { os << "No alignment to print.\n"; return; }
    if (seq1_aln.length() != seq2_aln.length()) { os << "Error: Aligned sequences have different lengths.\n"; return;}

    size_t pos1_count = 0, pos2_count = 0;
    for (size_t i = 0; i < aln_len; i += LINE_WIDTH) {
        size_t end_blk = std::min(i + LINE_WIDTH, aln_len);
        size_t blk_start_p1 = pos1_count + 1;
        size_t blk_start_p2 = pos2_count + 1;
        size_t current_blk_end_p1 = pos1_count;
        size_t current_blk_end_p2 = pos2_count;

        os << std::setw(6) << blk_start_p1 << " ";
        for (size_t j = i; j < end_blk; ++j) {
            if (seq1_aln[j] == seq2_aln[j]) os << GREEN << seq1_aln[j] << RESET;
            else if (seq1_aln[j] == '-' || seq2_aln[j] == '-') os << RED << seq1_aln[j] << RESET;
            else os << CYAN << seq1_aln[j] << RESET;
            if (seq1_aln[j] != '-') current_blk_end_p1++;
        }
        os << " " << current_blk_end_p1 << "\n";

        os << "       "; // Spacer for alignment bar or consensus
        for (size_t j = i; j < end_blk; ++j) {
            if (seq1_aln[j] == seq2_aln[j]) os << "|"; // Match
            else if (seq1_aln[j] == '-' || seq2_aln[j] == '-') os << " "; // Gap
            else os << "."; // Mismatch
        }
        os << "\n";

        os << std::setw(6) << blk_start_p2 << " ";
        for (size_t j = i; j < end_blk; ++j) {
            if (seq1_aln[j] == seq2_aln[j]) os << GREEN << seq2_aln[j] << RESET;
            else if (seq1_aln[j] == '-' || seq2_aln[j] == '-') os << RED << seq2_aln[j] << RESET;
            else os << CYAN << seq2_aln[j] << RESET;
            if (seq2_aln[j] != '-') current_blk_end_p2++;
        }
        os << " " << current_blk_end_p2 << "\n\n";
        pos1_count = current_blk_end_p1;
        pos2_count = current_blk_end_p2;
    }
}

void writeRawDPMatrix(const std::vector<std::vector<int>>& dp_matrix, const std::string& filename) {
    std::ofstream out_file(filename);
    if (!out_file) { std::cerr << "Error: Cannot write DP matrix to " << filename << "\n"; return; }
    for (const auto& row_vec : dp_matrix) {
        for (size_t j_col = 0; j_col < row_vec.size(); ++j_col) {
            out_file << std::setw(5) << row_vec[j_col];
            if (j_col != row_vec.size() - 1) out_file << " ";
        }
        out_file << "\n";
    }
    out_file.close();
}
void writeDPMatrix(const std::vector<std::vector<int>>& dp_matrix, const std::string& filename) {
    std::ofstream out_file(filename, std::ios::binary);
    if (!out_file) { std::cerr << "Error: Cannot write binary DP matrix to " << filename << "\n"; return; }
    int32_t num_rows = dp_matrix.empty() ? 0 : dp_matrix.size();
    int32_t num_cols = (num_rows > 0 && !dp_matrix[0].empty()) ? dp_matrix[0].size() : 0;
    out_file.write(reinterpret_cast<const char*>(&num_rows), sizeof(int32_t));
    out_file.write(reinterpret_cast<const char*>(&num_cols), sizeof(int32_t));
    if (num_rows > 0 && num_cols > 0) {
        for (const auto& row_vec : dp_matrix) {
            out_file.write(reinterpret_cast<const char*>(row_vec.data()), num_cols * sizeof(int32_t));
        }
    }
    out_file.close();
}
void writeRawCharMatrix(const std::vector<std::vector<char>>& char_matrix, const std::string& filename) {
    std::ofstream out_file(filename);
    if (!out_file) { std::cerr << "Error: Cannot open " << filename << "\n"; return; }
    size_t max_cols = 0; for (const auto& row_vec : char_matrix) max_cols = std::max(max_cols, row_vec.size());
    for (const auto& row_vec : char_matrix) {
        for (size_t i_col = 0; i_col < max_cols; ++i_col) {
            char ch_val = (i_col < row_vec.size()) ? row_vec[i_col] : ' ';
            out_file << ch_val;
            if (i_col + 1 < max_cols) out_file << ' ';
        }
        out_file << '\n';
    }
    out_file.close();
}
void writeCharMatrix(const std::vector<std::vector<char>>& char_matrix, const std::string& filename) {
    std::ofstream out_file(filename, std::ios::binary);
    if (!out_file) { std::cerr << "Error: Cannot open " << filename << " for writing.\n"; return; }
    int32_t num_rows = char_matrix.size();
    int32_t num_cols = 0; for (const auto& row_vec : char_matrix) num_cols = std::max<int32_t>(num_cols, static_cast<int32_t>(row_vec.size()));
    out_file.write(reinterpret_cast<const char*>(&num_rows), sizeof(int32_t));
    out_file.write(reinterpret_cast<const char*>(&num_cols), sizeof(int32_t));
    if (num_rows > 0 && num_cols > 0) {
        for (const auto& row_vec : char_matrix) {
            for (int32_t i_col = 0; i_col < num_cols; ++i_col) {
                char ch_val = (i_col < static_cast<int32_t>(row_vec.size())) ? row_vec[i_col] : ' ';
                out_file.write(&ch_val, sizeof(char));
            }
        }
    }
    out_file.close();
}

// -------- Original DP Helper Functions (for fallback full DP or if adapted) --------
struct AffineDPScores { int s_val = 0; int e_val = 0; int f_val = 0; char ptr = 'X'; }; // S, E(gap in X), F(gap in Y)

void initAffineDP(int n_len, std::vector<int>& prev_row_s, std::vector<int>& prev_row_e, // E is gap in X (horizontal)
                  std::vector<int>& prev_row_f, bool isGlobal) { // F is gap in Y (vertical)
    prev_row_s.assign(n_len + 1, isGlobal ? (INT_MIN / 2) : 0);
    prev_row_e.assign(n_len + 1, INT_MIN / 2); // E (gap in X)
    prev_row_f.assign(n_len + 1, INT_MIN / 2); // F (gap in Y)

    if (isGlobal) {
        prev_row_s[0] = 0; // S_0,0
        // For S_0,j (all gaps in X)
        for (int j = 1; j <= n_len; ++j) {
            prev_row_e[j] = (j == 1) ? (prev_row_s[j-1] + GAP_OPEN) : (prev_row_e[j-1] + GAP_EXTEND);
            prev_row_s[j] = prev_row_e[j];
        }
    } else { // Local
        // All initial scores are 0 for S, E, F for local alignment.
        std::fill(prev_row_s.begin(), prev_row_s.end(), 0);
        std::fill(prev_row_e.begin(), prev_row_e.end(), 0);
        std::fill(prev_row_f.begin(), prev_row_f.end(), 0);
    }
}

void computeAffineDPRow(int i_row, const std::string& x_str, const std::string& y_str,
                        std::vector<int>& prev_s_row, std::vector<int>& prev_e_row, std::vector<int>& prev_f_row,
                        std::vector<int>& curr_s_row, std::vector<int>& curr_e_row, std::vector<int>& curr_f_row,
                        std::vector<char>& curr_trace_s_row, ScoreFn score_fn_local, bool isGlobal) {
    int n_len = y_str.length();
    curr_s_row.assign(n_len + 1, 0); curr_e_row.assign(n_len + 1, 0); curr_f_row.assign(n_len + 1, 0);
    curr_trace_s_row.assign(n_len + 1, 'X');

    // Initialize S_i,0, E_i,0, F_i,0 (first column of current row)
    if (isGlobal) {
        // F_i,0 = S_i-1,0 + GAP_OPEN or F_i-1,0 + GAP_EXTEND (gap in Y)
        int f_i0_open = prev_s_row[0] + GAP_OPEN;
        int f_i0_extend = prev_f_row[0] + GAP_EXTEND; // prev_f_row is F for (i-1)th row
        curr_f_row[0] = std::max(f_i0_open, f_i0_extend);
        curr_s_row[0] = curr_f_row[0];
        curr_e_row[0] = INT_MIN / 2; // E_i,0 (gap in X before Y starts) is not typical for global start
        curr_trace_s_row[0] = (curr_f_row[0] == f_i0_open && curr_f_row[0] >= f_i0_extend) ? 'F' : 'f'; // From F state
    } else { // Local: S_i,0 = E_i,0 = F_i,0 = 0
        curr_s_row[0] = 0; curr_e_row[0] = 0; curr_f_row[0] = 0; curr_trace_s_row[0] = 'X';
    }

    for (int j_col = 1; j_col <= n_len; ++j_col) {
        // Match/Mismatch score M_i,j
        int s_diag_prev = prev_s_row[j_col-1];
        int e_diag_prev = prev_e_row[j_col-1]; // E from (i-1, j-1)
        int f_diag_prev = prev_f_row[j_col-1]; // F from (i-1, j-1)
        int match_pred_score = std::max({s_diag_prev, e_diag_prev, f_diag_prev});
        int m_val = match_pred_score + score_fn_local(x_str[i_row-1], y_str[j_col-1]);

        // E_i,j (gap in X, horizontal move from (i,j-1) )
        // E_i,j = max( S_i,j-1 + GAP_OPEN, E_i,j-1 + GAP_EXTEND )
        int e_open_score = curr_s_row[j_col-1] + GAP_OPEN;
        int e_extend_score = curr_e_row[j_col-1] + GAP_EXTEND;
        curr_e_row[j_col] = std::max(e_open_score, e_extend_score);

        // F_i,j (gap in Y, vertical move from (i-1,j) )
        // F_i,j = max( S_i-1,j + GAP_OPEN, F_i-1,j + GAP_EXTEND )
        int f_open_score = prev_s_row[j_col] + GAP_OPEN;
        int f_extend_score = prev_f_row[j_col] + GAP_EXTEND;
        curr_f_row[j_col] = std::max(f_open_score, f_extend_score);

        if (!isGlobal) { // For local, floor E and F at 0
            curr_e_row[j_col] = std::max(0, curr_e_row[j_col]);
            curr_f_row[j_col] = std::max(0, curr_f_row[j_col]);
            m_val = std::max(0, m_val); // also floor M before taking max for S
        }

        // S_i,j = max(M_i,j, E_i,j, F_i,j) (and 0 for local)
        curr_s_row[j_col] = std::max({m_val, curr_e_row[j_col], curr_f_row[j_col]});
        if (!isGlobal) curr_s_row[j_col] = std::max(0, curr_s_row[j_col]);


        // Determine pointer for S matrix
        if (curr_s_row[j_col] > 0 || (isGlobal && (curr_s_row[j_col] == m_val || curr_s_row[j_col] == curr_e_row[j_col] || curr_s_row[j_col] == curr_f_row[j_col]))) {
             if (curr_s_row[j_col] == m_val) curr_trace_s_row[j_col] = 'M';
             // Tie-breaking: M > E > F
             else if (curr_s_row[j_col] == curr_e_row[j_col]) {
                 curr_trace_s_row[j_col] = (curr_e_row[j_col] == e_open_score && curr_e_row[j_col] >= e_extend_score) ? 'E' : 'e';
             } else { // Must be F
                 curr_trace_s_row[j_col] = (curr_f_row[j_col] == f_open_score && curr_f_row[j_col] >= f_extend_score) ? 'F' : 'f';
             }
        } else { // Score is 0 for local, or some other condition for global if S becomes very negative
            curr_trace_s_row[j_col] = 'X';
        }
    }
}

// compute_local_affine_cell was integrated into perform_sw_in_window's logic
struct Loc { int score; int i; int j; };

// -------- Seed Management (Full Code) --------
struct Seed {
    int query_pos; int target_pos; int len;
    bool operator<(const Seed& other) const { if (query_pos != other.query_pos) return query_pos < other.query_pos; if (target_pos != other.target_pos) return target_pos < other.target_pos; return len < other.len; }
    int query_end() const { return query_pos + len - 1; }
    int target_end() const { return target_pos + len - 1; }
};
struct ChainedSeed {
    std::vector<Seed> seeds; double chain_score;
    int query_chain_start() const { return seeds.empty() ? -1 : seeds.front().query_pos; }
    int query_chain_end()   const { return seeds.empty() ? -1 : seeds.back().query_end(); }
    int target_chain_start() const { return seeds.empty() ? -1 : seeds.front().target_pos; }
    int target_chain_end()   const { return seeds.empty() ? -1 : seeds.back().target_end(); }
};

std::vector<Seed> generate_raw_seeds(const std::string& query_seq, const FMIndex& target_fm_index, int kmer_len, int mpi_rank_val = 0, int mpi_num_procs_val = 1) {
    std::vector<Seed> current_seeds;
    if (kmer_len <= 0) {
        if (mpi_rank_val == 0 && verbose) std::cout << "DEBUG [generate_raw_seeds]: kmer_len <= 0 (" << kmer_len << "), returning no seeds." << std::endl;
        return current_seeds;
    }
    size_t q_len = query_seq.length();
    if (static_cast<size_t>(kmer_len) > q_len) {
        if (mpi_rank_val == 0 && verbose) std::cout << "DEBUG [generate_raw_seeds]: kmer_len (" << kmer_len << ") > query_len (" << q_len << "), returning no seeds." << std::endl;
        return current_seeds;
    }

    // (MPI distribution logic for start/end kmer index - current setup uses mpi_rank_val=0, mpi_num_procs_val=1 for rank 0 generation)
    size_t num_kmers_total = q_len - kmer_len + 1;
    size_t kmers_per_rank_chunk = num_kmers_total / mpi_num_procs_val; // Should be num_kmers_total if mpi_num_procs_val is 1
    size_t remainder_kmers_val = num_kmers_total % mpi_num_procs_val;
    size_t my_start_kmer_idx_val = mpi_rank_val * kmers_per_rank_chunk + std::min(static_cast<size_t>(mpi_rank_val), remainder_kmers_val);
    size_t my_num_kmers_to_process_val = kmers_per_rank_chunk + (static_cast<size_t>(mpi_rank_val) < remainder_kmers_val ? 1 : 0);
    size_t my_end_kmer_idx_val = my_start_kmer_idx_val + my_num_kmers_to_process_val;

    if (mpi_rank_val == 0 && verbose) {
        std::cout << "DEBUG [generate_raw_seeds]: Query (first 30): '" << query_seq.substr(0, 30) << "'" << std::endl;
        std::cout << "DEBUG [generate_raw_seeds]: Target Text from FMIndex (first 30): '" << target_fm_index.text_with_sentinel.substr(0, 30) << "'" << std::endl;
        std::cout << "DEBUG [generate_raw_seeds]: kmer_len: " << kmer_len << ". Will iterate " << my_num_kmers_to_process_val << " kmers." << std::endl;
    }

    int found_count = 0;
    for (size_t i_kmer = my_start_kmer_idx_val; i_kmer < my_end_kmer_idx_val; ++i_kmer) {
        std::string kmer_str = query_seq.substr(i_kmer, kmer_len);
        // For debugging, print only the first few attempted k-mers to avoid excessive output
        if (mpi_rank_val == 0 && verbose && i_kmer < my_start_kmer_idx_val + 5) {
            std::cout << "DEBUG [generate_raw_seeds]: Searching for kmer #" << i_kmer << ": '" << kmer_str << "'" << std::endl;
        }

        std::vector<int> t_pos_list = target_fm_index.locate(kmer_str);

        if (!t_pos_list.empty()) {
            found_count++;
            if (mpi_rank_val == 0 && verbose && found_count <= 5) { // Print info for first 5 found k-mers
                 std::cout << "DEBUG [generate_raw_seeds]: ---> FOUND kmer '" << kmer_str << "' at "
                           << t_pos_list.size() << " target positions. First target pos: " << t_pos_list[0] << std::endl;
            }
            for (int t_p_val : t_pos_list) {
                current_seeds.push_back({(int)i_kmer, t_p_val, kmer_len});
            }
        } else if (mpi_rank_val == 0 && verbose && i_kmer < my_start_kmer_idx_val + 20 && (i_kmer % 5 ==0) ) { // Print not found for some early kmers
             std::cout << "DEBUG [generate_raw_seeds]: Did NOT find kmer: '" << kmer_str << "'" << std::endl;
        }
    }
    if (mpi_rank_val == 0 && verbose) {
        std::cout << "DEBUG [generate_raw_seeds]: Total k-mers found and added as seeds: " << current_seeds.size() << std::endl;
    }
    return current_seeds;
}

ChainedSeed find_best_seed_chain(std::vector<Seed>& seeds_vec, int min_diag_gap_val = 0, int max_diag_gap_val = 50000, int max_offset_dev_val = 50) {
    if (seeds_vec.empty()) { return {}; }
    std::sort(seeds_vec.begin(), seeds_vec.end());
    int n_s_val = seeds_vec.size();
    std::vector<double> dp_scores(n_s_val); std::vector<int> prev_indices(n_s_val, -1); double max_chain_score = 0; int best_chain_end_idx = -1;
    for (int i_s = 0; i_s < n_s_val; ++i_s) { dp_scores[i_s] = seeds_vec[i_s].len;
        for (int j_s = 0; j_s < i_s; ++j_s) {
            if (seeds_vec[j_s].query_end() + min_diag_gap_val < seeds_vec[i_s].query_pos && seeds_vec[j_s].target_end() + min_diag_gap_val < seeds_vec[i_s].target_pos) {
                if (std::abs(((long long)seeds_vec[i_s].query_pos - seeds_vec[i_s].target_pos) - ((long long)seeds_vec[j_s].query_pos - seeds_vec[j_s].target_pos)) > max_offset_dev_val) continue;
                if ((seeds_vec[i_s].query_pos - seeds_vec[j_s].query_end() > max_diag_gap_val) || (seeds_vec[i_s].target_pos - seeds_vec[j_s].target_end() > max_diag_gap_val)) continue;
                if (dp_scores[j_s] + seeds_vec[i_s].len > dp_scores[i_s]) { dp_scores[i_s] = dp_scores[j_s] + seeds_vec[i_s].len; prev_indices[i_s] = j_s; }
            }
        } if (dp_scores[i_s] > max_chain_score) { max_chain_score = dp_scores[i_s]; best_chain_end_idx = i_s; }
    } ChainedSeed result_chain; result_chain.chain_score = 0;
    if (best_chain_end_idx != -1) { result_chain.chain_score = max_chain_score; int current_s_idx = best_chain_end_idx; while (current_s_idx != -1) { result_chain.seeds.push_back(seeds_vec[current_s_idx]); current_s_idx = prev_indices[current_s_idx]; } std::reverse(result_chain.seeds.begin(), result_chain.seeds.end()); }
    return result_chain;
}

// -------- Segment/Window Alignment Helpers & Structs (Full Code) --------
struct AlignmentResult {
    std::string aligned_seq1; std::string aligned_seq2; int score = 0; long long time_ms = 0;
    int query_start_orig = -1, query_end_orig = -1; int target_start_orig = -1, target_end_orig = -1;
};
struct LcsSegmentResult {
    std::string lcs_string; int lcs_length = 0; std::string gapped_seq1; std::string gapped_seq2;
};

AlignmentResult perform_sw_in_window(const std::string& sub1, const std::string& sub2, /*ScoreMode sm,*/ ScoreFn sfn, double go, double ge, int q_off, int t_off) {    AlignmentResult res; res.score = 0; int m_len = sub1.length(); int n_len = sub2.length(); if (m_len == 0 || n_len == 0) return res;
    std::vector<std::vector<int>> s_mat(m_len + 1, std::vector<int>(n_len + 1, 0)); std::vector<std::vector<int>> e_mat(m_len + 1, std::vector<int>(n_len + 1, 0));
    std::vector<std::vector<int>> f_mat(m_len + 1, std::vector<int>(n_len + 1, 0));
    int max_s_val = 0; int max_i_s_val = 0, max_j_s_val = 0; // Use _val to avoid conflict with global rank
    for (int i = 1; i <= m_len; ++i) { for (int j = 1; j <= n_len; ++j) {
        int m_pred = std::max({s_mat[i-1][j-1], e_mat[i-1][j-1], f_mat[i-1][j-1]}); int m_v = m_pred + sfn(sub1[i-1], sub2[j-1]);
        int e_o = s_mat[i][j-1] + go; int e_e = e_mat[i][j-1] + ge; e_mat[i][j] = std::max({0, e_o, e_e});
        int f_o = s_mat[i-1][j] + go; int f_e = f_mat[i-1][j] + ge; f_mat[i][j] = std::max({0, f_o, f_e});
        s_mat[i][j] = std::max({0, m_v, e_mat[i][j], f_mat[i][j]});
        if (s_mat[i][j] > max_s_val) { max_s_val = s_mat[i][j]; max_i_s_val = i; max_j_s_val = j; }
    }} res.score = max_s_val;
    if (max_s_val > 0) { std::string r_a1, r_a2; int ci = max_i_s_val, cj = max_j_s_val; int current_state = 0; // 0=S, 1=E, 2=F
        // Determine starting state based on which matrix contributed to s_mat[ci][cj]
        int m_check = (std::max({(ci>0&&cj>0?s_mat[ci-1][cj-1]:INT_MIN/2), (ci>0&&cj>0?e_mat[ci-1][cj-1]:INT_MIN/2), (ci>0&&cj>0?f_mat[ci-1][cj-1]:INT_MIN/2)}) + (ci>0&&cj>0?sfn(sub1[ci-1],sub2[cj-1]):0) );
        if (ci>0 && cj>0 && s_mat[ci][cj] == m_check && s_mat[ci][cj] >= e_mat[ci][cj] && s_mat[ci][cj] >= f_mat[ci][cj]) current_state = 0; // Prefer M if equal
        else if (s_mat[ci][cj] == e_mat[ci][cj] && s_mat[ci][cj] >= f_mat[ci][cj]) current_state = 1;
        else if (s_mat[ci][cj] == f_mat[ci][cj]) current_state = 2;
        else if (s_mat[ci][cj] == 0) { /* No path */ } // Should not happen if max_s_val > 0

        while(s_mat[ci][cj] > 0 && (ci > 0 || cj > 0)) {
            if (current_state == 0) { // From S (Match/Mismatch)
                if (ci <=0 || cj <=0) break;
                r_a1 += sub1[ci-1]; r_a2 += sub2[cj-1];
                int prev_s = (ci>1&&cj>1?s_mat[ci-1][cj-1]:INT_MIN/2); int prev_e = (ci>1&&cj>1?e_mat[ci-1][cj-1]:INT_MIN/2); int prev_f = (ci>1&&cj>1?f_mat[ci-1][cj-1]:INT_MIN/2);
                ci--; cj--;
                if(ci<0 || cj<0) break; // Check bounds again after decrement
                // Update current_state based on which of S, E, F at (ci, cj) was max for the M_pred
                if(prev_s >= prev_e && prev_s >= prev_f) current_state = 0;
                else if (prev_e >= prev_f) current_state = 1;
                else current_state = 2;

            } else if (current_state == 1) { // From E (Gap in sub1, horizontal)
                if (ci <= 0) { break; }
                r_a1 += sub1[ci-1]; r_a2 += '-';
                if (e_mat[ci][cj] == (s_mat[ci][cj-1] + go) && e_mat[ci][cj] >= (e_mat[ci][cj-1] + ge) ) current_state = 0; // Opened from S
                // else: extended from E (stay in state 1, handled by cj--)
                cj--;
            } else { // current_state == 2, From F (Gap in sub2, vertical)
                if (ci <= 0) { break; }
                r_a1 += sub1[ci-1];
                r_a2 += '-';
                if (f_mat[ci][cj] == (s_mat[ci-1][cj] + go) && f_mat[ci][cj] >= (f_mat[ci-1][cj] + ge)) current_state = 0; // Opened from S
                // else: extended from F
                ci--;
            }
        } std::reverse(r_a1.begin(), r_a1.end()); std::reverse(r_a2.begin(), r_a2.end()); res.aligned_seq1 = r_a1; res.aligned_seq2 = r_a2;
        res.query_end_orig = q_off + max_i_s_val - 1; res.target_end_orig = t_off + max_j_s_val - 1;
        int q_chars_count = 0; for(char c_val : res.aligned_seq1) if(c_val != '-') q_chars_count++;
        int t_chars_count = 0; for(char c_val : res.aligned_seq2) if(c_val != '-') t_chars_count++;
        res.query_start_orig = q_off + (max_i_s_val - q_chars_count); res.target_start_orig = t_off + (max_j_s_val - t_chars_count);
    } return res;
}

AlignmentResult align_segment_globally(const std::string& seg1, const std::string& seg2, /*ScoreMode sm,*/ ScoreFn sfn, double go, double ge) {
    AlignmentResult res; int m_len = seg1.length(); int n_len = seg2.length();
    if (m_len==0&&n_len==0){res.score=0; return res;} if(m_len==0){res.aligned_seq1=std::string(n_len,'-');res.aligned_seq2=seg2;res.score=go+(n_len>1?(n_len-1)*ge:0);if(n_len==0)res.score=0;return res;} if(n_len==0){res.aligned_seq1=seg1;res.aligned_seq2=std::string(m_len,'-');res.score=go+(m_len>1?(m_len-1)*ge:0);if(m_len==0)res.score=0;return res;}
    std::vector<int> ps_r(n_len+1), pe_r(n_len+1,0), pf_r(n_len+1,0); std::vector<std::vector<char>> tr(m_len+1,std::vector<char>(n_len+1));
    ps_r[0]=0; pe_r[0]=INT_MIN/2; pf_r[0]=INT_MIN/2; tr[0][0]='S';
    for(int j=1;j<=n_len;++j){pe_r[j]=(j==1?(ps_r[j-1]+go):(pe_r[j-1]+ge)); ps_r[j]=pe_r[j]; pf_r[j]=INT_MIN/2; tr[0][j]=(j==1&&pe_r[j]==(ps_r[j-1]+go))?'E':'e';} // Tie break to 'e'
    std::vector<int> cs_r(n_len+1), ce_r(n_len+1), cf_r(n_len+1);
    for(int i=1;i<=m_len;++i){
        int f0_o=ps_r[0]+go; int f0_e=pf_r[0]+ge; cf_r[0]=std::max(f0_o,f0_e); cs_r[0]=cf_r[0]; ce_r[0]=INT_MIN/2; tr[i][0]=(cf_r[0]==f0_o && cf_r[0] >= f0_e)?'F':'f'; // Tie break to 'f'
        for(int j=1;j<=n_len;++j){
            int s_d_val=ps_r[j-1];int e_d_val=pe_r[j-1];int f_d_val=pf_r[j-1]; int m_p_val=std::max({s_d_val,e_d_val,f_d_val}); int m_v_val=m_p_val+sfn(seg1[i-1],seg2[j-1]);
            int e_o_val=cs_r[j-1]+go; int e_e_val=ce_r[j-1]+ge; ce_r[j]=std::max(e_o_val,e_e_val);
            int f_o_val=ps_r[j]+go; int f_e_val=pf_r[j]+ge; cf_r[j]=std::max(f_o_val,f_e_val);
            cs_r[j]=std::max({m_v_val,ce_r[j],cf_r[j]});
            if(cs_r[j]==m_v_val)tr[i][j]='M'; else if(cs_r[j]==ce_r[j])tr[i][j]=(ce_r[j]==e_o_val && ce_r[j] >= e_e_val)?'E':'e'; else tr[i][j]=(cf_r[j]==f_o_val && cf_r[j] >= f_e_val)?'F':'f';
        } ps_r.swap(cs_r); pe_r.swap(ce_r); pf_r.swap(cf_r); // Efficient row swap
    } res.score=ps_r[n_len]; std::string r_a1,r_a2; int ci_val=m_len,cj_val=n_len; // Renamed ci, cj
    while(ci_val>0||cj_val>0){ char move_char=tr[ci_val][cj_val];
        if(move_char=='M'){r_a1+=seg1[ci_val-1];r_a2+=seg2[cj_val-1];ci_val--;cj_val--;}
        else if(move_char=='E'||move_char=='e'){r_a1+='-';r_a2+=seg2[cj_val-1];cj_val--;}
        else if(move_char=='F'||move_char=='f'){r_a1+=seg1[ci_val-1];r_a2+='-';ci_val--;}
        else if(ci_val==0&&cj_val>0){r_a1+='-';r_a2+=seg2[cj_val-1];cj_val--;} else if(cj_val==0&&ci_val>0){r_a1+=seg1[ci_val-1];r_a2+='-';ci_val--;} else break;
    } std::reverse(r_a1.begin(),r_a1.end());std::reverse(r_a2.begin(),r_a2.end()); res.aligned_seq1=r_a1;res.aligned_seq2=r_a2; return res;
}

LcsSegmentResult compute_lcs_for_segment(const std::string& seg1, const std::string& seg2) {
    LcsSegmentResult res; int m_len=seg1.length(); int n_len=seg2.length();
    if(m_len==0||n_len==0){res.lcs_length=0;res.lcs_string=""; if(m_len==0&&n_len>0){res.gapped_seq1=std::string(n_len,'-');res.gapped_seq2=seg2;}else if(n_len==0&&m_len>0){res.gapped_seq1=seg1;res.gapped_seq2=std::string(m_len,'-');} return res;}
    std::vector<std::vector<int>>L_mat(m_len+1,std::vector<int>(n_len+1,0));std::vector<std::vector<char>>B_mat(m_len+1,std::vector<char>(n_len+1,' ')); // Renamed L, B
    for(int i=1;i<=m_len;++i){for(int j=1;j<=n_len;++j){if(seg1[i-1]==seg2[j-1]){L_mat[i][j]=L_mat[i-1][j-1]+1;B_mat[i][j]='D';}else if(L_mat[i-1][j]>=L_mat[i][j-1]){L_mat[i][j]=L_mat[i-1][j];B_mat[i][j]='U';}else{L_mat[i][j]=L_mat[i][j-1];B_mat[i][j]='L';}}}
    res.lcs_length=L_mat[m_len][n_len];std::string l_rev,g1_rev,g2_rev;int ci_lcs=m_len,cj_lcs=n_len; // Renamed ci, cj
    while(ci_lcs>0||cj_lcs>0){if(ci_lcs>0&&cj_lcs>0&&B_mat[ci_lcs][cj_lcs]=='D'){l_rev+=seg1[ci_lcs-1];g1_rev+=seg1[ci_lcs-1];g2_rev+=seg2[cj_lcs-1];ci_lcs--;cj_lcs--;}else if(ci_lcs>0&&(cj_lcs==0||B_mat[ci_lcs][cj_lcs]=='U')){g1_rev+=seg1[ci_lcs-1];g2_rev+='-';ci_lcs--;}else if(cj_lcs>0&&(ci_lcs==0||B_mat[ci_lcs][cj_lcs]=='L')){g1_rev+='-';g2_rev+=seg2[cj_lcs-1];cj_lcs--;}else break;}
    std::reverse(l_rev.begin(),l_rev.end());res.lcs_string=l_rev;std::reverse(g1_rev.begin(),g1_rev.end());res.gapped_seq1=g1_rev;std::reverse(g2_rev.begin(),g2_rev.end());res.gapped_seq2=g2_rev; return res;
}

// -------- Refactored Main Alignment Functions (with conceptual MPI) --------
void globalalign(const std::string &x_orig, const std::string &y_orig, const std::string &header1, const std::string &header2, const std::string &outdir, ScoreMode mode_val, ScoreFn score_fn_val, const FMIndex* target_fm_idx) {
    int m_o = x_orig.length(); int n_o = y_orig.length(); int rnk_val = ::rank; int n_procs_val = 0; MPI_Comm_size(MPI_COMM_WORLD, &n_procs_val); // Use global rank
    auto t_start_ga = std::chrono::high_resolution_clock::now(); AlignmentResult final_aln_ga; final_aln_ga.score = 0; ChainedSeed best_chain_ga;
    bool use_anchoring = false;

    if(target_fm_idx){
        if(rnk_val==0 && verbose)std::cout<<"Global alignment: Using FM-Index for anchoring."<<std::endl;
        if(rnk_val==0){
            int k_anchor_ga=std::min({8, m_o/12, n_o/12}); if(k_anchor_ga<8 && std::min(m_o,n_o)>=8)k_anchor_ga=std::min(8,std::min(m_o,n_o)); if(k_anchor_ga<=0 && std::min(m_o,n_o)>0) k_anchor_ga=std::min(std::min(m_o,n_o), 8); // ensure k is at least 1 if possible
            if(k_anchor_ga > 0) {
                std::vector<Seed> raw_sds_ga=generate_raw_seeds(x_orig,*target_fm_idx,k_anchor_ga,0,1);
                if(!raw_sds_ga.empty()){best_chain_ga=find_best_seed_chain(raw_sds_ga,1); if(verbose && !best_chain_ga.seeds.empty())std::cout<<"Rank 0: Global anchor chain found: " << best_chain_ga.seeds.size() << " anchors, score " << best_chain_ga.chain_score << std::endl;}
                else if(verbose) std::cout << "Rank 0: No raw seeds for global anchoring." << std::endl;
            } else if(verbose) std::cout << "Rank 0: k-mer for global anchors too small." << std::endl;
        }
        int num_anchors_ga=(rnk_val==0)?best_chain_ga.seeds.size():0; MPI_Bcast(&num_anchors_ga,1,MPI_INT,0,MPI_COMM_WORLD);
        if(num_anchors_ga > 0){
            use_anchoring = true;
            if(rnk_val!=0) { best_chain_ga.seeds.resize(num_anchors_ga); }
            MPI_Bcast(best_chain_ga.seeds.data(),num_anchors_ga*sizeof(Seed),MPI_BYTE,0,MPI_COMM_WORLD);
            // MPI_Bcast(&best_chain_ga.chain_score, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD); // Bcast score if needed

            std::vector<std::pair<std::string,std::string>>segs_ga; int cur_x_ga=0,cur_y_ga=0;
            for(const auto&anc_ga:best_chain_ga.seeds){segs_ga.push_back({x_orig.substr(cur_x_ga,anc_ga.query_pos-cur_x_ga),y_orig.substr(cur_y_ga,anc_ga.target_pos-cur_y_ga)});cur_x_ga=anc_ga.query_pos+anc_ga.len;cur_y_ga=anc_ga.target_pos+anc_ga.len;}
            segs_ga.push_back({x_orig.substr(cur_x_ga),y_orig.substr(cur_y_ga)});
            std::vector<AlignmentResult>seg_res_ga(segs_ga.size());

            // TODO: Implement robust MPI for distributing segment alignments and gathering AlignmentResult (strings!)
            // Simplified: each rank does its share, rank 0 reconstructs. This requires seg_res_ga to be gathered.
            for(size_t i_ga=0; i_ga < segs_ga.size(); ++i_ga){
                if (i_ga % static_cast<size_t>(n_procs_val) == static_cast<size_t>(rnk_val)) { // Static distribution
                     if (verbose) std::cout << "Rank " << rnk_val << " aligning global segment " << i_ga << std::endl;
                    seg_res_ga[i_ga]=align_segment_globally(segs_ga[i_ga].first,segs_ga[i_ga].second, /* mode_val, */ score_fn_val,GAP_OPEN,GAP_EXTEND);
                }
            }
            // Gather results to rank 0
            if(rnk_val==0){
                for(size_t i_seg_ga = 0; i_seg_ga < segs_ga.size(); ++i_seg_ga){
                    int owner_rank = static_cast<int>(i_seg_ga % static_cast<size_t>(n_procs_val)); // Corrected calculation of owner_rank
                    if (owner_rank != 0) { // Receive from worker ranks
                        // MPI_Recv score, lengths of aligned_seq1, aligned_seq2, then the strings.
                        // This is a placeholder for complex MPI string communication.
                        MPI_Status recv_status_ga;
                        MPI_Recv(&seg_res_ga[i_seg_ga].score, 1, MPI_INT, owner_rank, i_seg_ga * 10 + 0, MPI_COMM_WORLD, &recv_status_ga);
                        int len1_ga, len2_ga;
                        MPI_Recv(&len1_ga, 1, MPI_INT, owner_rank, i_seg_ga * 10 + 1, MPI_COMM_WORLD, &recv_status_ga);
                        seg_res_ga[i_seg_ga].aligned_seq1.resize(len1_ga);
                        if (len1_ga > 0) MPI_Recv(&seg_res_ga[i_seg_ga].aligned_seq1[0], len1_ga, MPI_CHAR, owner_rank, i_seg_ga * 10 + 2, MPI_COMM_WORLD, &recv_status_ga);
                        MPI_Recv(&len2_ga, 1, MPI_INT, owner_rank, i_seg_ga * 10 + 3, MPI_COMM_WORLD, &recv_status_ga);
                        seg_res_ga[i_seg_ga].aligned_seq2.resize(len2_ga);
                        if (len2_ga > 0) MPI_Recv(&seg_res_ga[i_seg_ga].aligned_seq2[0], len2_ga, MPI_CHAR, owner_rank, i_seg_ga * 10 + 4, MPI_COMM_WORLD, &recv_status_ga);
                    }
                }
                // Rank 0 reconstructs
                for(size_t i_ga_rec=0;i_ga_rec<best_chain_ga.seeds.size();++i_ga_rec){
                    final_aln_ga.aligned_seq1+=seg_res_ga[i_ga_rec].aligned_seq1; final_aln_ga.aligned_seq2+=seg_res_ga[i_ga_rec].aligned_seq2; final_aln_ga.score+=seg_res_ga[i_ga_rec].score;
                    const auto&anc_ga_rec=best_chain_ga.seeds[i_ga_rec];std::string as_ga=x_orig.substr(anc_ga_rec.query_pos,anc_ga_rec.len);
                    final_aln_ga.aligned_seq1+=as_ga;final_aln_ga.aligned_seq2+=as_ga;for(char c_ga:as_ga)final_aln_ga.score+=score_fn_val(c_ga,c_ga);
                }
                final_aln_ga.aligned_seq1+=seg_res_ga.back().aligned_seq1;final_aln_ga.aligned_seq2+=seg_res_ga.back().aligned_seq2;final_aln_ga.score+=seg_res_ga.back().score;
            } else { // Worker ranks send their computed segments
                 for(size_t i_ga=0; i_ga < segs_ga.size(); ++i_ga){
                     if (i_ga % static_cast<size_t>(n_procs_val) == static_cast<size_t>(rnk_val)) {
                        MPI_Send(&seg_res_ga[i_ga].score, 1, MPI_INT, 0, i_ga * 10 + 0, MPI_COMM_WORLD);
                        int len1_ga_send = seg_res_ga[i_ga].aligned_seq1.length(); int len2_ga_send = seg_res_ga[i_ga].aligned_seq2.length();
                        MPI_Send(&len1_ga_send, 1, MPI_INT, 0, i_ga * 10 + 1, MPI_COMM_WORLD);
                        if (len1_ga_send > 0) MPI_Send(seg_res_ga[i_ga].aligned_seq1.data(), len1_ga_send, MPI_CHAR, 0, i_ga * 10 + 2, MPI_COMM_WORLD);
                        MPI_Send(&len2_ga_send, 1, MPI_INT, 0, i_ga * 10 + 3, MPI_COMM_WORLD);
                        if (len2_ga_send > 0) MPI_Send(seg_res_ga[i_ga].aligned_seq2.data(), len2_ga_send, MPI_CHAR, 0, i_ga * 10 + 4, MPI_COMM_WORLD);
                     }
                 }
            }
        } else { use_anchoring = false; } // No anchors, will fallback
    }

    if(!use_anchoring){ // Fallback to original full DP
        if(rnk_val==0 && verbose)std::cout<<"Global alignment: No anchors or FM-Index not used. Fallback to full DP (conceptual)."<<std::endl;
        // TODO: Implement the original full Needleman-Wunsch MPI logic here.
        // This would fill final_aln_ga.
        // For now, just a placeholder indicating it would run.
        if (rnk_val == 0) {
            final_aln_ga = align_segment_globally(x_orig, y_orig, score_fn_val, GAP_OPEN, GAP_EXTEND); // Example of CPU full version
        }
        // MPI_Bcast score, and string lengths, then strings if rank 0 did it.
    }

    auto t_end_ga=std::chrono::high_resolution_clock::now();final_aln_ga.time_ms=std::chrono::duration_cast<std::chrono::milliseconds>(t_end_ga-t_start_ga).count();
    if(rnk_val==0){ if(verbose){ std::cout<<"\n--- Global Alignment Final Score: "<<final_aln_ga.score<<"\n"; printColoredAlignment(final_aln_ga.aligned_seq1,final_aln_ga.aligned_seq2); std::cout<<"Time: "<<final_aln_ga.time_ms<<"ms\n";}
        // Save files (plain alignment, JSON stats)
        std::string acc1_ga = getAccession(header1, mode_val); std::string acc2_ga = getAccession(header2, mode_val);
        std::ofstream outf_ga(outdir + "/global_alignment.fasta"); if (outf_ga) savePlainAlignment(acc1_ga, acc2_ga, final_aln_ga.aligned_seq1, final_aln_ga.aligned_seq2, outf_ga);
        // ... JSON stats save ...
    } MPI_Barrier(MPI_COMM_WORLD);
}

void localalign(const std::string &x, const std::string &y, const std::string &h1, const std::string &h2, const std::string &odir, ScoreMode mval, ScoreFn sfn_val, const FMIndex* tfm_idx) {
    int m_len=x.length();int n_len=y.length();int r_la=::rank;int np_la=0;MPI_Comm_size(MPI_COMM_WORLD, &np_la); auto ts_la=std::chrono::high_resolution_clock::now();
    AlignmentResult g_best_aln_la; g_best_aln_la.score = 0; bool use_fmindex_la = false;

    if(tfm_idx){
        if(r_la==0&&verbose)std::cout<<"Local alignment: Using FM-Index seed-and-extend."<<std::endl;
        std::vector<Seed>all_sds_la; int k_s_la=std::min({11,m_len/20,n_len/20}); if(k_s_la<7&&std::min(m_len,n_len)>=7)k_s_la=std::min(7,std::min(m_len,n_len)); if(k_s_la<=0&&std::min(m_len,n_len)>0) k_s_la=std::min(std::min(m_len,n_len),7);
        if(k_s_la > 0) {
            if(r_la==0) { all_sds_la=generate_raw_seeds(x,*tfm_idx,k_s_la,0,1); if(verbose) std::cout << "Rank 0: Local generated " << all_sds_la.size() << " raw seeds, k=" << k_s_la << std::endl; }
            int tot_sds_la=(r_la==0)?all_sds_la.size():0; MPI_Bcast(&tot_sds_la,1,MPI_INT,0,MPI_COMM_WORLD);
            if(tot_sds_la>0){
                use_fmindex_la = true;
                if(r_la!=0) {
                    all_sds_la.resize(tot_sds_la);
                }
                MPI_Bcast(all_sds_la.data(), tot_sds_la * sizeof(Seed), MPI_BYTE, 0, MPI_COMM_WORLD);
                AlignmentResult r_best_aln_la; r_best_aln_la.score=0;
                int s_start_la=(tot_sds_la/np_la*r_la)+std::min(r_la,tot_sds_la%np_la); int s_end_la=(tot_sds_la/np_la*(r_la+1))+std::min(r_la+1,tot_sds_la%np_la);
                long num_seeds_this_rank = static_cast<long>(s_end_la - s_start_la); // Calculate the actual number of seeds
                if (verbose && num_seeds_this_rank > 0) {
                    std::cout << "Rank " << r_la << " processing " << num_seeds_this_rank << " seeds for local." << std::endl;
                }
                for(int si_la=s_start_la;si_la<s_end_la;++si_la){const auto&cs_la=all_sds_la[si_la]; int wrq_la=std::max(100,cs_la.len*3);int wrt_la=std::max(100,cs_la.len*3); // Increased window
                    int qws_la=std::max(0,cs_la.query_pos-wrq_la);int qwe_la=std::min(m_len,cs_la.query_pos+cs_la.len+wrq_la); int tws_la=std::max(0,cs_la.target_pos-wrt_la);int twe_la=std::min(n_len,cs_la.target_pos+cs_la.len+wrt_la);
                    std::string sx_la=x.substr(qws_la,qwe_la-qws_la); std::string sy_la=y.substr(tws_la,twe_la-tws_la); if(sx_la.empty()||sy_la.empty())continue;
                    AlignmentResult cwa_la=perform_sw_in_window(sx_la,sy_la, sfn_val, GAP_OPEN,GAP_EXTEND,qws_la,tws_la); if(cwa_la.score>r_best_aln_la.score)r_best_aln_la=cwa_la;
                }
                struct SR_la {int s; int rnk_val;}; SR_la my_sr_la={r_best_aln_la.score,r_la}; SR_la g_max_sr_la={0,0}; MPI_Allreduce(&my_sr_la,&g_max_sr_la,1,MPI_2INT,MPI_MAXLOC,MPI_COMM_WORLD);
                g_best_aln_la.score=g_max_sr_la.s; int best_r_la=g_max_sr_la.rnk_val;
                // Rank best_r_la sends its full AlignmentResult to rank 0
                if(g_best_aln_la.score>0){
//                    char* aln1_buf = nullptr; char* aln2_buf = nullptr;
                    int aln1_len=0, aln2_len=0;
                    if(r_la==best_r_la){g_best_aln_la=r_best_aln_la; aln1_len=g_best_aln_la.aligned_seq1.length(); aln2_len=g_best_aln_la.aligned_seq2.length();}
                    MPI_Bcast(&aln1_len, 1, MPI_INT, best_r_la, MPI_COMM_WORLD); MPI_Bcast(&aln2_len, 1, MPI_INT, best_r_la, MPI_COMM_WORLD);
                    MPI_Bcast(&g_best_aln_la.query_start_orig, 1, MPI_INT, best_r_la, MPI_COMM_WORLD); MPI_Bcast(&g_best_aln_la.query_end_orig, 1, MPI_INT, best_r_la, MPI_COMM_WORLD);
                    MPI_Bcast(&g_best_aln_la.target_start_orig, 1, MPI_INT, best_r_la, MPI_COMM_WORLD); MPI_Bcast(&g_best_aln_la.target_end_orig, 1, MPI_INT, best_r_la, MPI_COMM_WORLD);
                    if(r_la!=best_r_la && g_best_aln_la.score > 0){ g_best_aln_la.aligned_seq1.resize(aln1_len); g_best_aln_la.aligned_seq2.resize(aln2_len);}
                    if(g_best_aln_la.score > 0){
                        if(aln1_len>0) MPI_Bcast(&g_best_aln_la.aligned_seq1[0], aln1_len, MPI_CHAR, best_r_la, MPI_COMM_WORLD);
                        if(aln2_len>0) MPI_Bcast(&g_best_aln_la.aligned_seq2[0], aln2_len, MPI_CHAR, best_r_la, MPI_COMM_WORLD);
                    }
                }
            } else {if(r_la==0&&verbose) std::cout<<"Local: No seeds generated."<<std::endl; use_fmindex_la = false;}
        } else {if(r_la==0&&verbose) std::cout<<"Local: k-mer for seeding too small."<<std::endl; use_fmindex_la = false;}
    }

    if(!use_fmindex_la){ if(r_la==0&&verbose)std::cout<<"Local alignment: FM-Index not used. Fallback to full DP (conceptual)."<<std::endl;
        // TODO: Implement the original full Smith-Waterman MPI logic here.
        // This would populate g_best_aln_la on rank 0.
    }
    auto te_la=std::chrono::high_resolution_clock::now();g_best_aln_la.time_ms=std::chrono::duration_cast<std::chrono::milliseconds>(te_la-ts_la).count();
    if(r_la==0){if(verbose){std::cout<<"\n--- Local Alignment Best Score: "<<g_best_aln_la.score<<"\n"; if(g_best_aln_la.score>0)printColoredAlignment(g_best_aln_la.aligned_seq1,g_best_aln_la.aligned_seq2); std::cout<<"Time: "<<g_best_aln_la.time_ms<<"ms\n";}
        // Save files for local alignment
        std::string acc1_la = getAccession(h1, mval); std::string acc2_la = getAccession(h2, mval);
        if (g_best_aln_la.score > 0) { std::ofstream outf_la(odir + "/local_alignment.fasta"); if (outf_la) savePlainAlignment(acc1_la + "_local", acc2_la + "_local", g_best_aln_la.aligned_seq1, g_best_aln_la.aligned_seq2, outf_la); }
        // ... JSON stats save for local...
    } MPI_Barrier(MPI_COMM_WORLD);
}

void lcs(const std::string &x_o, const std::string &y_o, const std::string &h1_lcs, const std::string &h2_lcs, const std::string &odir_lcs, ScoreMode mode_lcs, const FMIndex* tfm_idx_lcs) {
    int m_l=x_o.length();int n_l=y_o.length(); int r_lcs=::rank;int np_lcs=0;MPI_Comm_size(MPI_COMM_WORLD, &np_lcs); auto ts_lcs=std::chrono::high_resolution_clock::now();
    LcsSegmentResult final_lcs_overall_res; ChainedSeed bc_lcs_val; bool use_anchoring_lcs = false;

    if(tfm_idx_lcs){if(r_lcs==0&&verbose)std::cout<<"LCS: Using FM-Index for anchoring."<<std::endl;
        if(r_lcs==0){int k_lcs=std::min({10,m_l/15,n_l/15});if(k_lcs<5&&std::min(m_l,n_l)>=5)k_lcs=std::min(5,std::min(m_l,n_l)); if(k_lcs<=0&&std::min(m_l,n_l)>0)k_lcs=std::min(std::min(m_l,n_l),5);
            if(k_lcs>0){std::vector<Seed>rs_lcs=generate_raw_seeds(x_o,*tfm_idx_lcs,k_lcs,0,1);if(!rs_lcs.empty()){bc_lcs_val=find_best_seed_chain(rs_lcs,1); if(verbose && !bc_lcs_val.seeds.empty())std::cout<<"Rank 0: LCS anchor chain found: "<<bc_lcs_val.seeds.size()<<" anchors, score "<<bc_lcs_val.chain_score<<std::endl;}}
            else if(verbose) std::cout << "Rank 0: k-mer for LCS anchors too small." << std::endl;
        }
        int n_anc_lcs=(r_lcs==0)?bc_lcs_val.seeds.size():0; MPI_Bcast(&n_anc_lcs,1,MPI_INT,0,MPI_COMM_WORLD);
        if(n_anc_lcs>0){
            use_anchoring_lcs = true;
            if(r_lcs!=0) {
                bc_lcs_val.seeds.resize(n_anc_lcs);
            }
            MPI_Bcast(bc_lcs_val.seeds.data(), n_anc_lcs * sizeof(Seed), MPI_BYTE, 0, MPI_COMM_WORLD);
            std::vector<std::pair<std::string,std::string>>segs_lcs;int cx_lcs=0,cy_lcs=0;
            for(const auto&al_lcs:bc_lcs_val.seeds){segs_lcs.push_back({x_o.substr(cx_lcs,al_lcs.query_pos-cx_lcs),y_o.substr(cy_lcs,al_lcs.target_pos-cy_lcs)});cx_lcs=al_lcs.query_pos+al_lcs.len;cy_lcs=al_lcs.target_pos+al_lcs.len;}
            segs_lcs.push_back({x_o.substr(cx_lcs),y_o.substr(cy_lcs)});
            std::vector<LcsSegmentResult>seg_res_lcs_vec(segs_lcs.size());
            // TODO: Implement robust MPI for distributing LCS segment computations and gathering LcsSegmentResult (multiple strings!)
            for(size_t i_lcs=0; i_lcs < segs_lcs.size(); ++i_lcs){
                if(i_lcs % static_cast<size_t>(np_lcs) == static_cast<size_t>(r_lcs)){ // Static distribution
                    if(verbose) std::cout << "Rank " << r_lcs << " computing LCS for segment " << i_lcs << std::endl;
                    seg_res_lcs_vec[i_lcs]=compute_lcs_for_segment(segs_lcs[i_lcs].first,segs_lcs[i_lcs].second);
                }
            }
            // Gather to rank 0 (conceptual)
            if(r_lcs==0){
                for(size_t i_lcs_seg = 0; i_lcs_seg < segs_lcs.size(); ++i_lcs_seg) {
                    int owner = static_cast<int>(i_lcs_seg % static_cast<size_t>(np_lcs));
                    if (owner != 0) { // Receive from worker
                        // MPI_Recv LcsSegmentResult components (lengths, then strings)
                        // This is a placeholder for complex MPI string communication.
                        MPI_Status st_lcs;
                        MPI_Recv(&seg_res_lcs_vec[i_lcs_seg].lcs_length, 1, MPI_INT, owner, i_lcs_seg * 10 + 0, MPI_COMM_WORLD, &st_lcs);
                        int len_s1, len_s2, len_s3; // lcs_string, gapped1, gapped2
                        MPI_Recv(&len_s1, 1, MPI_INT, owner, i_lcs_seg * 10 + 1, MPI_COMM_WORLD, &st_lcs); seg_res_lcs_vec[i_lcs_seg].lcs_string.resize(len_s1); if(len_s1>0) MPI_Recv(&seg_res_lcs_vec[i_lcs_seg].lcs_string[0], len_s1, MPI_CHAR, owner, i_lcs_seg*10+2,MPI_COMM_WORLD, &st_lcs);
                        MPI_Recv(&len_s2, 1, MPI_INT, owner, i_lcs_seg * 10 + 3, MPI_COMM_WORLD, &st_lcs); seg_res_lcs_vec[i_lcs_seg].gapped_seq1.resize(len_s2); if(len_s2>0) MPI_Recv(&seg_res_lcs_vec[i_lcs_seg].gapped_seq1[0], len_s2, MPI_CHAR, owner, i_lcs_seg*10+4,MPI_COMM_WORLD, &st_lcs);
                        MPI_Recv(&len_s3, 1, MPI_INT, owner, i_lcs_seg * 10 + 5, MPI_COMM_WORLD, &st_lcs); seg_res_lcs_vec[i_lcs_seg].gapped_seq2.resize(len_s3); if(len_s3>0) MPI_Recv(&seg_res_lcs_vec[i_lcs_seg].gapped_seq2[0], len_s3, MPI_CHAR, owner, i_lcs_seg*10+6,MPI_COMM_WORLD, &st_lcs);
                    }
                }
                // Rank 0 reconstructs
                for(size_t i_lcs_rec=0;i_lcs_rec<bc_lcs_val.seeds.size();++i_lcs_rec){
                    final_lcs_overall_res.lcs_string+=seg_res_lcs_vec[i_lcs_rec].lcs_string; final_lcs_overall_res.lcs_length+=seg_res_lcs_vec[i_lcs_rec].lcs_length; final_lcs_overall_res.gapped_seq1+=seg_res_lcs_vec[i_lcs_rec].gapped_seq1; final_lcs_overall_res.gapped_seq2+=seg_res_lcs_vec[i_lcs_rec].gapped_seq2;
                    const auto&al_lcs_rec=bc_lcs_val.seeds[i_lcs_rec];std::string as_lcs=x_o.substr(al_lcs_rec.query_pos,al_lcs_rec.len);
                    final_lcs_overall_res.lcs_string+=as_lcs;final_lcs_overall_res.lcs_length+=al_lcs_rec.len;final_lcs_overall_res.gapped_seq1+=as_lcs;final_lcs_overall_res.gapped_seq2+=as_lcs;
                }
                final_lcs_overall_res.lcs_string+=seg_res_lcs_vec.back().lcs_string;final_lcs_overall_res.lcs_length+=seg_res_lcs_vec.back().lcs_length;final_lcs_overall_res.gapped_seq1+=seg_res_lcs_vec.back().gapped_seq1;final_lcs_overall_res.gapped_seq2+=seg_res_lcs_vec.back().gapped_seq2;
            } else { // Worker ranks send their computed LCS segments
                for(size_t i_lcs=0; i_lcs < segs_lcs.size(); ++i_lcs){
                    if(i_lcs % static_cast<size_t>(np_lcs) == static_cast<size_t>(r_lcs)){
                        MPI_Send(&seg_res_lcs_vec[i_lcs].lcs_length, 1, MPI_INT, 0, i_lcs * 10 + 0, MPI_COMM_WORLD);
                        int len_s1_send = seg_res_lcs_vec[i_lcs].lcs_string.length(); int len_s2_send = seg_res_lcs_vec[i_lcs].gapped_seq1.length(); int len_s3_send = seg_res_lcs_vec[i_lcs].gapped_seq2.length();
                        MPI_Send(&len_s1_send, 1, MPI_INT, 0, i_lcs*10+1, MPI_COMM_WORLD); if(len_s1_send>0) MPI_Send(seg_res_lcs_vec[i_lcs].lcs_string.data(), len_s1_send, MPI_CHAR, 0, i_lcs*10+2, MPI_COMM_WORLD);
                        MPI_Send(&len_s2_send, 1, MPI_INT, 0, i_lcs*10+3, MPI_COMM_WORLD); if(len_s2_send>0) MPI_Send(seg_res_lcs_vec[i_lcs].gapped_seq1.data(), len_s2_send, MPI_CHAR, 0, i_lcs*10+4, MPI_COMM_WORLD);
                        MPI_Send(&len_s3_send, 1, MPI_INT, 0, i_lcs*10+5, MPI_COMM_WORLD); if(len_s3_send>0) MPI_Send(seg_res_lcs_vec[i_lcs].gapped_seq2.data(), len_s3_send, MPI_CHAR, 0, i_lcs*10+6, MPI_COMM_WORLD);
                    }
                }
            }
        } else { use_anchoring_lcs = false; }
    }

    if(!use_anchoring_lcs){ if(r_lcs==0&&verbose)std::cout<<"LCS: No anchors or FM-Index not used. Fallback to full DP."<<std::endl;
        if(r_lcs==0) final_lcs_overall_res=compute_lcs_for_segment(x_o,y_o);
        // Bcast results if rank 0 did the fallback
    }
    auto te_lcs=std::chrono::high_resolution_clock::now();long long tms_lcs=std::chrono::duration_cast<std::chrono::milliseconds>(te_lcs-ts_lcs).count();
    if(r_lcs==0){if(verbose){std::cout<<"\n--- LCS Final Length: "<<final_lcs_overall_res.lcs_length<<"\nLCS: "<<final_lcs_overall_res.lcs_string.substr(0,100)<<"...\n"; printColoredAlignment(final_lcs_overall_res.gapped_seq1,final_lcs_overall_res.gapped_seq2);std::cout<<"Time: "<<tms_lcs<<"ms\n";}
        std::string acc1_lcs = getAccession(h1_lcs, mode_lcs); std::string acc2_lcs = getAccession(h2_lcs, mode_lcs);
        std::ofstream outf_lcs_str(odir_lcs + "/lcs.fasta"); if (outf_lcs_str) saveLCS(acc1_lcs + "_" + acc2_lcs, final_lcs_overall_res.lcs_string, outf_lcs_str);
        std::ofstream outf_lcs_aln(odir_lcs + "/lcs_alignment.fasta"); if (outf_lcs_aln) savePlainAlignment(acc1_lcs + "_LCS_aligned", acc2_lcs + "_LCS_aligned", final_lcs_overall_res.gapped_seq1, final_lcs_overall_res.gapped_seq2, outf_lcs_aln);
    } MPI_Barrier(MPI_COMM_WORLD);
}

// -------- Main Function --------
int main(int argc, char** argv) {
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &::rank); // Initialize global rank

    std::unique_ptr<FMIndex> fm_index_target_obj_ptr = nullptr;
    try {
        std::string query_file, target_file, output_dir = ".";
        std::string fm_idx_path_for_target = "";
        int align_choice = -1; ScoreMode current_mode = MODE_DNA;
        for (int i = 1; i < argc; ++i) {
            std::string current_arg = argv[i];
            if (current_arg == "--query" && i + 1 < argc) query_file = argv[++i];
            else if (current_arg == "--target" && i + 1 < argc) target_file = argv[++i];
            else if (current_arg == "--choice" && i + 1 < argc) align_choice = std::stoi(argv[++i]);
            else if (current_arg == "--mode" && i + 1 < argc) { std::string mode_str = argv[++i]; if (mode_str == "dna") current_mode = MODE_DNA; else if (mode_str == "protein") current_mode = MODE_PROTEIN; else { if (::rank == 0) std::cerr << "Unknown mode: " << mode_str << "\n"; MPI_Abort(MPI_COMM_WORLD, 1); return 1; }}
            else if (current_arg == "--outdir" && i + 1 < argc) output_dir = argv[++i];
            else if (current_arg == "--fmindex" && i + 1 < argc) fm_idx_path_for_target = argv[++i];
            else if (current_arg == "--verbose") verbose = true; else if (current_arg == "--binary") binary = true; else if (current_arg == "--txt") txt = true;
            else if (current_arg == "--gap_open" && i + 1 < argc) GAP_OPEN = std::stod(argv[++i]); else if (current_arg == "--gap_extend" && i + 1 < argc) GAP_EXTEND = std::stod(argv[++i]);
            else if (current_arg == "--help") { if (::rank == 0) {std::cout << "Usage: ... (Full help message)\n";} MPI_Finalize(); return 0;}
            else { if (::rank == 0) std::cerr << "Unknown option: " << current_arg << "\n"; MPI_Abort(MPI_COMM_WORLD, 1); return 1; }
        }
        ScoreFn current_score_fn = (current_mode == MODE_DNA) ? &edna_score : &blosum62_score;
        if (query_file.empty() || target_file.empty() || align_choice == -1) { if (::rank == 0) std::cerr << "Missing required arguments...\n"; MPI_Abort(MPI_COMM_WORLD, 1); return 1; }
        if (::rank == 0) { try { std::filesystem::create_directories(output_dir); } catch (const std::exception& e) { std::cerr << "Error creating output dir " << output_dir << ": " << e.what() << std::endl; MPI_Abort(MPI_COMM_WORLD, 1); return 1;} }
        std::string seq1_data, seq2_data, h1_data, h2_data;
        if (::rank == 0) {
            try { processFasta(query_file, h1_data, seq1_data); processFasta(target_file, h2_data, seq2_data); } catch (const std::runtime_error& e) { std::cerr << "FASTA error: " << e.what() << std::endl; MPI_Abort(MPI_COMM_WORLD, 1); return 1; }
            if (!fm_idx_path_for_target.empty()) {
                std::ifstream fm_ifs(fm_idx_path_for_target, std::ios::binary);
                if (fm_ifs) { fm_index_target_obj_ptr = std::make_unique<FMIndex>(); if (!fm_index_target_obj_ptr->load(fm_ifs)) { std::cerr << "Rank 0: Error! Failed to load FM-Index from " << fm_idx_path_for_target << "\n"; fm_index_target_obj_ptr.reset(); } else if (verbose) { std::cout << "Rank 0: Loaded FM-Index for target from " << fm_idx_path_for_target << std::endl; }
                } else { std::cerr << "Rank 0: Error! Cannot open FM-Index file: " << fm_idx_path_for_target << "\n"; }
            }
        }
        int s1_len_val = (::rank == 0) ? seq1_data.length() : 0; int s2_len_val = (::rank == 0) ? seq2_data.length() : 0; // Renamed
        MPI_Bcast(&s1_len_val, 1, MPI_INT, 0, MPI_COMM_WORLD); MPI_Bcast(&s2_len_val, 1, MPI_INT, 0, MPI_COMM_WORLD);
        if (::rank != 0) { seq1_data.resize(s1_len_val); seq2_data.resize(s2_len_val); }
        if (s1_len_val > 0) MPI_Bcast(&seq1_data[0], s1_len_val, MPI_CHAR, 0, MPI_COMM_WORLD); 
        if (s2_len_val > 0) MPI_Bcast(&seq2_data[0], s2_len_val, MPI_CHAR, 0, MPI_COMM_WORLD);
        // Broadcast headers if needed by all, for now assume Rank 0 uses them mainly for output.
        int h1_len = (::rank == 0) ? h1_data.length() : 0; int h2_len = (::rank == 0) ? h2_data.length() : 0;
        MPI_Bcast(&h1_len, 1, MPI_INT, 0, MPI_COMM_WORLD); MPI_Bcast(&h2_len, 1, MPI_INT, 0, MPI_COMM_WORLD);
        if (::rank != 0) { h1_data.resize(h1_len); h2_data.resize(h2_len); }
        if (h1_len > 0) MPI_Bcast(&h1_data[0], h1_len, MPI_CHAR, 0, MPI_COMM_WORLD);
        if (h2_len > 0) MPI_Bcast(&h2_data[0], h2_len, MPI_CHAR, 0, MPI_COMM_WORLD);


        const FMIndex* fm_target_raw_ptr_val = fm_index_target_obj_ptr.get(); // Rank 0 has it, others get NULL if not broadcasted.
        // To use FMIndex on all ranks, it must be broadcasted or loaded by each.
        // For this conceptual version, rank 0 does FMIndex queries and broadcasts seeds/anchors.
        // If fm_target_raw_ptr_val is to be used by worker ranks directly in generate_raw_seeds,
        // then the FMIndex object needs to be available on those ranks.
        // This would require broadcasting the entire FMIndex object or having each rank load it.
        // A simpler pattern for MPI is: rank 0 uses its FMIndex to find seeds, then broadcasts seeds.

        if (align_choice == 1) globalalign(seq1_data, seq2_data, h1_data, h2_data, output_dir, current_mode, current_score_fn, fm_target_raw_ptr_val);
        else if (align_choice == 2) localalign(seq1_data, seq2_data, h1_data, h2_data, output_dir, current_mode, current_score_fn, fm_target_raw_ptr_val);
        else if (align_choice == 3) lcs(seq1_data, seq2_data, h1_data, h2_data, output_dir, current_mode, fm_target_raw_ptr_val);
        else if (align_choice == 4) {
            globalalign(seq1_data, seq2_data, h1_data, h2_data, output_dir, current_mode, current_score_fn, fm_target_raw_ptr_val); MPI_Barrier(MPI_COMM_WORLD);
            localalign(seq1_data, seq2_data, h1_data, h2_data, output_dir, current_mode, current_score_fn, fm_target_raw_ptr_val); MPI_Barrier(MPI_COMM_WORLD);
            lcs(seq1_data, seq2_data, h1_data, h2_data, output_dir, current_mode, fm_target_raw_ptr_val);
        } else if (::rank == 0) { std::cerr << "Invalid choice...\n"; }
    } catch (const std::exception& e) { if (::rank == 0) std::cerr << "Runtime Exception: " << e.what() << "\n"; MPI_Abort(MPI_COMM_WORLD, 1);}
      catch (...) { if (::rank == 0) std::cerr << "Unknown exception.\n"; MPI_Abort(MPI_COMM_WORLD, 1);}
    MPI_Finalize(); return 0;
}