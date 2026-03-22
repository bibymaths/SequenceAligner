#include <immintrin.h>
#include <mpi.h>
#include <omp.h>

#include <algorithm>
#include <array>
#include <chrono>
#include <climits>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <memory>
#include <numeric>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>

#ifndef EDNAFULL_MATRIX_DEFINED
#define EDNAFULL_MATRIX_DEFINED
const int EDNAFULL_SIZE                                 = 15;
int       EDNAFULL_matrix[EDNAFULL_SIZE][EDNAFULL_SIZE] = {
    {5, -4, -4, -4, 1, -4, 1, 1, -4, 1, -4, 1, 1, 1, -2},
    {-4, 5, -4, -4, -4, 1, 1, -4, 1, -4, 1, 1, -4, 1, -2},
    {-4, -4, 5, -4, 1, 1, -4, -4, 1, -4, 1, -4, 1, 1, -2},
    {-4, -4, -4, 5, -4, 1, -4, 1, 1, -4, 1, 1, -4, 1, -2},
    {1, -4, 1, -4, -1, -4, -2, -2, -2, -2, -3, -2, -2, -2, -1},
    {-4, 1, 1, 1, -4, -1, -2, -2, -2, -2, -2, -3, -2, -2, -1},
    {1, 1, -4, -4, -2, -2, -1, -4, -2, -4, -2, -2, -2, -2, -1},
    {1, -4, -4, 1, -2, -2, -4, -1, -4, -2, -2, -2, -2, -2, -1},
    {-4, 1, 1, 1, -2, -2, -2, -4, -1, -4, -2, -2, -2, -2, -1},
    {1, -4, 1, -4, -2, -4, -2, -2, -4, -1, -2, -2, -2, -2, -1},
    {-4, 1, 1, 1, -3, -2, -2, -2, -2, -2, -1, -2, -3, -3, -1},
    {1, 1, -4, 1, -2, -3, -2, -2, -2, -2, -2, -1, -3, -3, -1},
    {1, 1, 1, -4, -2, -2, -2, -2, -2, -2, -3, -3, -1, -3, -1},
    {1, 1, 1, 1, -2, -2, -2, -2, -2, -2, -3, -3, -3, -1, -1},
    {-2, -2, -2, -2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1}
};
#endif

#ifndef EBLOSUM62_MATRIX_DEFINED
#define EBLOSUM62_MATRIX_DEFINED
const int EBLOSUM62_SIZE = 24;
int       EBLOSUM62_matrix[EBLOSUM62_SIZE][EBLOSUM62_SIZE] = {
    {4,  -1, -2, -2, 0, -1, -1, 0, -2, -1, -1, -1, -1, -2, -1, 1,  0, -3, -2, 0, -2, -1, -0, -4},
    {-1, 5,  0,  -2, -3, 1,  0,  -2, 0,  -3, -2, 2, -1, -3, -2, -1, -1, -3, -2, -3, -1, 0,  -1, -4},
    {-2, 0,  6,  1, -3, 0,  0,  0,  1, -3, -3, 0, -2, -3, -2, 1, 0,  -4, -2, -3, 3, 0,  -1, -4},
    {-2, -2, 1,  6, -3, 0,  2,  -1, -1, -3, -4, -1, -3, -3, -1, 0, -1, -4, -3, -3, 4,  1,  -1, -4},
    {0,  -3, -3, -3, 9,  -3, -4, -3, -3, -1, -1, -3, -1, -2, -3, -1, -1, -2, -2, -1, -3, -3, -2, -4},
    {-1, 1,  0,  0, -3, 5,  2,  -2, 0, -3, -2, 1, 0,  -3, -1, 0, -1, -2, -1, -2, 0, 3,  -1, -4},
    {-1, 0,  0,  2, -4, 2,  5,  -2, 0, -3, -3, 1, -2, -3, -1, 0, -1, -3, -2, -2, 1, 4,  -1, -4},
    {0,  -2, 0,  -1, -3, -2, -2, 6,  -2, -4, -4, -2, -3, -3, -2, 0,  -2, -2, -3, -3, -1, -2, -1, -4},
    {-2, 0,  1,  -1, -3, 0,  0, -2, 8, -3, -3, -1, -2, -1, -2, -1, -2, -2, 2, -3, 0, 0,  -1, -4},
    {-1, -3, -3, -3, -1, -3, -3, -4, -3, 4,  2,  -3, 1,  0,  -3, -2, -1, -3, -1, 3,  -3, -3, -1, -4},
    {-1, -2, -3, -4, -1, -2, -3, -4, -3, 2,  4,  -2, 2,  0,  -3, -2, -1, -2, -1, 1,  -4, -3, -1, -4},
    {-1, 2,  0,  -1, -3, 1,  1,  -2, -1, -3, -2, 5, -1, -3, -1, 0,  -1, -3, -2, -2, 0,  1,  -1, -4},
    {-1, -1, -2, -3, -1, 0,  -2, -3, -2, 1,  2,  -1, 5,  0,  -2, -1, -1, -1, -1, 1,  -3, -1, -1, -4},
    {-2, -3, -3, -3, -2, -3, -3, -3, -1, 0,  0,  -3, 0,  6,  -4, -2, -2, 1,  3,  -1, -3, -3, -1, -4},
    {-1, -2, -2, -1, -3, -1, -1, -2, -2, -3, -3, -1, -2, -4, 7,  -1, -1, -4, -3, -2, -2, -1, -2, -4},
    {1,  -1, 1,  0, -1, 0,  0,  0,  -1, -2, -2, 0, -1, -2, -1, 4, 1,  -3, -2, -2, 0,  0,  0,  -4},
    {0,  -1, 0,  -1, -1, -1, -1, -2, -2, -1, -1, -1, -1, -2, -1, 1,  5,  -2, -2, 0,  -1, -1, 0,  -4},
    {-3, -3, -4, -4, -2, -2, -3, -2, -2, -3, -2, -3, -1, 1,  -4, -3, -2, 11, 2,  -3, -4, -3, -2, -4},
    {-2, -2, -2, -3, -2, -1, -2, -3, 2,  -1, -1, -2, -1, 3,  -3, -2, -2, 2,  7,  -1, -3, -2, -1, -4},
    {0, -3, -3, -3, -1, -2, -2, -3, -3, 3,  1,  -2, 1, -1, -2, -2, 0,  -3, -1, 4,  -3, -2, -1, -4},
    {-2, -1, 3,  4, -3, 0,  1,  -1, 0, -3, -4, 0, -3, -3, -2, 0, -1, -4, -3, -3, 4, 1,  -1, -4},
    {-1, 0,  0,  1, -3, 3,  4,  -2, 0, -3, -3, 1, -1, -3, -1, 0, -1, -3, -2, -2, 1, 4,  -1, -4},
    {0,  -1, -1, -1, -2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -2, 0,  0,  -2, -1, -1, -1, -1, -1, -4},
    {-4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, 1}
};
#endif

// Global Variables
bool verbose = false;
bool binary  = false;
bool txt     = false;
int  rank_val; // Re-named to avoid conflict

enum ScoreMode { MODE_DNA, MODE_PROTEIN };
using ScoreFn = int (*)(char, char);

double           GAP_OPEN   = -5.0;
double           GAP_EXTEND = -1.0;
static const int LINE_WIDTH = 80;

#define RESET "\033[0m"
#define GREEN "\033[32m"
#define RED "\033[31m"
#define CYAN "\033[36m"

// Suffix Array Construction
std::vector<int> suffix_array_construction(const std::string& s) {
  int              n = s.length();
  std::vector<int> sa_arr(n);
  if (n == 0) return sa_arr;
  std::iota(sa_arr.begin(), sa_arr.end(), 0);

  std::vector<int> rank_arr(n);
  for (int i = 0; i < n; ++i) {
    rank_arr[i] = static_cast<unsigned char>(s[i]);
  }

  std::vector<int> tmp_rank_arr(n);
  for (int k_iter = 1;; k_iter <<= 1) {
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
      int  prev         = sa_arr[i - 1];
      int  curr         = sa_arr[i];
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
class FMIndex {
 public:
  std::string                      text_with_sentinel;
  std::vector<int>                 sa;
  std::string                      bwt;
  char                             sentinel_char;
  std::map<char, int>              C;
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
    if (!this->sa.empty()) {
      for (size_t i = 0; i < this->sa.size(); ++i) {
        if (this->sa[i] == 0) {
          this->bwt += text_with_sentinel.back();
        } else {
          this->bwt += text_with_sentinel[this->sa[i] - 1];
        }
      }
    } else if (text_with_sentinel.length() == 1 &&
               text_with_sentinel[0] == sentinel_char) {
      this->bwt = text_with_sentinel;
    }
    build_c_table();
    build_occ_table();
  }

  FMIndex() : sentinel_char('$') {}

  void build_c_table() {
    std::map<char, int> counts;
    for (char ch : this->bwt) {
      counts[ch]++;
    }
    this->C.clear();
    int total = 0;
    for (const auto& pair : counts) {
      this->C[pair.first] = total;
      total += pair.second;
    }
  }

  void build_occ_table() {
    this->Occ.clear();
    for (const auto& pair : this->C) {
      this->Occ[pair.first].push_back(0);
    }
    for (char ch_in_bwt : this->bwt) {
      if (this->Occ.find(ch_in_bwt) == this->Occ.end()) {
        this->Occ[ch_in_bwt].push_back(0);
      }
    }
    for (size_t i = 0; i < this->bwt.length(); ++i) {
      char current_char_in_bwt = this->bwt[i];
      for (auto& pair : this->Occ) {
        char c          = pair.first;
        int  prev_count = pair.second.back();
        pair.second.push_back(prev_count + (c == current_char_in_bwt ? 1 : 0));
      }
    }
  }

  std::pair<int, int> backward_search(const std::string& pattern) const {
    if (bwt.empty() || pattern.empty()) return {0, 0};
    int l = 0;
    int r = bwt.length();
    for (int i = pattern.length() - 1; i >= 0; --i) {
      char ch     = pattern[i];
      auto c_it   = C.find(ch);
      auto occ_it = Occ.find(ch);
      if (c_it == C.end() || occ_it == Occ.end() || occ_it->second.empty())
        return {0, 0};
      if (static_cast<size_t>(l) >= occ_it->second.size() ||
          static_cast<size_t>(r) >= occ_it->second.size())
        return {0, 0};
      l = c_it->second + occ_it->second[l];
      r = c_it->second + occ_it->second[r];
      if (l >= r) return {0, 0};
    }
    return {l, r};
  }

  std::vector<int> locate(const std::string& pattern) const {
    std::pair<int, int> sa_range = backward_search(pattern);
    std::vector<int>    positions;
    if (sa_range.first < sa_range.second) {
      for (int i = sa_range.first; i < sa_range.second; ++i) {
        if (static_cast<size_t>(i) < sa.size()) positions.push_back(sa[i]);
      }
    }
    std::sort(positions.begin(), positions.end());
    return positions;
  }

  void save(std::ostream& os) const {
    size_t len = text_with_sentinel.length();
    os.write(reinterpret_cast<const char*>(&len), sizeof(len));
    if (len > 0) os.write(text_with_sentinel.data(), len);
    os.write(reinterpret_cast<const char*>(&sentinel_char),
             sizeof(sentinel_char));
    len = sa.size();
    os.write(reinterpret_cast<const char*>(&len), sizeof(len));
    if (len > 0)
      os.write(reinterpret_cast<const char*>(sa.data()), len * sizeof(int));
    len = bwt.length();
    os.write(reinterpret_cast<const char*>(&len), sizeof(len));
    if (len > 0) os.write(bwt.data(), len);
    len = C.size();
    os.write(reinterpret_cast<const char*>(&len), sizeof(len));
    for (const auto& pair : C) {
      os.write(reinterpret_cast<const char*>(&pair.first), sizeof(char));
      os.write(reinterpret_cast<const char*>(&pair.second), sizeof(int));
    }
    len = Occ.size();
    os.write(reinterpret_cast<const char*>(&len), sizeof(len));
    for (const auto& pair : Occ) {
      os.write(reinterpret_cast<const char*>(&pair.first), sizeof(char));
      size_t vec_len = pair.second.size();
      os.write(reinterpret_cast<const char*>(&vec_len), sizeof(vec_len));
      if (vec_len > 0)
        os.write(reinterpret_cast<const char*>(pair.second.data()),
                 vec_len * sizeof(int));
    }
  }

  bool load(std::istream& is) {
    try {
      size_t len;
      is.read(reinterpret_cast<char*>(&len), sizeof(len));
      if (!is || len > 2000000000) return false;
      text_with_sentinel.resize(len);
      if (len > 0) is.read(&text_with_sentinel[0], len);
      is.read(reinterpret_cast<char*>(&sentinel_char), sizeof(sentinel_char));
      is.read(reinterpret_cast<char*>(&len), sizeof(len));
      if (!is || len > 2000000000) return false;
      sa.resize(len);
      if (len > 0)
        is.read(reinterpret_cast<char*>(sa.data()), len * sizeof(int));
      is.read(reinterpret_cast<char*>(&len), sizeof(len));
      if (!is || len > 2000000000) return false;
      bwt.resize(len);
      if (len > 0) is.read(&bwt[0], len);
      is.read(reinterpret_cast<char*>(&len), sizeof(len));
      if (!is) return false;
      C.clear();
      for (size_t i = 0; i < len; ++i) {
        char ch_c;
        int  val_c;
        is.read(reinterpret_cast<char*>(&ch_c), sizeof(char));
        is.read(reinterpret_cast<char*>(&val_c), sizeof(int));
        if (!is) return false;
        C[ch_c] = val_c;
      }
      is.read(reinterpret_cast<char*>(&len), sizeof(len));
      if (!is) return false;
      Occ.clear();
      for (size_t i = 0; i < len; ++i) {
        char   ch_o;
        size_t vec_len;
        is.read(reinterpret_cast<char*>(&ch_o), sizeof(char));
        is.read(reinterpret_cast<char*>(&vec_len), sizeof(vec_len));
        if (!is || vec_len > 2000000000) return false;
        std::vector<int> occ_row(vec_len);
        if (vec_len > 0)
          is.read(reinterpret_cast<char*>(occ_row.data()),
                  vec_len * sizeof(int));
        if (!is && vec_len > 0) return false;
        Occ[ch_o] = occ_row;
      }
    } catch (const std::ios_base::failure& e) {
      if (verbose && rank_val == 0)
        std::cerr << "FMIndex load I/O exception: " << e.what() << std::endl;
      return false;
    } catch (const std::bad_alloc& e) {
      if (verbose && rank_val == 0)
        std::cerr << "FMIndex load memory exception: " << e.what() << std::endl;
      return false;
    }
    return is.good() && !is.eof();
  }
};

// -------- Scoring Lookups & Functions --------
static const std::array<uint8_t, 256> char2idx = []() {
  std::array<uint8_t, 256> m{};
  m.fill(255);
  m[static_cast<unsigned char>('A')] = 0;
  m[static_cast<unsigned char>('C')] = 1;
  m[static_cast<unsigned char>('G')] = 2;
  m[static_cast<unsigned char>('T')] = 3;
  m[static_cast<unsigned char>('U')] = 3;
  m[static_cast<unsigned char>('R')] = 4;
  m[static_cast<unsigned char>('Y')] = 5;
  m[static_cast<unsigned char>('S')] = 6;
  m[static_cast<unsigned char>('W')] = 7;
  m[static_cast<unsigned char>('K')] = 8;
  m[static_cast<unsigned char>('M')] = 9;
  m[static_cast<unsigned char>('B')] = 10;
  m[static_cast<unsigned char>('D')] = 11;
  m[static_cast<unsigned char>('H')] = 12;
  m[static_cast<unsigned char>('V')] = 13;
  m[static_cast<unsigned char>('N')] = 14;
  m[static_cast<unsigned char>('X')] = 14;
  return m;
}();
static const std::array<uint8_t, 256> prot_idx = []() {
  std::array<uint8_t, 256> m{};
  m.fill(255);
  m[static_cast<unsigned char>('A')] = 0;
  m[static_cast<unsigned char>('R')] = 1;
  m[static_cast<unsigned char>('N')] = 2;
  m[static_cast<unsigned char>('D')] = 3;
  m[static_cast<unsigned char>('C')] = 4;
  m[static_cast<unsigned char>('Q')] = 5;
  m[static_cast<unsigned char>('E')] = 6;
  m[static_cast<unsigned char>('G')] = 7;
  m[static_cast<unsigned char>('H')] = 8;
  m[static_cast<unsigned char>('I')] = 9;
  m[static_cast<unsigned char>('L')] = 10;
  m[static_cast<unsigned char>('K')] = 11;
  m[static_cast<unsigned char>('M')] = 12;
  m[static_cast<unsigned char>('F')] = 13;
  m[static_cast<unsigned char>('P')] = 14;
  m[static_cast<unsigned char>('S')] = 15;
  m[static_cast<unsigned char>('T')] = 16;
  m[static_cast<unsigned char>('W')] = 17;
  m[static_cast<unsigned char>('Y')] = 18;
  m[static_cast<unsigned char>('V')] = 19;
  m[static_cast<unsigned char>('B')] = 20;
  m[static_cast<unsigned char>('Z')] = 21;
  m[static_cast<unsigned char>('X')] = 22;
  m[static_cast<unsigned char>('*')] = 23;
  return m;
}();

inline int edna_score(char x, char y) {
  uint8_t ix = char2idx[static_cast<uint8_t>(x)];
  uint8_t iy = char2idx[static_cast<uint8_t>(y)];
  if (ix == 255 || iy == 255)
    throw std::runtime_error(std::string("Invalid DNA code in edna_score: '") +
                             x + "','" + y + "'");
  return EDNAFULL_matrix[ix][iy];
}
inline int blosum62_score(char x, char y) {
  uint8_t ix = prot_idx[static_cast<uint8_t>(x)];
  uint8_t iy = prot_idx[static_cast<uint8_t>(y)];
  if (ix == 255 || iy == 255)
    throw std::runtime_error(
        std::string("Invalid protein code in blosum62_score: '") + x + "','" +
        y + "'");
  return EBLOSUM62_matrix[ix][iy];
}
inline int score(char x, char y, ScoreMode mode) {
  if (mode == MODE_DNA) return edna_score(x, y);
  return blosum62_score(x, y);
}

// -------- Utilities --------
void showProgressBar(int progress, int total) {
  using namespace std::chrono;
  static auto start_time = steady_clock::now();
  if (rank_val != 0) return;

  auto now       = steady_clock::now();
  auto elapsed_s = duration_cast<seconds>(now - start_time).count();
  long eta_s     = 0;
  if (progress > 0 && progress < total) {
    eta_s = elapsed_s * (total - progress) / progress;
  }

  auto format_hms_func = [](long secs_val) {
    std::ostringstream os_formatter;
    long               h_val = secs_val / 3600;
    long               m_val = (secs_val % 3600) / 60;
    long               s_val = secs_val % 60;
    if (h_val > 0) os_formatter << h_val << "h";
    if (m_val > 0 || h_val > 0)
      os_formatter << std::setw(h_val > 0 ? 2 : 1)
                   << std::setfill(h_val > 0 ? '0' : ' ') << m_val << "m";
    os_formatter << std::setw(2) << std::setfill('0') << s_val << "s";
    return os_formatter.str();
  };

  constexpr int bar_width = 50;
  float current_ratio = total > 0 ? static_cast<float>(progress) / total : 0.0f;
  int   filled_pos    = static_cast<int>(bar_width * current_ratio);

  std::cout << "\r[";
  for (int i_bar = 0; i_bar < bar_width; ++i_bar) {
    if (i_bar < filled_pos) std::cout << "=";
    else if (i_bar == filled_pos) std::cout << ">";
    else std::cout << " ";
  }
  std::cout << "] " << std::setw(3) << static_cast<int>(current_ratio * 100.0f)
            << "% " << progress << "/" << total
            << " | Elapsed: " << format_hms_func(elapsed_s)
            << " | ETA: " << format_hms_func(eta_s) << "   ";
  std::cout << std::flush;
  if (progress == total) {
    std::cout << std::endl;
    start_time = steady_clock::now();
  }
}

std::string getAccession(const std::string& header, ScoreMode mode) {
  if (mode == MODE_PROTEIN) {
    size_t firstPipe = header.find('|');
    if (firstPipe != std::string::npos) {
      size_t secondPipe = header.find('|', firstPipe + 1);
      if (secondPipe != std::string::npos) {
        return header.substr(firstPipe + 1, secondPipe - firstPipe - 1);
      }
    }
  }
  std::istringstream iss(header);
  std::string        accession_val;
  iss >> accession_val;
  return accession_val;
}

std::string getGeneSymbol(const std::string& header, ScoreMode mode) {
  if (mode == MODE_DNA) {
    size_t open_p  = header.find('(');
    size_t close_p = (open_p != std::string::npos) ? header.find(')', open_p + 1) : std::string::npos;
    if (open_p != std::string::npos && close_p != std::string::npos && close_p > open_p + 1) {
      return header.substr(open_p + 1, close_p - open_p - 1);
    }
  } else if (mode == MODE_PROTEIN) {
    size_t gn_pos = header.find("GN=");
    if (gn_pos != std::string::npos) {
      size_t start_gn = gn_pos + 3;
      size_t end_gn   = header.find_first_of(" \t", start_gn);
      if (end_gn == std::string::npos) end_gn = header.length();
      if (end_gn > start_gn) return header.substr(start_gn, end_gn - start_gn);
    }
    size_t firstPipe  = header.find('|');
    size_t secondPipe = (firstPipe != std::string::npos) ? header.find('|', firstPipe + 1) : std::string::npos;
    if (secondPipe != std::string::npos) {
      size_t protein_name_start = secondPipe + 1;
      size_t underscore_pos     = header.find('_', protein_name_start);
      if (underscore_pos != std::string::npos && underscore_pos > protein_name_start) {
        return header.substr(protein_name_start, underscore_pos - protein_name_start);
      }
      size_t first_space_after_pipes = header.find(' ', protein_name_start);
      if (first_space_after_pipes != std::string::npos) {
        return header.substr(protein_name_start, first_space_after_pipes - protein_name_start);
      } else {
        return header.substr(protein_name_start);
      }
    }
  }
  return "";
}

void processFasta(const std::string& filename, std::string& header_out, std::string& sequence_out) {
  std::ifstream file_in(filename);
  if (!file_in) {
    throw std::runtime_error("Error: Unable to open FASTA file " + filename);
  }
  std::string line_buffer;
  header_out              = "";
  sequence_out            = "";
  bool first_header_found = false;

  while (getline(file_in, line_buffer)) {
    if (line_buffer.empty()) continue;
    if (line_buffer[0] == '>') {
      if (!first_header_found) {
        header_out = line_buffer.substr(1);
        if (!header_out.empty() && header_out.back() == '\r') {
          header_out.pop_back();
        }
        first_header_found = true;
      } else {
        break;
      }
    } else if (first_header_found) {
      if (!line_buffer.empty() && line_buffer.back() == '\r') {
        line_buffer.pop_back();
      }
      sequence_out += line_buffer;
    }
  }
}

void savePlainAlignment(const std::string& h1, const std::string& h2,
                        const std::string& a1, const std::string& a2,
                        std::ostream& os) {
  os << '>' << h1 << '\n';
  for (size_t i = 0; i < a1.length(); i += LINE_WIDTH) {
    os << a1.substr(i, LINE_WIDTH) << '\n';
  }
  os << '>' << h2 << '\n';
  for (size_t i = 0; i < a2.length(); i += LINE_WIDTH) {
    os << a2.substr(i, LINE_WIDTH) << '\n';
  }
}

void saveLCS(const std::string& id, const std::string& lcs_str_val, std::ostream& os) {
  os << '>' << id << "_LCS_len=" << lcs_str_val.size() << "\n";
  for (size_t i = 0; i < lcs_str_val.size(); i += LINE_WIDTH) {
    os << lcs_str_val.substr(i, LINE_WIDTH) << "\n";
  }
}

void printColoredAlignment(const std::string& seq1_aln,
                           const std::string& seq2_aln,
                           std::ostream&      os = std::cout) {
  size_t aln_len = seq1_aln.length();
  if (aln_len == 0) {
    os << "No alignment to print.\n";
    return;
  }
  if (seq1_aln.length() != seq2_aln.length()) {
    os << "Error: Aligned sequences have different lengths.\n";
    return;
  }

  size_t pos1_count = 0, pos2_count = 0;
  for (size_t i = 0; i < aln_len; i += LINE_WIDTH) {
    size_t end_blk            = std::min(i + LINE_WIDTH, aln_len);
    size_t blk_start_p1       = pos1_count + 1;
    size_t blk_start_p2       = pos2_count + 1;
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

    os << "       ";
    for (size_t j = i; j < end_blk; ++j) {
      if (seq1_aln[j] == seq2_aln[j]) os << "|";
      else if (seq1_aln[j] == '-' || seq2_aln[j] == '-') os << " ";
      else os << ".";
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
  if (!out_file) {
    std::cerr << "Error: Cannot write DP matrix to " << filename << "\n";
    return;
  }
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
  if (!out_file) {
    std::cerr << "Error: Cannot write binary DP matrix to " << filename << "\n";
    return;
  }
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
  if (!out_file) {
    std::cerr << "Error: Cannot open " << filename << "\n";
    return;
  }
  size_t max_cols = 0;
  for (const auto& row_vec : char_matrix) max_cols = std::max(max_cols, row_vec.size());
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
  if (!out_file) {
    std::cerr << "Error: Cannot open " << filename << " for writing.\n";
    return;
  }
  int32_t num_rows = char_matrix.size();
  int32_t num_cols = 0;
  for (const auto& row_vec : char_matrix)
    num_cols = std::max<int32_t>(num_cols, static_cast<int32_t>(row_vec.size()));
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

// -------- Original DP Helper Functions --------
struct AffineDPScores {
  int  s_val = 0;
  int  e_val = 0;
  int  f_val = 0;
  char ptr   = 'X';
};

void initAffineDP(int n_len, std::vector<int>& prev_row_s,
                  std::vector<int>& prev_row_e,
                  std::vector<int>& prev_row_f,
                  bool              isGlobal) {
  prev_row_s.assign(n_len + 1, isGlobal ? (INT_MIN / 2) : 0);
  prev_row_e.assign(n_len + 1, INT_MIN / 2);
  prev_row_f.assign(n_len + 1, INT_MIN / 2);

  if (isGlobal) {
    prev_row_s[0] = 0;
    for (int j = 1; j <= n_len; ++j) {
      prev_row_e[j] = (j == 1) ? (prev_row_s[j - 1] + GAP_OPEN) : (prev_row_e[j - 1] + GAP_EXTEND);
      prev_row_s[j] = prev_row_e[j];
    }
  } else {
    std::fill(prev_row_s.begin(), prev_row_s.end(), 0);
    std::fill(prev_row_e.begin(), prev_row_e.end(), 0);
    std::fill(prev_row_f.begin(), prev_row_f.end(), 0);
  }
}

void computeAffineDPRow(int i_row, const std::string& x_str, const std::string& y_str,
                        std::vector<int>& prev_s_row, std::vector<int>&  prev_e_row, std::vector<int>&  prev_f_row,
                        std::vector<int>& curr_s_row, std::vector<int>&  curr_e_row, std::vector<int>&  curr_f_row,
                        std::vector<char>& curr_trace_s_row, ScoreFn score_fn_local, bool isGlobal) {
  int n_len = y_str.length();
  curr_s_row.assign(n_len + 1, 0);
  curr_e_row.assign(n_len + 1, 0);
  curr_f_row.assign(n_len + 1, 0);
  curr_trace_s_row.assign(n_len + 1, 'X');

  if (isGlobal) {
    int f_i0_open = prev_s_row[0] + GAP_OPEN;
    int f_i0_extend = prev_f_row[0] + GAP_EXTEND;
    curr_f_row[0] = std::max(f_i0_open, f_i0_extend);
    curr_s_row[0] = curr_f_row[0];
    curr_e_row[0] = INT_MIN / 2;
    curr_trace_s_row[0] = (curr_f_row[0] == f_i0_open && curr_f_row[0] >= f_i0_extend) ? 'F' : 'f';
  } else {
    curr_s_row[0]       = 0;
    curr_e_row[0]       = 0;
    curr_f_row[0]       = 0;
    curr_trace_s_row[0] = 'X';
  }

  for (int j_col = 1; j_col <= n_len; ++j_col) {
    int s_diag_prev      = prev_s_row[j_col - 1];
    int e_diag_prev      = prev_e_row[j_col - 1];
    int f_diag_prev      = prev_f_row[j_col - 1];
    int match_pred_score = std::max(s_diag_prev, std::max(e_diag_prev, f_diag_prev));

    int m_val = match_pred_score + score_fn_local(x_str[i_row - 1], y_str[j_col - 1]);

    int e_open_score   = curr_s_row[j_col - 1] + GAP_OPEN;
    int e_extend_score = curr_e_row[j_col - 1] + GAP_EXTEND;
    curr_e_row[j_col]  = std::max(e_open_score, e_extend_score);

    int f_open_score   = prev_s_row[j_col] + GAP_OPEN;
    int f_extend_score = prev_f_row[j_col] + GAP_EXTEND;
    curr_f_row[j_col]  = std::max(f_open_score, f_extend_score);

    if (!isGlobal) {
      curr_e_row[j_col] = std::max(0, curr_e_row[j_col]);
      curr_f_row[j_col] = std::max(0, curr_f_row[j_col]);
      m_val = std::max(0, m_val);
    }

    curr_s_row[j_col] = std::max(m_val, std::max(curr_e_row[j_col], curr_f_row[j_col]));
    if (!isGlobal) curr_s_row[j_col] = std::max(0, curr_s_row[j_col]);

    if (curr_s_row[j_col] > 0 ||
        (isGlobal && (curr_s_row[j_col] == m_val ||
                      curr_s_row[j_col] == curr_e_row[j_col] ||
                      curr_s_row[j_col] == curr_f_row[j_col]))) {
      if (curr_s_row[j_col] == m_val) curr_trace_s_row[j_col] = 'M';
      else if (curr_s_row[j_col] == curr_e_row[j_col]) {
        curr_trace_s_row[j_col] = (curr_e_row[j_col] == e_open_score && curr_e_row[j_col] >= e_extend_score) ? 'E' : 'e';
      } else {
        curr_trace_s_row[j_col] = (curr_f_row[j_col] == f_open_score && curr_f_row[j_col] >= f_extend_score) ? 'F' : 'f';
      }
    } else {
      curr_trace_s_row[j_col] = 'X';
    }
  }
}

struct Loc {
  int score;
  int i;
  int j;
};

// -------- Seed Management --------
struct Seed {
  int  query_pos;
  int  target_pos;
  int  len;
  bool operator<(const Seed& other) const {
    if (query_pos != other.query_pos) return query_pos < other.query_pos;
    if (target_pos != other.target_pos) return target_pos < other.target_pos;
    return len < other.len;
  }
  int query_end() const { return query_pos + len - 1; }
  int target_end() const { return target_pos + len - 1; }
};
struct ChainedSeed {
  std::vector<Seed> seeds;
  double            chain_score;
  int               query_chain_start() const { return seeds.empty() ? -1 : seeds.front().query_pos; }
  int query_chain_end() const { return seeds.empty() ? -1 : seeds.back().query_end(); }
  int target_chain_start() const { return seeds.empty() ? -1 : seeds.front().target_pos; }
  int target_chain_end() const { return seeds.empty() ? -1 : seeds.back().target_end(); }
};

std::vector<Seed> generate_raw_seeds(const std::string& query_seq, const FMIndex& target_fm_index,
                                     int kmer_len, int mpi_rank_val = 0, int mpi_num_procs_val = 1) {
  std::vector<Seed> current_seeds;
  if (kmer_len <= 0) return current_seeds;
  size_t q_len = query_seq.length();
  if (static_cast<size_t>(kmer_len) > q_len) return current_seeds;

  size_t num_kmers_total = q_len - kmer_len + 1;
  size_t kmers_per_rank_chunk = num_kmers_total / mpi_num_procs_val;
  size_t remainder_kmers_val = num_kmers_total % mpi_num_procs_val;
  size_t my_start_kmer_idx_val = mpi_rank_val * kmers_per_rank_chunk + std::min(static_cast<size_t>(mpi_rank_val), remainder_kmers_val);
  size_t my_num_kmers_to_process_val = kmers_per_rank_chunk + (static_cast<size_t>(mpi_rank_val) < remainder_kmers_val ? 1 : 0);
  size_t my_end_kmer_idx_val = my_start_kmer_idx_val + my_num_kmers_to_process_val;

  for (size_t i_kmer = my_start_kmer_idx_val; i_kmer < my_end_kmer_idx_val; ++i_kmer) {
    std::string kmer_str = query_seq.substr(i_kmer, kmer_len);
    std::vector<int> t_pos_list = target_fm_index.locate(kmer_str);
    if (!t_pos_list.empty()) {
      for (int t_p_val : t_pos_list) {
        current_seeds.push_back({(int)i_kmer, t_p_val, kmer_len});
      }
    }
  }
  return current_seeds;
}

ChainedSeed find_best_seed_chain(std::vector<Seed>& seeds_vec, int min_diag_gap_val = 0,
                                 int max_diag_gap_val = 50000, int max_offset_dev_val = 50) {
  if (seeds_vec.empty()) return {};
  std::sort(seeds_vec.begin(), seeds_vec.end());
  int                 n = seeds_vec.size();
  std::vector<double> dp(n, 0.0);
  std::vector<int>    prev(n, -1);
  double best_score = 0.0;
  int    best_idx   = -1;

  for (int i = 0; i < n; ++i) {
    double seed_score = static_cast<double>(seeds_vec[i].len);
    dp[i]             = seed_score;
    for (int j = i - 1; j >= 0; --j) {
      if (seeds_vec[j].query_end() + min_diag_gap_val >= seeds_vec[i].query_pos) continue;
      if (seeds_vec[j].target_end() + min_diag_gap_val >= seeds_vec[i].target_pos) continue;
      int dq = seeds_vec[i].query_pos - seeds_vec[j].query_end() - 1;
      int dt = seeds_vec[i].target_pos - seeds_vec[j].target_end() - 1;
      if (dq < 0 || dt < 0) continue;
      if (dq > max_diag_gap_val || dt > max_diag_gap_val) continue;
      int diag_j = seeds_vec[j].query_pos - seeds_vec[j].target_pos;
      int diag_i = seeds_vec[i].query_pos - seeds_vec[i].target_pos;
      if (std::abs(diag_i - diag_j) > max_offset_dev_val) continue;

      double cost_q   = dq > 0 ? (GAP_OPEN + (dq - 1) * GAP_EXTEND) : 0.0;
      double cost_t   = dt > 0 ? (GAP_OPEN + (dt - 1) * GAP_EXTEND) : 0.0;
      double gap_cost = cost_q + cost_t;
      double cand = dp[j] + seed_score - gap_cost;
      if (cand > dp[i]) { dp[i]   = cand; prev[i] = j; }
    }
    if (dp[i] > best_score) { best_score = dp[i]; best_idx   = i; }
  }
  ChainedSeed chain;
  chain.chain_score = best_score;
  for (int cur = best_idx; cur != -1; cur = prev[cur]) {
    chain.seeds.push_back(seeds_vec[cur]);
  }
  std::reverse(chain.seeds.begin(), chain.seeds.end());
  return chain;
}

// -------- Segment/Window Alignment Helpers & Structs --------
struct AlignmentResult {
  std::string aligned_seq1;
  std::string aligned_seq2;
  int         score            = 0;
  long long   time_ms          = 0;
  int         query_start_orig = -1, query_end_orig = -1;
  int         target_start_orig = -1, target_end_orig = -1;
};
struct LcsSegmentResult {
  std::string lcs_string;
  int         lcs_length = 0;
  std::string gapped_seq1;
  std::string gapped_seq2;
};

AlignmentResult perform_sw_in_window(const std::string& sub1, const std::string& sub2,
                                     ScoreFn sfn, double go, double ge, int q_off, int t_off) {
  AlignmentResult res;
  res.score = 0;
  int m_len = sub1.length();
  int n_len = sub2.length();
  if (m_len == 0 || n_len == 0) return res;
  std::vector<std::vector<int>> s_mat(m_len + 1, std::vector<int>(n_len + 1, 0));
  std::vector<std::vector<int>> e_mat(m_len + 1, std::vector<int>(n_len + 1, 0));
  std::vector<std::vector<int>> f_mat(m_len + 1, std::vector<int>(n_len + 1, 0));
  int max_s_val = 0, max_i_s_val = 0, max_j_s_val = 0;

  for (int i = 1; i <= m_len; ++i) {
    for (int j = 1; j <= n_len; ++j) {
      int m_pred = std::max(s_mat[i - 1][j - 1], std::max(e_mat[i - 1][j - 1], f_mat[i - 1][j - 1]));
      int m_v     = m_pred + sfn(sub1[i - 1], sub2[j - 1]);
      int e_o     = s_mat[i][j - 1] + go;
      int e_e     = e_mat[i][j - 1] + ge;
      e_mat[i][j] = std::max(0, std::max(e_o, e_e));
      int f_o     = s_mat[i - 1][j] + go;
      int f_e     = f_mat[i - 1][j] + ge;
      f_mat[i][j] = std::max(0, std::max(f_o, f_e));
      s_mat[i][j] = std::max(0, std::max(m_v, std::max(e_mat[i][j], f_mat[i][j])));
      if (s_mat[i][j] > max_s_val) {
        max_s_val   = s_mat[i][j];
        max_i_s_val = i;
        max_j_s_val = j;
      }
    }
  }
  res.score = max_s_val;
  if (max_s_val > 0) {
    std::string r_a1, r_a2;
    int         ci = max_i_s_val, cj = max_j_s_val;
    int         current_state = 0;
    int m_check = std::max((ci > 0 && cj > 0 ? s_mat[ci - 1][cj - 1] : INT_MIN / 2),
                  std::max((ci > 0 && cj > 0 ? e_mat[ci - 1][cj - 1] : INT_MIN / 2),
                           (ci > 0 && cj > 0 ? f_mat[ci - 1][cj - 1] : INT_MIN / 2))) +
                  (ci > 0 && cj > 0 ? sfn(sub1[ci - 1], sub2[cj - 1]) : 0);

    if (ci > 0 && cj > 0 && s_mat[ci][cj] == m_check && s_mat[ci][cj] >= e_mat[ci][cj] && s_mat[ci][cj] >= f_mat[ci][cj])
      current_state = 0;
    else if (s_mat[ci][cj] == e_mat[ci][cj] && s_mat[ci][cj] >= f_mat[ci][cj])
      current_state = 1;
    else if (s_mat[ci][cj] == f_mat[ci][cj])
      current_state = 2;

    while (s_mat[ci][cj] > 0 && (ci > 0 || cj > 0)) {
      if (current_state == 0) {
        if (ci <= 0 || cj <= 0) break;
        r_a1 += sub1[ci - 1]; r_a2 += sub2[cj - 1];
        int prev_s = (ci > 1 && cj > 1 ? s_mat[ci - 1][cj - 1] : INT_MIN / 2);
        int prev_e = (ci > 1 && cj > 1 ? e_mat[ci - 1][cj - 1] : INT_MIN / 2);
        int prev_f = (ci > 1 && cj > 1 ? f_mat[ci - 1][cj - 1] : INT_MIN / 2);
        ci--; cj--;
        if (ci < 0 || cj < 0) break;
        if (prev_s >= prev_e && prev_s >= prev_f) current_state = 0;
        else if (prev_e >= prev_f) current_state = 1;
        else current_state = 2;
      } else if (current_state == 1) {
        if (ci <= 0) break;
        r_a1 += sub1[ci - 1]; r_a2 += '-';
        if (e_mat[ci][cj] == (s_mat[ci][cj - 1] + go) && e_mat[ci][cj] >= (e_mat[ci][cj - 1] + ge)) current_state = 0;
        cj--;
      } else {
        if (ci <= 0) break;
        r_a1 += sub1[ci - 1]; r_a2 += '-';
        if (f_mat[ci][cj] == (s_mat[ci - 1][cj] + go) && f_mat[ci][cj] >= (f_mat[ci - 1][cj] + ge)) current_state = 0;
        ci--;
      }
    }
    std::reverse(r_a1.begin(), r_a1.end()); std::reverse(r_a2.begin(), r_a2.end());
    res.aligned_seq1    = r_a1; res.aligned_seq2    = r_a2;
    res.query_end_orig  = q_off + max_i_s_val - 1; res.target_end_orig = t_off + max_j_s_val - 1;
    int q_chars_count   = 0; for (char c_val : res.aligned_seq1) if (c_val != '-') q_chars_count++;
    int t_chars_count = 0; for (char c_val : res.aligned_seq2) if (c_val != '-') t_chars_count++;
    res.query_start_orig  = q_off + (max_i_s_val - q_chars_count);
    res.target_start_orig = t_off + (max_j_s_val - t_chars_count);
  }
  return res;
}

AlignmentResult align_segment_globally(const std::string& seg1, const std::string& seg2,
                                       ScoreFn sfn, double go, double ge) {
  AlignmentResult res;
  int m_len = seg1.length(), n_len = seg2.length();
  if (m_len == 0 && n_len == 0) { res.score = 0; return res; }
  if (m_len == 0) {
    res.aligned_seq1 = std::string(n_len, '-'); res.aligned_seq2 = seg2;
    res.score        = go + (n_len > 1 ? (n_len - 1) * ge : 0); return res;
  }
  if (n_len == 0) {
    res.aligned_seq1 = seg1; res.aligned_seq2 = std::string(m_len, '-');
    res.score        = go + (m_len > 1 ? (m_len - 1) * ge : 0); return res;
  }
  std::vector<int> ps_r(n_len + 1), pe_r(n_len + 1, 0), pf_r(n_len + 1, 0);
  std::vector<std::vector<char>> tr(m_len + 1, std::vector<char>(n_len + 1));
  ps_r[0]  = 0; pe_r[0]  = INT_MIN / 2; pf_r[0]  = INT_MIN / 2; tr[0][0] = 'S';
  for (int j = 1; j <= n_len; ++j) {
    pe_r[j]  = (j == 1 ? (ps_r[j - 1] + go) : (pe_r[j - 1] + ge));
    ps_r[j]  = pe_r[j]; pf_r[j]  = INT_MIN / 2; tr[0][j] = (j == 1 && pe_r[j] == (ps_r[j - 1] + go)) ? 'E' : 'e';
  }
  std::vector<int> cs_r(n_len + 1), ce_r(n_len + 1), cf_r(n_len + 1);
  for (int i = 1; i <= m_len; ++i) {
    int f0_o = ps_r[0] + go, f0_e = pf_r[0] + ge;
    cf_r[0]  = std::max(f0_o, f0_e); cs_r[0]  = cf_r[0]; ce_r[0]  = INT_MIN / 2;
    tr[i][0] = (cf_r[0] == f0_o && cf_r[0] >= f0_e) ? 'F' : 'f';
    for (int j = 1; j <= n_len; ++j) {
      int s_d_val = ps_r[j - 1], e_d_val = pe_r[j - 1], f_d_val = pf_r[j - 1];
      int m_p_val = std::max(s_d_val, std::max(e_d_val, f_d_val));
      int m_v_val = m_p_val + sfn(seg1[i - 1], seg2[j - 1]);
      int e_o_val = cs_r[j - 1] + go, e_e_val = ce_r[j - 1] + ge;
      ce_r[j]     = std::max(e_o_val, e_e_val);
      int f_o_val = ps_r[j] + go, f_e_val = pf_r[j] + ge;
      cf_r[j]     = std::max(f_o_val, f_e_val);
      cs_r[j]     = std::max(m_v_val, std::max(ce_r[j], cf_r[j]));
      if (cs_r[j] == m_v_val) tr[i][j] = 'M';
      else if (cs_r[j] == ce_r[j]) tr[i][j] = (ce_r[j] == e_o_val && ce_r[j] >= e_e_val) ? 'E' : 'e';
      else tr[i][j] = (cf_r[j] == f_o_val && cf_r[j] >= f_e_val) ? 'F' : 'f';
    }
    ps_r.swap(cs_r); pe_r.swap(ce_r); pf_r.swap(cf_r);
  }
  res.score = ps_r[n_len];
  std::string r_a1, r_a2;
  int ci_val = m_len, cj_val = n_len;
  while (ci_val > 0 || cj_val > 0) {
    char move_char = tr[ci_val][cj_val];
    if (move_char == 'M') { r_a1 += seg1[ci_val - 1]; r_a2 += seg2[cj_val - 1]; ci_val--; cj_val--; }
    else if (move_char == 'E' || move_char == 'e') { r_a1 += '-'; r_a2 += seg2[cj_val - 1]; cj_val--; }
    else if (move_char == 'F' || move_char == 'f') { r_a1 += seg1[ci_val - 1]; r_a2 += '-'; ci_val--; }
    else if (ci_val == 0 && cj_val > 0) { r_a1 += '-'; r_a2 += seg2[cj_val - 1]; cj_val--; }
    else if (cj_val == 0 && ci_val > 0) { r_a1 += seg1[ci_val - 1]; r_a2 += '-'; ci_val--; }
    else break;
  }
  std::reverse(r_a1.begin(), r_a1.end()); std::reverse(r_a2.begin(), r_a2.end());
  res.aligned_seq1 = r_a1; res.aligned_seq2 = r_a2;
  return res;
}

LcsSegmentResult compute_lcs_for_segment(const std::string& seg1, const std::string& seg2) {
  LcsSegmentResult res;
  int m_len = seg1.length(), n_len = seg2.length();
  if (m_len == 0 || n_len == 0) {
    res.lcs_length = 0; res.lcs_string = "";
    if (m_len == 0 && n_len > 0) { res.gapped_seq1 = std::string(n_len, '-'); res.gapped_seq2 = seg2; }
    else if (n_len == 0 && m_len > 0) { res.gapped_seq1 = seg1; res.gapped_seq2 = std::string(m_len, '-'); }
    return res;
  }
  std::vector<std::vector<int>>  L_mat(m_len + 1, std::vector<int>(n_len + 1, 0));
  std::vector<std::vector<char>> B_mat(m_len + 1, std::vector<char>(n_len + 1, ' '));
  for (int i = 1; i <= m_len; ++i) {
    for (int j = 1; j <= n_len; ++j) {
      if (seg1[i - 1] == seg2[j - 1]) { L_mat[i][j] = L_mat[i - 1][j - 1] + 1; B_mat[i][j] = 'D'; }
      else if (L_mat[i - 1][j] >= L_mat[i][j - 1]) { L_mat[i][j] = L_mat[i - 1][j]; B_mat[i][j] = 'U'; }
      else { L_mat[i][j] = L_mat[i][j - 1]; B_mat[i][j] = 'L'; }
    }
  }
  res.lcs_length = L_mat[m_len][n_len];
  std::string l_rev, g1_rev, g2_rev;
  int ci_lcs = m_len, cj_lcs = n_len;
  while (ci_lcs > 0 || cj_lcs > 0) {
    if (ci_lcs > 0 && cj_lcs > 0 && B_mat[ci_lcs][cj_lcs] == 'D') {
      l_rev += seg1[ci_lcs - 1]; g1_rev += seg1[ci_lcs - 1]; g2_rev += seg2[cj_lcs - 1]; ci_lcs--; cj_lcs--;
    } else if (ci_lcs > 0 && (cj_lcs == 0 || B_mat[ci_lcs][cj_lcs] == 'U')) {
      g1_rev += seg1[ci_lcs - 1]; g2_rev += '-'; ci_lcs--;
    } else if (cj_lcs > 0 && (ci_lcs == 0 || B_mat[ci_lcs][cj_lcs] == 'L')) {
      g1_rev += '-'; g2_rev += seg2[cj_lcs - 1]; cj_lcs--;
    } else break;
  }
  std::reverse(l_rev.begin(), l_rev.end()); std::reverse(g1_rev.begin(), g1_rev.end()); std::reverse(g2_rev.begin(), g2_rev.end());
  res.lcs_string = l_rev; res.gapped_seq1 = g1_rev; res.gapped_seq2 = g2_rev;
  return res;
}

// ===========================
// ADD THESE HELPERS
// ===========================

static inline void get_row_partition(int total_rows, int world_size, int world_rank,
                                     int& start_row_0based, int& row_count) {
  int base = total_rows / world_size;
  int rem  = total_rows % world_size;
  row_count = base + (world_rank < rem ? 1 : 0);
  start_row_0based = world_rank * base + std::min(world_rank, rem);
}

static std::vector<std::vector<int>> gather_score_matrix_rows(
    const std::vector<std::vector<int>>& local_block_with_boundary,
    int total_rows_without_row0, int n_cols_minus1, bool is_global, MPI_Comm comm) {
  int world_rank = 0, world_size = 1;
  MPI_Comm_rank(comm, &world_rank);
  MPI_Comm_size(comm, &world_size);

  const int ncols = n_cols_minus1 + 1;
  int start_row = 0, local_rows = 0;
  get_row_partition(total_rows_without_row0, world_size, world_rank, start_row, local_rows);

  std::vector<int> sendbuf(local_rows * ncols, 0);
  for (int i = 0; i < local_rows; ++i) {
    std::memcpy(&sendbuf[i * ncols], local_block_with_boundary[i + 1].data(), ncols * sizeof(int));
  }

  std::vector<int> recvcounts, displs, recvbuf;
  if (world_rank == 0) {
    recvcounts.resize(world_size, 0); displs.resize(world_size, 0);
    int total_elems = 0;
    for (int r = 0; r < world_size; ++r) {
      int rs = 0, rr = 0;
      get_row_partition(total_rows_without_row0, world_size, r, rs, rr);
      recvcounts[r] = rr * ncols; displs[r] = total_elems;
      total_elems += recvcounts[r];
    }
    recvbuf.resize(total_elems, 0);
  }

  MPI_Gatherv(sendbuf.data(), static_cast<int>(sendbuf.size()), MPI_INT,
              world_rank == 0 ? recvbuf.data() : nullptr,
              world_rank == 0 ? recvcounts.data() : nullptr,
              world_rank == 0 ? displs.data() : nullptr, MPI_INT, 0, comm);

  std::vector<std::vector<int>> full;
  if (world_rank == 0) {
    full.assign(total_rows_without_row0 + 1, std::vector<int>(ncols, 0));
    if (is_global) {
      std::vector<int> row0_s, row0_e, row0_f;
      initAffineDP(n_cols_minus1, row0_s, row0_e, row0_f, true);
      full[0] = row0_s;
    }
    for (int r = 0; r < world_size; ++r) {
      int rs = 0, rr = 0;
      get_row_partition(total_rows_without_row0, world_size, r, rs, rr);
      for (int i = 0; i < rr; ++i) {
        std::memcpy(full[rs + 1 + i].data(), &recvbuf[displs[r] + i * ncols], ncols * sizeof(int));
      }
    }
  }
  return full;
}

static void save_path_file(const std::string& filename, const std::vector<std::pair<int, int>>& path) {
  std::ofstream pf(filename);
  if (!pf) { std::cerr << "Error: Cannot open " << filename << " for writing.\n"; return; }
  for (const auto& p : path) pf << p.first << " " << p.second << "\n";
}

static AlignmentResult traceback_global_from_full(
    const std::string& x, const std::string& y, ScoreFn score_fn,
    std::vector<std::pair<int, int>>& global_path_out) {
  const int m = static_cast<int>(x.size()), n = static_cast<int>(y.size());
  std::vector<int> prev_s, prev_e, prev_f, curr_s, curr_e, curr_f;
  std::vector<char> curr_tr;
  std::vector<std::vector<char>> trace(m + 1, std::vector<char>(n + 1, 'X'));

  prev_s.assign(n + 1, 0); prev_e.assign(n + 1, INT_MIN / 2); prev_f.assign(n + 1, INT_MIN / 2);
  trace[0][0] = 'S';
  for (int j = 1; j <= n; ++j) {
    prev_e[j] = (j == 1) ? (prev_s[j - 1] + GAP_OPEN) : (prev_e[j - 1] + GAP_EXTEND);
    prev_s[j] = prev_e[j]; prev_f[j] = INT_MIN / 2;
    trace[0][j] = (j == 1 && prev_e[j] == prev_s[j - 1] + GAP_OPEN) ? 'E' : 'e';
  }

  for (int i = 1; i <= m; ++i) {
    computeAffineDPRow(i, x, y, prev_s, prev_e, prev_f, curr_s, curr_e, curr_f, curr_tr, score_fn, true);
    for (int j = 0; j <= n; ++j) trace[i][j] = curr_tr[j];
    prev_s.swap(curr_s); prev_e.swap(curr_e); prev_f.swap(curr_f);
  }

  AlignmentResult res; res.score = prev_s[n];
  int ci = m, cj = n;
  global_path_out.clear(); global_path_out.emplace_back(cj, ci);
  std::string ax, ay;

  while (ci > 0 || cj > 0) {
    char t = trace[ci][cj];
    if (t == 'M') { ax += x[ci - 1]; ay += y[cj - 1]; --ci; --cj; }
    else if (t == 'F' || t == 'f') { ax += x[ci - 1]; ay += '-'; --ci; }
    else if (t == 'E' || t == 'e') { ax += '-'; ay += y[cj - 1]; --cj; }
    else {
      if (ci == 0 && cj > 0) { ax += '-'; ay += y[cj - 1]; --cj; }
      else if (cj == 0 && ci > 0) { ax += x[ci - 1]; ay += '-'; --ci; }
      else break;
    }
    global_path_out.emplace_back(cj, ci);
  }
  std::reverse(ax.begin(), ax.end()); std::reverse(ay.begin(), ay.end());
  res.aligned_seq1 = ax; res.aligned_seq2 = ay; return res;
}

static AlignmentResult traceback_local_from_full(
    const std::string& x, const std::string& y, ScoreFn score_fn,
    int best_i, int best_j, std::vector<std::pair<int, int>>& local_path_out) {
  const int m = static_cast<int>(x.size()), n = static_cast<int>(y.size());
  std::vector<int> prev_s, prev_e, prev_f, curr_s, curr_e, curr_f;
  std::vector<char> curr_tr;

  prev_s.assign(n + 1, 0); prev_e.assign(n + 1, 0); prev_f.assign(n + 1, 0);
  std::vector<std::vector<int>> S(m + 1, std::vector<int>(n + 1, 0));
  std::vector<std::vector<char>> trace(m + 1, std::vector<char>(n + 1, 'X'));

  for (int i = 1; i <= m; ++i) {
    computeAffineDPRow(i, x, y, prev_s, prev_e, prev_f, curr_s, curr_e, curr_f, curr_tr, score_fn, false);
    for (int j = 0; j <= n; ++j) { S[i][j] = curr_s[j]; trace[i][j] = curr_tr[j]; }
    prev_s.swap(curr_s); prev_e.swap(curr_e); prev_f.swap(curr_f);
  }

  AlignmentResult res; res.score = S[best_i][best_j];
  int ci = best_i, cj = best_j; std::string ax, ay;
  local_path_out.clear(); local_path_out.emplace_back(cj, ci);

  while (ci > 0 || cj > 0) {
    if (S[ci][cj] == 0 && trace[ci][cj] == 'X') break;
    char t = trace[ci][cj];
    if (t == 'M') { ax += x[ci - 1]; ay += y[cj - 1]; --ci; --cj; }
    else if (t == 'F' || t == 'f') { ax += x[ci - 1]; ay += '-'; --ci; }
    else if (t == 'E' || t == 'e') { ax += '-'; ay += y[cj - 1]; --cj; }
    else break;
    local_path_out.emplace_back(cj, ci);
  }

  std::reverse(ax.begin(), ax.end()); std::reverse(ay.begin(), ay.end()); std::reverse(local_path_out.begin(), local_path_out.end());
  res.aligned_seq1 = ax; res.aligned_seq2 = ay; return res;
}

static LcsSegmentResult traceback_lcs_from_full(
    const std::string& x, const std::string& y, const std::vector<std::vector<int>>& L,
    std::vector<std::vector<char>>& B, std::vector<std::pair<int, int>>& path_out) {
  const int m = static_cast<int>(x.size()), n = static_cast<int>(y.size());
  B.assign(m + 1, std::vector<char>(n + 1, ' '));
  for (int i = 1; i <= m; ++i) {
    for (int j = 1; j <= n; ++j) {
      if (x[i - 1] == y[j - 1] && L[i][j] == L[i - 1][j - 1] + 1) B[i][j] = 'D';
      else if (L[i - 1][j] >= L[i][j - 1]) B[i][j] = 'U';
      else B[i][j] = 'L';
    }
  }

  LcsSegmentResult res; res.lcs_length = L[m][n];
  int ci = m, cj = n; std::string lcs_rev, ax_rev, ay_rev; path_out.clear();

  while (ci > 0 && cj > 0) {
    path_out.emplace_back(cj, ci);
    if (B[ci][cj] == 'D') { lcs_rev += x[ci - 1]; ax_rev += x[ci - 1]; ay_rev += y[cj - 1]; --ci; --cj; }
    else if (B[ci][cj] == 'U') { ax_rev += x[ci - 1]; ay_rev += '-'; --ci; }
    else { ax_rev += '-'; ay_rev += y[cj - 1]; --cj; }
  }
  while (ci > 0) { path_out.emplace_back(cj, ci); ax_rev += x[ci - 1]; ay_rev += '-'; --ci; }
  while (cj > 0) { path_out.emplace_back(cj, ci); ax_rev += '-'; ay_rev += y[cj - 1]; --cj; }
  path_out.emplace_back(0, 0);

  std::reverse(lcs_rev.begin(), lcs_rev.end()); std::reverse(ax_rev.begin(), ax_rev.end());
  std::reverse(ay_rev.begin(), ay_rev.end()); std::reverse(path_out.begin(), path_out.end());
  res.lcs_string = lcs_rev; res.gapped_seq1 = ax_rev; res.gapped_seq2 = ay_rev; return res;
}

// ===========================
// ALIGNMENT FUNCTIONS
// ===========================

void globalalign(const std::string& x_orig, const std::string& y_orig,
                 const std::string& header1, const std::string& header2,
                 const std::string& outdir, ScoreMode mode_val,
                 ScoreFn score_fn_val, const FMIndex* target_fm_idx) {
  const int m = static_cast<int>(x_orig.size()), n = static_cast<int>(y_orig.size());
  int world_rank = rank_val, world_size = 1; // using global rank_val variable
  MPI_Comm_size(MPI_COMM_WORLD, &world_size);

  auto t_start = std::chrono::high_resolution_clock::now();

  bool use_anchoring = false;
  ChainedSeed best_chain;

  if (world_rank == 0 && target_fm_idx != nullptr) {
    int k_anchor = std::min(8, std::min(m / 12, n / 12));
    if (std::min(m, n) < k_anchor) k_anchor = std::min(m, n);

    if (k_anchor > 0) {
      std::vector<Seed> raw_seeds = generate_raw_seeds(x_orig, *target_fm_idx, k_anchor, 0, 1);
      if (!raw_seeds.empty()) {
        best_chain = find_best_seed_chain(raw_seeds, 1);
        use_anchoring = !best_chain.seeds.empty();
      }
    }
  }

  int anchor_count = (world_rank == 0 && use_anchoring) ? static_cast<int>(best_chain.seeds.size()) : 0;
  MPI_Bcast(&anchor_count, 1, MPI_INT, 0, MPI_COMM_WORLD);
  use_anchoring = anchor_count > 0;

  if (use_anchoring && world_rank != 0) best_chain.seeds.resize(anchor_count);
  if (use_anchoring) {
    MPI_Bcast(best_chain.seeds.data(), anchor_count * static_cast<int>(sizeof(Seed)), MPI_BYTE, 0, MPI_COMM_WORLD);
  }

  AlignmentResult final_aln;
  std::vector<std::pair<int, int>> global_path;

  if (use_anchoring) {
    if (world_rank == 0 && verbose) std::cout << "Global alignment: FM-index anchors found (" << anchor_count << "). Using anchored segmentation.\n";

    std::vector<std::pair<std::string, std::string>> segments;
    int cx = 0, cy = 0;
    for (const auto& anc : best_chain.seeds) {
      segments.push_back({ x_orig.substr(cx, anc.query_pos - cx), y_orig.substr(cy, anc.target_pos - cy) });
      cx = anc.query_pos + anc.len; cy = anc.target_pos + anc.len;
    }
    segments.push_back({x_orig.substr(cx), y_orig.substr(cy)});

    std::vector<AlignmentResult> seg_results(segments.size());
    for (size_t i = 0; i < segments.size(); ++i) {
      if (static_cast<int>(i % world_size) == world_rank) {
        seg_results[i] = align_segment_globally(segments[i].first, segments[i].second, score_fn_val, GAP_OPEN, GAP_EXTEND);
      }
    }

    if (world_rank == 0) {
      for (size_t i = 0; i < segments.size(); ++i) {
        int owner = static_cast<int>(i % world_size);
        if (owner != 0) {
          MPI_Status st;
          MPI_Recv(&seg_results[i].score, 1, MPI_INT, owner, static_cast<int>(i) * 10 + 0, MPI_COMM_WORLD, &st);
          int len1 = 0, len2 = 0;
          MPI_Recv(&len1, 1, MPI_INT, owner, static_cast<int>(i) * 10 + 1, MPI_COMM_WORLD, &st);
          seg_results[i].aligned_seq1.resize(len1);
          if (len1 > 0) MPI_Recv(&seg_results[i].aligned_seq1[0], len1, MPI_CHAR, owner, static_cast<int>(i) * 10 + 2, MPI_COMM_WORLD, &st);
          MPI_Recv(&len2, 1, MPI_INT, owner, static_cast<int>(i) * 10 + 3, MPI_COMM_WORLD, &st);
          seg_results[i].aligned_seq2.resize(len2);
          if (len2 > 0) MPI_Recv(&seg_results[i].aligned_seq2[0], len2, MPI_CHAR, owner, static_cast<int>(i) * 10 + 4, MPI_COMM_WORLD, &st);
        }
      }

      for (size_t i = 0; i < best_chain.seeds.size(); ++i) {
        final_aln.aligned_seq1 += seg_results[i].aligned_seq1; final_aln.aligned_seq2 += seg_results[i].aligned_seq2;
        final_aln.score += seg_results[i].score;
        const auto& anc = best_chain.seeds[i];
        std::string exact = x_orig.substr(anc.query_pos, anc.len);
        final_aln.aligned_seq1 += exact; final_aln.aligned_seq2 += exact;
        for (char c : exact) final_aln.score += score_fn_val(c, c);
      }
      final_aln.aligned_seq1 += seg_results.back().aligned_seq1; final_aln.aligned_seq2 += seg_results.back().aligned_seq2;
      final_aln.score += seg_results.back().score;
    } else {
      for (size_t i = 0; i < segments.size(); ++i) {
        if (static_cast<int>(i % world_size) == world_rank) {
          MPI_Send(&seg_results[i].score, 1, MPI_INT, 0, static_cast<int>(i) * 10 + 0, MPI_COMM_WORLD);
          int len1 = static_cast<int>(seg_results[i].aligned_seq1.size()), len2 = static_cast<int>(seg_results[i].aligned_seq2.size());
          MPI_Send(&len1, 1, MPI_INT, 0, static_cast<int>(i) * 10 + 1, MPI_COMM_WORLD);
          if (len1 > 0) MPI_Send(seg_results[i].aligned_seq1.data(), len1, MPI_CHAR, 0, static_cast<int>(i) * 10 + 2, MPI_COMM_WORLD);
          MPI_Send(&len2, 1, MPI_INT, 0, static_cast<int>(i) * 10 + 3, MPI_COMM_WORLD);
          if (len2 > 0) MPI_Send(seg_results[i].aligned_seq2.data(), len2, MPI_CHAR, 0, static_cast<int>(i) * 10 + 4, MPI_COMM_WORLD);
        }
      }
    }
  } else {
    if (world_rank == 0 && verbose) std::cout << "Global alignment: FM-index anchoring unavailable/failed. Falling back to MPI full DP.\n";
    int start_row = 0, local_rows = 0;
    get_row_partition(m, world_size, world_rank, start_row, local_rows);
    std::vector<std::vector<int>> s_block(local_rows + 1, std::vector<int>(n + 1, 0)), e_block(local_rows + 1, std::vector<int>(n + 1, 0)), f_block(local_rows + 1, std::vector<int>(n + 1, 0));

    if (world_rank == 0) {
      std::vector<int> row0s, row0e, row0f;
      initAffineDP(n, row0s, row0e, row0f, true);
      s_block[0] = row0s; e_block[0] = row0e; f_block[0] = row0f;
    } else {
      MPI_Recv(s_block[0].data(), n + 1, MPI_INT, world_rank - 1, 100, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      MPI_Recv(e_block[0].data(), n + 1, MPI_INT, world_rank - 1, 101, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      MPI_Recv(f_block[0].data(), n + 1, MPI_INT, world_rank - 1, 102, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }

    std::vector<int> curr_s, curr_e, curr_f; std::vector<char> curr_tr;
    for (int ii = 1; ii <= local_rows; ++ii) {
      int global_i = start_row + ii;
      computeAffineDPRow(global_i, x_orig, y_orig, s_block[ii - 1], e_block[ii - 1], f_block[ii - 1], curr_s, curr_e, curr_f, curr_tr, score_fn_val, true);
      s_block[ii] = curr_s; e_block[ii] = curr_e; f_block[ii] = curr_f;
      if (verbose && world_rank == 0 && (global_i % 100 == 0 || global_i == m)) showProgressBar(global_i, m);
    }

    if (world_rank + 1 < world_size) {
      MPI_Send(s_block[local_rows].data(), n + 1, MPI_INT, world_rank + 1, 100, MPI_COMM_WORLD);
      MPI_Send(e_block[local_rows].data(), n + 1, MPI_INT, world_rank + 1, 101, MPI_COMM_WORLD);
      MPI_Send(f_block[local_rows].data(), n + 1, MPI_INT, world_rank + 1, 102, MPI_COMM_WORLD);
    }

    std::vector<std::vector<int>> fullS = gather_score_matrix_rows(s_block, m, n, true, MPI_COMM_WORLD);

    if (world_rank == 0) {
      if (binary) writeDPMatrix(fullS, outdir + "/global_dp_matrix.bin");
      else if (txt) writeRawDPMatrix(fullS, outdir + "/global_dp_matrix.txt");
      final_aln = traceback_global_from_full(x_orig, y_orig, score_fn_val, global_path);
    }
  }

  auto t_end = std::chrono::high_resolution_clock::now();
  final_aln.time_ms = std::chrono::duration_cast<std::chrono::milliseconds>(t_end - t_start).count();

  if (world_rank == 0) {
    if (use_anchoring) {
      if (txt || binary) std::cout << "\nNotice: DP Matrix skipped during anchored Global Alignment.\n";
      int tr_cx = m, tr_cy = n;
      global_path.emplace_back(tr_cy, tr_cx);
      for (int i = final_aln.aligned_seq1.size() - 1; i >= 0; --i) {
        if (final_aln.aligned_seq1[i] != '-') tr_cx--;
        if (final_aln.aligned_seq2[i] != '-') tr_cy--;
        global_path.emplace_back(tr_cy, tr_cx);
      }
      std::reverse(global_path.begin(), global_path.end());
    }

    size_t total = final_aln.aligned_seq1.size(), gaps = 0, matches = 0;
    for (size_t i = 0; i < total; ++i) {
      if (final_aln.aligned_seq1[i] == '-' || final_aln.aligned_seq2[i] == '-') ++gaps;
      else if (final_aln.aligned_seq1[i] == final_aln.aligned_seq2[i]) ++matches;
    }
    double identity = total > 0 ? static_cast<double>(matches) / total : 0.0;
    double coverage = total > 0 ? static_cast<double>(total - gaps) / total : 0.0;

    std::string acc1 = getAccession(header1, mode_val), acc2 = getAccession(header2, mode_val);
    std::string gene1 = getGeneSymbol(header1, mode_val), gene2 = getGeneSymbol(header2, mode_val);

    save_path_file(outdir + "/global_path.txt", global_path);

    if (verbose) {
      std::cout << "\n\nGlobal Alignment Score: " << final_aln.score << "\n";
      std::cout << "Matches: " << matches << " | Gaps: " << gaps << " | Total: " << total << "\n";
      std::cout << "Identity: " << identity * 100.0 << "% | Coverage: " << coverage * 100.0 << "%\n";
      std::cout << "Time: " << final_aln.time_ms << " ms\n\n";
      printColoredAlignment(final_aln.aligned_seq1, final_aln.aligned_seq2);
    }

    std::ofstream fa(outdir + "/global_alignment.fasta");
    if (fa) savePlainAlignment(acc1, acc2, final_aln.aligned_seq1, final_aln.aligned_seq2, fa);

    std::ofstream js(outdir + "/global_stats.json");
    if (js) {
      js << std::fixed << std::setprecision(6) << "{\n"
         << "  \"method\": \"global\",\n"
         << "  \"score\": " << final_aln.score << ",\n"
         << "  \"matches\": " << matches << ",\n"
         << "  \"gaps\": " << gaps << ",\n"
         << "  \"total\": " << total << ",\n"
         << "  \"identity\": " << identity << ",\n"
         << "  \"coverage\": " << coverage << ",\n"
         << "  \"time_ms\": " << final_aln.time_ms << ",\n"
         << R"(  "query": ")" << acc1 << "\",\n"
         << R"(  "target": ")" << acc2 << "\",\n"
         << R"(  "queryid": ")" << gene1 << "\",\n"
         << R"(  "targetid": ")" << gene2 << "\",\n"
         << "  \"query_length\": " << m << ",\n"
         << "  \"target_length\": " << n << "\n"
         << "}\n";
    }
  }
  MPI_Barrier(MPI_COMM_WORLD);
}

void localalign(const std::string& x, const std::string& y,
                const std::string& h1, const std::string& h2,
                const std::string& odir, ScoreMode mval, ScoreFn sfn_val,
                const FMIndex* tfm_idx) {
  const int m = static_cast<int>(x.size()), n = static_cast<int>(y.size());
  int world_rank = rank_val, world_size = 1;
  MPI_Comm_size(MPI_COMM_WORLD, &world_size);
  auto t_start = std::chrono::high_resolution_clock::now();

  bool use_fmindex = false; AlignmentResult best_seed_extend; best_seed_extend.score = 0;

  if (world_rank == 0 && tfm_idx != nullptr) {
    int k = std::min(11, std::min(m / 20, n / 20));
    if (std::min(m, n) < k) k = std::min(m, n);

    if (k > 0) {
      std::vector<Seed> all_seeds = generate_raw_seeds(x, *tfm_idx, k, 0, 1);
      if (!all_seeds.empty()) {
        use_fmindex = true;
        int total_seeds = static_cast<int>(all_seeds.size());
        MPI_Bcast(&total_seeds, 1, MPI_INT, 0, MPI_COMM_WORLD);
        MPI_Bcast(all_seeds.data(), total_seeds * static_cast<int>(sizeof(Seed)), MPI_BYTE, 0, MPI_COMM_WORLD);

        int s_start = (total_seeds / world_size * world_rank) + std::min(world_rank, total_seeds % world_size);
        int s_end   = (total_seeds / world_size * (world_rank + 1)) + std::min(world_rank + 1, total_seeds % world_size);
        AlignmentResult local_best; local_best.score = 0;

        for (int si = s_start; si < s_end; ++si) {
          const auto& seed = all_seeds[si];
          int win_q = std::max(100, seed.len * 3), win_t = std::max(100, seed.len * 3);
          int qws = std::max(0, seed.query_pos - win_q), qwe = std::min(m, seed.query_pos + seed.len + win_q);
          int tws = std::max(0, seed.target_pos - win_t), twe = std::min(n, seed.target_pos + seed.len + win_t);
          AlignmentResult cand = perform_sw_in_window(x.substr(qws, qwe - qws), y.substr(tws, twe - tws), sfn_val, GAP_OPEN, GAP_EXTEND, qws, tws);
          if (cand.score > local_best.score) local_best = cand;
        }

        struct { int score; int rank; } mine{local_best.score, world_rank}, best{0, 0};
        MPI_Allreduce(&mine, &best, 1, MPI_2INT, MPI_MAXLOC, MPI_COMM_WORLD);

        if (best.score > 0) {
          int best_rank = best.rank, len1 = 0, len2 = 0;
          if (world_rank == best_rank) {
            best_seed_extend = local_best; best_seed_extend.score = best.score;
            len1 = static_cast<int>(local_best.aligned_seq1.size()); len2 = static_cast<int>(local_best.aligned_seq2.size());
          }
          MPI_Bcast(&len1, 1, MPI_INT, best_rank, MPI_COMM_WORLD);
          MPI_Bcast(&len2, 1, MPI_INT, best_rank, MPI_COMM_WORLD);
          MPI_Bcast(&best_seed_extend.query_start_orig, 1, MPI_INT, best_rank, MPI_COMM_WORLD);
          MPI_Bcast(&best_seed_extend.query_end_orig,   1, MPI_INT, best_rank, MPI_COMM_WORLD);
          MPI_Bcast(&best_seed_extend.target_start_orig, 1, MPI_INT, best_rank, MPI_COMM_WORLD);
          MPI_Bcast(&best_seed_extend.target_end_orig,   1, MPI_INT, best_rank, MPI_COMM_WORLD);

          if (world_rank != best_rank) { best_seed_extend.aligned_seq1.resize(len1); best_seed_extend.aligned_seq2.resize(len2); }
          if (len1 > 0) MPI_Bcast(&best_seed_extend.aligned_seq1[0], len1, MPI_CHAR, best_rank, MPI_COMM_WORLD);
          if (len2 > 0) MPI_Bcast(&best_seed_extend.aligned_seq2[0], len2, MPI_CHAR, best_rank, MPI_COMM_WORLD);
        }
      } else { int zero = 0; MPI_Bcast(&zero, 1, MPI_INT, 0, MPI_COMM_WORLD); }
    } else { int zero = 0; MPI_Bcast(&zero, 1, MPI_INT, 0, MPI_COMM_WORLD); }
  } else {
    int total_seeds = 0;
    MPI_Bcast(&total_seeds, 1, MPI_INT, 0, MPI_COMM_WORLD);
    if (total_seeds > 0) {
        // ... (worker ranks logic matches above if rank != 0) ...
        // In the interest of ensuring the worker execution functions when seeded
        use_fmindex = true;
        std::vector<Seed> all_seeds(total_seeds);
        MPI_Bcast(all_seeds.data(), total_seeds * static_cast<int>(sizeof(Seed)), MPI_BYTE, 0, MPI_COMM_WORLD);

        int s_start = (total_seeds / world_size * world_rank) + std::min(world_rank, total_seeds % world_size);
        int s_end   = (total_seeds / world_size * (world_rank + 1)) + std::min(world_rank + 1, total_seeds % world_size);

        AlignmentResult local_best; local_best.score = 0;

        for (int si = s_start; si < s_end; ++si) {
          const auto& seed = all_seeds[si];
          int win_q = std::max(100, seed.len * 3), win_t = std::max(100, seed.len * 3);
          int qws = std::max(0, seed.query_pos - win_q), qwe = std::min(m, seed.query_pos + seed.len + win_q);
          int tws = std::max(0, seed.target_pos - win_t), twe = std::min(n, seed.target_pos + seed.len + win_t);
          AlignmentResult cand = perform_sw_in_window(x.substr(qws, qwe - qws), y.substr(tws, twe - tws), sfn_val, GAP_OPEN, GAP_EXTEND, qws, tws);
          if (cand.score > local_best.score) local_best = cand;
        }

        struct { int score; int rank; } mine{local_best.score, world_rank}, best{0, 0};
        MPI_Allreduce(&mine, &best, 1, MPI_2INT, MPI_MAXLOC, MPI_COMM_WORLD);

        if (best.score > 0) {
          int best_rank = best.rank, len1 = 0, len2 = 0;
          if (world_rank == best_rank) {
            best_seed_extend = local_best; best_seed_extend.score = best.score;
            len1 = static_cast<int>(local_best.aligned_seq1.size()); len2 = static_cast<int>(local_best.aligned_seq2.size());
          }
          MPI_Bcast(&len1, 1, MPI_INT, best_rank, MPI_COMM_WORLD);
          MPI_Bcast(&len2, 1, MPI_INT, best_rank, MPI_COMM_WORLD);
          MPI_Bcast(&best_seed_extend.query_start_orig, 1, MPI_INT, best_rank, MPI_COMM_WORLD);
          MPI_Bcast(&best_seed_extend.query_end_orig,   1, MPI_INT, best_rank, MPI_COMM_WORLD);
          MPI_Bcast(&best_seed_extend.target_start_orig, 1, MPI_INT, best_rank, MPI_COMM_WORLD);
          MPI_Bcast(&best_seed_extend.target_end_orig,   1, MPI_INT, best_rank, MPI_COMM_WORLD);

          if (world_rank != best_rank) { best_seed_extend.aligned_seq1.resize(len1); best_seed_extend.aligned_seq2.resize(len2); }
          if (len1 > 0) MPI_Bcast(&best_seed_extend.aligned_seq1[0], len1, MPI_CHAR, best_rank, MPI_COMM_WORLD);
          if (len2 > 0) MPI_Bcast(&best_seed_extend.aligned_seq2[0], len2, MPI_CHAR, best_rank, MPI_COMM_WORLD);
        }
    }
  }

  AlignmentResult final_aln; final_aln.score = 0;
  std::vector<std::pair<int, int>> local_path;

  if (use_fmindex && best_seed_extend.score > 0) {
    final_aln = best_seed_extend;
  } else {
    if (world_rank == 0 && verbose) std::cout << "Local alignment: FM-index anchoring unavailable/failed. Fallback to MPI full DP.\n";
    int start_row = 0, local_rows = 0; get_row_partition(m, world_size, world_rank, start_row, local_rows);
    std::vector<std::vector<int>> s_block(local_rows + 1, std::vector<int>(n + 1, 0)), e_block(local_rows + 1, std::vector<int>(n + 1, 0)), f_block(local_rows + 1, std::vector<int>(n + 1, 0));

    if (world_rank > 0) {
      MPI_Recv(s_block[0].data(), n + 1, MPI_INT, world_rank - 1, 200, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      MPI_Recv(e_block[0].data(), n + 1, MPI_INT, world_rank - 1, 201, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      MPI_Recv(f_block[0].data(), n + 1, MPI_INT, world_rank - 1, 202, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }

    int local_best_score = 0, local_best_i = 0, local_best_j = 0;
    std::vector<int> curr_s, curr_e, curr_f; std::vector<char> curr_tr;
    for (int ii = 1; ii <= local_rows; ++ii) {
      int global_i = start_row + ii;
      computeAffineDPRow(global_i, x, y, s_block[ii - 1], e_block[ii - 1], f_block[ii - 1], curr_s, curr_e, curr_f, curr_tr, sfn_val, false);
      s_block[ii] = curr_s; e_block[ii] = curr_e; f_block[ii] = curr_f;
      for (int j = 1; j <= n; ++j) { if (s_block[ii][j] > local_best_score) { local_best_score = s_block[ii][j]; local_best_i = global_i; local_best_j = j; } }
      if (verbose && world_rank == 0 && (global_i % 100 == 0 || global_i == m)) showProgressBar(global_i, m);
    }

    if (world_rank + 1 < world_size) {
      MPI_Send(s_block[local_rows].data(), n + 1, MPI_INT, world_rank + 1, 200, MPI_COMM_WORLD);
      MPI_Send(e_block[local_rows].data(), n + 1, MPI_INT, world_rank + 1, 201, MPI_COMM_WORLD);
      MPI_Send(f_block[local_rows].data(), n + 1, MPI_INT, world_rank + 1, 202, MPI_COMM_WORLD);
    }

    int mine[4] = {local_best_score, world_rank, local_best_i, local_best_j};
    std::vector<int> gathered; if (world_rank == 0) gathered.resize(4 * world_size, 0);
    MPI_Gather(mine, 4, MPI_INT, world_rank == 0 ? gathered.data() : nullptr, 4, MPI_INT, 0, MPI_COMM_WORLD);

    int best_score = 0, best_i = 0, best_j = 0;
    if (world_rank == 0) {
      for (int r = 0; r < world_size; ++r) {
        if (gathered[4 * r + 0] > best_score) { best_score = gathered[4 * r + 0]; best_i = gathered[4 * r + 2]; best_j = gathered[4 * r + 3]; }
      }
    }

    MPI_Bcast(&best_score, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&best_i, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&best_j, 1, MPI_INT, 0, MPI_COMM_WORLD);

    std::vector<std::vector<int>> fullS = gather_score_matrix_rows(s_block, m, n, false, MPI_COMM_WORLD);
    if (world_rank == 0) {
      if (binary) writeDPMatrix(fullS, odir + "/local_dp_matrix.bin");
      else if (txt) writeRawDPMatrix(fullS, odir + "/local_dp_matrix.txt");
      final_aln = traceback_local_from_full(x, y, sfn_val, best_i, best_j, local_path);
      final_aln.score = best_score;
    }
  }

  auto t_end = std::chrono::high_resolution_clock::now();
  final_aln.time_ms = std::chrono::duration_cast<std::chrono::milliseconds>(t_end - t_start).count();

  if (world_rank == 0) {
    if (use_fmindex) {
        if (txt || binary) std::cout << "\nNotice: DP Matrix skipped during anchored Local Alignment.\n";
        int tr_cx = final_aln.query_end_orig + 1, tr_cy = final_aln.target_end_orig + 1;
        local_path.emplace_back(tr_cy, tr_cx);
        for (int i = final_aln.aligned_seq1.size() - 1; i >= 0; --i) {
            if (final_aln.aligned_seq1[i] != '-') tr_cx--;
            if (final_aln.aligned_seq2[i] != '-') tr_cy--;
            local_path.emplace_back(tr_cy, tr_cx);
        }
        std::reverse(local_path.begin(), local_path.end());
    }

    size_t total = final_aln.aligned_seq1.size(), gaps = 0, matches = 0;
    for (size_t i = 0; i < total; ++i) {
      if (final_aln.aligned_seq1[i] == '-' || final_aln.aligned_seq2[i] == '-') ++gaps;
      else if (final_aln.aligned_seq1[i] == final_aln.aligned_seq2[i]) ++matches;
    }
    double identity = total > 0 ? static_cast<double>(matches) / total : 0.0;
    double coverage = total > 0 ? static_cast<double>(total - gaps) / total : 0.0;

    std::string acc1 = getAccession(h1, mval), acc2 = getAccession(h2, mval);
    std::string gene1 = getGeneSymbol(h1, mval), gene2 = getGeneSymbol(h2, mval);

    save_path_file(odir + "/local_path.txt", local_path);

    if (verbose) {
      std::cout << "\n\nLocal Alignment Score: " << final_aln.score << "\n";
      std::cout << "Matches: " << matches << " | Gaps: " << gaps << " | Aligned Length: " << total << "\n";
      std::cout << "Time: " << final_aln.time_ms << " ms\n\n";
      if (final_aln.score > 0) printColoredAlignment(final_aln.aligned_seq1, final_aln.aligned_seq2);
    }

    std::ofstream fa(odir + "/local_alignment.fasta");
    if (fa) savePlainAlignment(acc1 + "_local", acc2 + "_local", final_aln.aligned_seq1, final_aln.aligned_seq2, fa);

    std::ofstream js(odir + "/local_stats.json");
    if (js) {
      js << std::fixed << std::setprecision(6) << "{\n"
         << "  \"method\": \"local\",\n"
         << "  \"score\": " << final_aln.score << ",\n"
         << "  \"matches\": " << matches << ",\n"
         << "  \"gaps\": " << gaps << ",\n"
         << "  \"aligned_length\": " << total << ",\n"
         << "  \"identity\": " << identity << ",\n"
         << "  \"coverage_aligned\": " << coverage << ",\n"
         << "  \"time_ms\": " << final_aln.time_ms << ",\n"
         << R"(  "query": ")" << acc1 << "\",\n"
         << R"(  "target": ")" << acc2 << "\",\n"
         << R"(  "queryid": ")" << gene1 << "\",\n"
         << R"(  "targetid": ")" << gene2 << "\",\n"
         << "  \"query_length\": " << m << ",\n"
         << "  \"target_length\": " << n << "\n"
         << "}\n";
    }
  }
  MPI_Barrier(MPI_COMM_WORLD);
}

void lcs(const std::string& x_o, const std::string& y_o,
         const std::string& h1_lcs, const std::string& h2_lcs,
         const std::string& odir_lcs, ScoreMode mode_lcs,
         const FMIndex* tfm_idx_lcs) {
  const int m = static_cast<int>(x_o.size()), n = static_cast<int>(y_o.size());
  int world_rank = rank_val, world_size = 1;
  MPI_Comm_size(MPI_COMM_WORLD, &world_size);
  auto ts_lcs = std::chrono::high_resolution_clock::now();

  bool use_anchoring = false; ChainedSeed best_chain; LcsSegmentResult final_res;
  std::vector<std::pair<int, int>> lcs_path_coords;

  if (world_rank == 0 && tfm_idx_lcs != nullptr) {
    int k = std::min(10, std::min(m / 15, n / 15));
    if (std::min(m, n) < k) k = std::min(m, n);
    if (k > 0) {
      std::vector<Seed> raw_seeds = generate_raw_seeds(x_o, *tfm_idx_lcs, k, 0, 1);
      if (!raw_seeds.empty()) { best_chain = find_best_seed_chain(raw_seeds, 1); use_anchoring = !best_chain.seeds.empty(); }
    }
  }

  int anchor_count = (world_rank == 0 && use_anchoring) ? static_cast<int>(best_chain.seeds.size()) : 0;
  MPI_Bcast(&anchor_count, 1, MPI_INT, 0, MPI_COMM_WORLD);
  use_anchoring = anchor_count > 0;
  if (use_anchoring && world_rank != 0) best_chain.seeds.resize(anchor_count);
  if (use_anchoring) MPI_Bcast(best_chain.seeds.data(), anchor_count * static_cast<int>(sizeof(Seed)), MPI_BYTE, 0, MPI_COMM_WORLD);

  if (use_anchoring) {
    std::vector<std::pair<std::string, std::string>> segments;
    int cx = 0, cy = 0;
    for (const auto& anc : best_chain.seeds) {
      segments.push_back({ x_o.substr(cx, anc.query_pos - cx), y_o.substr(cy, anc.target_pos - cy) });
      cx = anc.query_pos + anc.len; cy = anc.target_pos + anc.len;
    }
    segments.push_back({x_o.substr(cx), y_o.substr(cy)});
    std::vector<LcsSegmentResult> seg_results(segments.size());

    for (size_t i = 0; i < segments.size(); ++i) {
      if (static_cast<int>(i % world_size) == world_rank) seg_results[i] = compute_lcs_for_segment(segments[i].first, segments[i].second);
    }

    if (world_rank == 0) {
      for (size_t i = 0; i < segments.size(); ++i) {
        int owner = static_cast<int>(i % world_size);
        if (owner != 0) {
          MPI_Status st; MPI_Recv(&seg_results[i].lcs_length, 1, MPI_INT, owner, static_cast<int>(i) * 10 + 0, MPI_COMM_WORLD, &st);
          int len1 = 0, len2 = 0, len3 = 0;
          MPI_Recv(&len1, 1, MPI_INT, owner, static_cast<int>(i) * 10 + 1, MPI_COMM_WORLD, &st);
          seg_results[i].lcs_string.resize(len1);
          if (len1 > 0) MPI_Recv(&seg_results[i].lcs_string[0], len1, MPI_CHAR, owner, static_cast<int>(i) * 10 + 2, MPI_COMM_WORLD, &st);
          MPI_Recv(&len2, 1, MPI_INT, owner, static_cast<int>(i) * 10 + 3, MPI_COMM_WORLD, &st);
          seg_results[i].gapped_seq1.resize(len2);
          if (len2 > 0) MPI_Recv(&seg_results[i].gapped_seq1[0], len2, MPI_CHAR, owner, static_cast<int>(i) * 10 + 4, MPI_COMM_WORLD, &st);
          MPI_Recv(&len3, 1, MPI_INT, owner, static_cast<int>(i) * 10 + 5, MPI_COMM_WORLD, &st);
          seg_results[i].gapped_seq2.resize(len3);
          if (len3 > 0) MPI_Recv(&seg_results[i].gapped_seq2[0], len3, MPI_CHAR, owner, static_cast<int>(i) * 10 + 6, MPI_COMM_WORLD, &st);
        }
      }
      for (size_t i = 0; i < best_chain.seeds.size(); ++i) {
        final_res.lcs_string += seg_results[i].lcs_string; final_res.lcs_length += seg_results[i].lcs_length;
        final_res.gapped_seq1 += seg_results[i].gapped_seq1; final_res.gapped_seq2 += seg_results[i].gapped_seq2;
        const auto& anc = best_chain.seeds[i];
        std::string exact = x_o.substr(anc.query_pos, anc.len);
        final_res.lcs_string += exact; final_res.lcs_length += anc.len;
        final_res.gapped_seq1 += exact; final_res.gapped_seq2 += exact;
      }
      final_res.lcs_string += seg_results.back().lcs_string; final_res.lcs_length += seg_results.back().lcs_length;
      final_res.gapped_seq1 += seg_results.back().gapped_seq1; final_res.gapped_seq2 += seg_results.back().gapped_seq2;
    } else {
      for (size_t i = 0; i < segments.size(); ++i) {
        if (static_cast<int>(i % world_size) == world_rank) {
          MPI_Send(&seg_results[i].lcs_length, 1, MPI_INT, 0, static_cast<int>(i) * 10 + 0, MPI_COMM_WORLD);
          int len1 = seg_results[i].lcs_string.size(), len2 = seg_results[i].gapped_seq1.size(), len3 = seg_results[i].gapped_seq2.size();
          MPI_Send(&len1, 1, MPI_INT, 0, static_cast<int>(i) * 10 + 1, MPI_COMM_WORLD);
          if (len1 > 0) MPI_Send(seg_results[i].lcs_string.data(), len1, MPI_CHAR, 0, static_cast<int>(i) * 10 + 2, MPI_COMM_WORLD);
          MPI_Send(&len2, 1, MPI_INT, 0, static_cast<int>(i) * 10 + 3, MPI_COMM_WORLD);
          if (len2 > 0) MPI_Send(seg_results[i].gapped_seq1.data(), len2, MPI_CHAR, 0, static_cast<int>(i) * 10 + 4, MPI_COMM_WORLD);
          MPI_Send(&len3, 1, MPI_INT, 0, static_cast<int>(i) * 10 + 5, MPI_COMM_WORLD);
          if (len3 > 0) MPI_Send(seg_results[i].gapped_seq2.data(), len3, MPI_CHAR, 0, static_cast<int>(i) * 10 + 6, MPI_COMM_WORLD);
        }
      }
    }
  } else {
    int start_row = 0, local_rows = 0; get_row_partition(m, world_size, world_rank, start_row, local_rows);
    std::vector<std::vector<int>> Lblock(local_rows + 1, std::vector<int>(n + 1, 0));
    if (world_rank > 0) MPI_Recv(Lblock[0].data(), n + 1, MPI_INT, world_rank - 1, 300, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

    for (int ii = 1; ii <= local_rows; ++ii) {
      int global_i = start_row + ii; Lblock[ii][0] = 0;
      for (int j = 1; j <= n; ++j) {
        if (x_o[global_i - 1] == y_o[j - 1]) Lblock[ii][j] = Lblock[ii - 1][j - 1] + 1;
        else Lblock[ii][j] = std::max(Lblock[ii - 1][j], Lblock[ii][j - 1]);
      }
    }

    if (world_rank + 1 < world_size) MPI_Send(Lblock[local_rows].data(), n + 1, MPI_INT, world_rank + 1, 300, MPI_COMM_WORLD);
    std::vector<std::vector<int>> fullL = gather_score_matrix_rows(Lblock, m, n, false, MPI_COMM_WORLD);

    if (world_rank == 0) {
      std::vector<std::vector<char>> B;
      final_res = traceback_lcs_from_full(x_o, y_o, fullL, B, lcs_path_coords);
      if (binary) writeDPMatrix(fullL, odir_lcs + "/lcs_dp_lengths.bin"); else if (txt) writeRawDPMatrix(fullL, odir_lcs + "/lcs_dp_lengths.txt");
      if (binary) writeCharMatrix(B, odir_lcs + "/lcs_traceback_pointers.bin"); else if (txt) writeRawCharMatrix(B, odir_lcs + "/lcs_traceback_pointers.txt");
    }
  }

  auto te_lcs = std::chrono::high_resolution_clock::now();
  long long tms_lcs = std::chrono::duration_cast<std::chrono::milliseconds>(te_lcs - ts_lcs).count();

  if (world_rank == 0) {
    if (use_anchoring) {
        if (txt || binary) std::cout << "\nNotice: DP Matrix skipped during anchored LCS.\n";
        int tr_cx = m, tr_cy = n;
        lcs_path_coords.emplace_back(tr_cy, tr_cx);
        for (int i = final_res.gapped_seq1.size() - 1; i >= 0; --i) {
            if (final_res.gapped_seq1[i] != '-') tr_cx--;
            if (final_res.gapped_seq2[i] != '-') tr_cy--;
            lcs_path_coords.emplace_back(tr_cy, tr_cx);
        }
        std::reverse(lcs_path_coords.begin(), lcs_path_coords.end());
    }

    std::string acc1 = getAccession(h1_lcs, mode_lcs), acc2 = getAccession(h2_lcs, mode_lcs);
    std::ofstream out_lcs(odir_lcs + "/lcs.fasta");
    if (out_lcs) saveLCS(acc1 + "_" + acc2, final_res.lcs_string, out_lcs);

    std::ofstream out_aln(odir_lcs + "/lcs_alignment.fasta");
    if (out_aln) savePlainAlignment(acc1 + "_LCS_aligned", acc2 + "_LCS_aligned", final_res.gapped_seq1, final_res.gapped_seq2, out_aln);

    save_path_file(odir_lcs + "/lcs_path.txt", lcs_path_coords);

    if (verbose) {
      std::cout << "\n--- LCS Final Length: " << final_res.lcs_length << "\n";
      std::cout << "Time: " << tms_lcs << " ms\n";
    }
  }
  MPI_Barrier(MPI_COMM_WORLD);
}

// -------- Main Function --------
int main(int argc, char** argv) {
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank_val);

  std::unique_ptr<FMIndex> fm_index_target_obj_ptr = nullptr;
  try {
    std::string query_file, target_file, output_dir = ".";
    std::string fm_idx_path_for_target = "";
    int         align_choice           = -1;
    ScoreMode   current_mode           = MODE_DNA;

    for (int i = 1; i < argc; ++i) {
      std::string current_arg = argv[i];
      if (current_arg == "--query" && i + 1 < argc) query_file = argv[++i];
      else if (current_arg == "--target" && i + 1 < argc) target_file = argv[++i];
      else if (current_arg == "--choice" && i + 1 < argc) align_choice = std::stoi(argv[++i]);
      else if (current_arg == "--mode" && i + 1 < argc) {
        std::string mode_str = argv[++i];
        if (mode_str == "dna") current_mode = MODE_DNA;
        else if (mode_str == "protein") current_mode = MODE_PROTEIN;
        else {
          if (rank_val == 0) std::cerr << "Unknown mode: " << mode_str << "\n";
          MPI_Abort(MPI_COMM_WORLD, 1);
          return 1;
        }
      } else if (current_arg == "--outdir" && i + 1 < argc) output_dir = argv[++i];
      else if (current_arg == "--fmindex" && i + 1 < argc) fm_idx_path_for_target = argv[++i];
      else if (current_arg == "--verbose") verbose = true;
      else if (current_arg == "--binary") binary = true;
      else if (current_arg == "--txt") txt = true;
      else if (current_arg == "--gap_open" && i + 1 < argc) GAP_OPEN = std::stod(argv[++i]);
      else if (current_arg == "--gap_extend" && i + 1 < argc) GAP_EXTEND = std::stod(argv[++i]);
      else if (current_arg == "--help") {
        if (rank_val == 0) std::cout << "Usage: ... (Full help message)\n";
        MPI_Finalize();
        return 0;
      } else {
        if (rank_val == 0) std::cerr << "Unknown option: " << current_arg << "\n";
        MPI_Abort(MPI_COMM_WORLD, 1);
        return 1;
      }
    }

    ScoreFn current_score_fn = (current_mode == MODE_DNA) ? &edna_score : &blosum62_score;
    if (query_file.empty() || target_file.empty() || align_choice == -1) {
      if (rank_val == 0) std::cerr << "Missing required arguments...\n";
      MPI_Abort(MPI_COMM_WORLD, 1);
      return 1;
    }
    if (rank_val == 0) {
      try {
        std::filesystem::create_directories(output_dir);
      } catch (const std::exception& e) {
        std::cerr << "Error creating output dir " << output_dir << ": " << e.what() << std::endl;
        MPI_Abort(MPI_COMM_WORLD, 1);
        return 1;
      }
    }

    std::string seq1_data, seq2_data, h1_data, h2_data;
    if (rank_val == 0) {
      try {
        processFasta(query_file, h1_data, seq1_data);
        processFasta(target_file, h2_data, seq2_data);
      } catch (const std::runtime_error& e) {
        std::cerr << "FASTA error: " << e.what() << std::endl;
        MPI_Abort(MPI_COMM_WORLD, 1);
        return 1;
      }
      if (!fm_idx_path_for_target.empty()) {
        std::ifstream fm_ifs(fm_idx_path_for_target, std::ios::binary);
        if (fm_ifs) {
          fm_index_target_obj_ptr = std::make_unique<FMIndex>();
          if (!fm_index_target_obj_ptr->load(fm_ifs)) {
            std::cerr << "Rank 0: Error! Failed to load FM-Index from " << fm_idx_path_for_target << "\n";
            fm_index_target_obj_ptr.reset();
          } else if (verbose) {
            std::cout << "Rank 0: Loaded FM-Index for target from " << fm_idx_path_for_target << std::endl;
          }
        }
      }
    }

    int s1_len_val = (rank_val == 0) ? seq1_data.length() : 0;
    int s2_len_val = (rank_val == 0) ? seq2_data.length() : 0;
    MPI_Bcast(&s1_len_val, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&s2_len_val, 1, MPI_INT, 0, MPI_COMM_WORLD);
    if (rank_val != 0) {
      seq1_data.resize(s1_len_val);
      seq2_data.resize(s2_len_val);
    }
    if (s1_len_val > 0) MPI_Bcast(&seq1_data[0], s1_len_val, MPI_CHAR, 0, MPI_COMM_WORLD);
    if (s2_len_val > 0) MPI_Bcast(&seq2_data[0], s2_len_val, MPI_CHAR, 0, MPI_COMM_WORLD);

    int h1_len = (rank_val == 0) ? h1_data.length() : 0;
    int h2_len = (rank_val == 0) ? h2_data.length() : 0;
    MPI_Bcast(&h1_len, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&h2_len, 1, MPI_INT, 0, MPI_COMM_WORLD);
    if (rank_val != 0) { h1_data.resize(h1_len); h2_data.resize(h2_len); }
    if (h1_len > 0) MPI_Bcast(&h1_data[0], h1_len, MPI_CHAR, 0, MPI_COMM_WORLD);
    if (h2_len > 0) MPI_Bcast(&h2_data[0], h2_len, MPI_CHAR, 0, MPI_COMM_WORLD);

    const FMIndex* fm_target_raw_ptr_val = (rank_val == 0) ? fm_index_target_obj_ptr.get() : nullptr;

    if (align_choice == 1) {
      globalalign(seq1_data, seq2_data, h1_data, h2_data, output_dir, current_mode, current_score_fn, fm_target_raw_ptr_val);
    } else if (align_choice == 2) {
      localalign(seq1_data, seq2_data, h1_data, h2_data, output_dir, current_mode, current_score_fn, fm_target_raw_ptr_val);
    } else if (align_choice == 3) {
      lcs(seq1_data, seq2_data, h1_data, h2_data, output_dir, current_mode, fm_target_raw_ptr_val);
    } else if (align_choice == 4) {
      globalalign(seq1_data, seq2_data, h1_data, h2_data, output_dir, current_mode, current_score_fn, fm_target_raw_ptr_val);
      MPI_Barrier(MPI_COMM_WORLD);
      localalign(seq1_data, seq2_data, h1_data, h2_data, output_dir, current_mode, current_score_fn, fm_target_raw_ptr_val);
      MPI_Barrier(MPI_COMM_WORLD);
      lcs(seq1_data, seq2_data, h1_data, h2_data, output_dir, current_mode, fm_target_raw_ptr_val);
    } else if (rank_val == 0) {
      std::cerr << "Invalid choice. Use --choice 1/2/3/4.\n";
    }

  } catch (const std::exception& e) {
    if (rank_val == 0) std::cerr << "Runtime Exception: " << e.what() << "\n";
    MPI_Abort(MPI_COMM_WORLD, 1);
  } catch (...) {
    if (rank_val == 0) std::cerr << "Unknown exception.\n";
    MPI_Abort(MPI_COMM_WORLD, 1);
  }

  MPI_Finalize();
  return 0;
}