/*
* @file main.cpp
* @brief Global and local sequence alignment using MPI and OpenMP.
* @author Abhinav Mishra
* @date 2023-10-01
*
* This program performs global and local sequence alignment using the
* Needleman-Wunsch and Smith-Waterman algorithms, respectively. It uses
* MPI for parallel processing across multiple nodes and OpenMP for
* parallelization within each node.
*/
#include <zlib.h>
#include <array>
#include <cstring>
#include <iostream>
#include <chrono>
#include <iomanip>
#include <fstream>
#include <vector>
#include <algorithm>
#include <stdexcept>
#include <omp.h>
#include <immintrin.h>
#include <mpi.h>
#include <sstream>
#include "EDNAFULL.h"
#include "EBLOSUM62.h"
#include <unordered_map>
#include <filesystem>
#include <climits>

bool verbose = false;
bool binary = false;
bool txt = false;
using namespace std;

using score_t = int;                           // the cell‐type in your DP
static std::vector<score_t> dp_prev, dp_curr;  // two rolling rows

static void allocThreadBuffers(int n) {
    // allocate (n+1) entries in each row, zero-initialized
    dp_prev.assign(n+1, 0);
    dp_curr.assign(n+1, 0);
}

enum ScoreMode { MODE_DNA, MODE_PROTEIN };
using ScoreFn = int(*)(char,char);

/// Gap penalty
double GAP_OPEN   = -5.0;  // penalty to open a gap
double GAP_EXTEND = -1.0;  // penalty to extend an existing gap in global alignment
/// Number of columns per line when printing alignments
static const int LINE_WIDTH = 120;

/// ANSI escape code: reset color
#define RESET "\033[0m"
/// ANSI escape code: green
#define GREEN "\033[32m"
/// ANSI escape code: red
#define RED   "\033[31m"
/// ANSI escape code: cyan
#define CYAN  "\033[36m"

// DNA / EDNAFULL lookup
static const std::array<uint8_t,256> char2idx = [](){
    std::array<uint8_t,256> m{};
    m.fill(255);  // mark “invalid”
    m['A'] =  0;  m['C'] =  1;  m['G'] =  2;  m['T'] =  3;  m['U'] =  3; // T=U
    m['R'] =  4;  m['Y'] =  5;  m['S'] =  6;  m['W'] =  7;
    m['K'] =  8;  m['M'] =  9;  m['B'] = 10;  m['D'] = 11;
    m['H'] = 12;  m['V'] = 13;  m['N'] = 14;  m['X'] = 14;
    return m;
}();

// Protein / BLOSUM62 lookup
static const std::array<uint8_t,256> prot_idx = [](){
    std::array<uint8_t,256> m{};
    m.fill(255);
    m['A']=0;   m['R']=1;   m['N']=2;   m['D']=3;   m['C']=4;
    m['Q']=5;   m['E']=6;   m['G']=7;   m['H']=8;   m['I']=9;
    m['L']=10;  m['K']=11;  m['M']=12;  m['F']=13;  m['P']=14;
    m['S']=15;  m['T']=16;  m['W']=17;  m['Y']=18;  m['V']=19;
    m['B']=20;  m['Z']=21;  m['X']=22;  m['*']=23;
    return m;
}();

/**
 * @brief Lookup the score between two characters based on the selected mode.
 * @param x First character (base or amino acid).
 * @param y Second character (base or amino acid).
 * @return Score based on the selected scoring matrix.
**/
inline int edna_score(char x, char y) {
    uint8_t ix = char2idx[uint8_t(x)];
    uint8_t iy = char2idx[uint8_t(y)];
    if (ix == 255 || iy == 255)
        throw std::runtime_error(std::string("Invalid DNA code: ") + x + "," + y);
    return EDNAFULL_matrix[ix][iy];
}

/**
 * @brief Lookup the score between two characters based on the selected mode.
 * @param x First character (base or amino acid).
 * @param y Second character (base or amino acid).
 * @return Score based on the selected scoring matrix.
**/
inline int blosum62_score(char x, char y) {
    uint8_t ix = prot_idx[uint8_t(x)];
    uint8_t iy = prot_idx[uint8_t(y)];
    if (ix == 255 || iy == 255)
        throw std::runtime_error(std::string("Invalid protein code: ") + x + "," + y);
    return EBLOSUM62_matrix[ix][iy];
}

/**
 * @brief Lookup the score between two characters based on the selected mode.
 * @param x First character (base or amino acid).
 * @param y Second character (base or amino acid).
 * @param mode Scoring mode (MODE_DNA or MODE_PROTEIN).
 * @return Score based on the selected scoring matrix.
**/
inline int score(char x, char y, ScoreMode mode) {
    if (mode == MODE_DNA) return edna_score(x, y);
    else return blosum62_score(x, y);
}


/** * @brief Show a progress bar in the console.
 *
 * @param progress Current progress (0 to total).
 * @param total    Total number of steps.
 */
void showProgressBar(int progress, int total) {
    using namespace std::chrono;

    static auto start_time = steady_clock::now();
    auto now     = steady_clock::now();
    auto elapsed = duration_cast<seconds>(now - start_time).count();

    // compute ETA in seconds
    long eta = 0;
    if (progress > 0 && progress < total) {
        eta = elapsed * (total - progress) / progress;
    }

    // format elapsed and ETA as H:MM:SS
    auto format_hms = [](long secs) {
        long h = secs / 3600;
        long m = (secs % 3600) / 60;
        long s = secs % 60;
        std::ostringstream os;
        if (h) os << h << ":";
        os << std::setw(2) << std::setfill('0') << m << ":"
           << std::setw(2) << std::setfill('0') << s;
        return os.str();
    };

    constexpr int barWidth = 100;
    float ratio = float(progress) / total;
    int pos = int(barWidth * ratio);

    // render bar
    std::cout << "\r[";
    for (int i = 0; i < barWidth; ++i) {
        if      (i < pos)  std::cout << "=";
        else if (i == pos) std::cout << ">";
        else               std::cout << " ";
    }
    std::cout << "] "
              << std::setw(3) << int(ratio * 100) << "% "
              << progress << "/" << total
              << " Elapsed: " << format_hms(elapsed)
              << " ETA: "     << format_hms(eta)
              << std::flush;
}

/**
 * @brief Extract the accession number from a FASTA header line.
 *
 * The accession number is the first word in the header line, which
 * is expected to start with a `>`.
 *
 * @param header The header line from a FASTA file.
 * @return The accession number (first word).
 */
std::string getAccession(const std::string& header, ScoreMode mode) {
    if (mode == MODE_DNA) {
        std::istringstream iss(header);
        std::string accession;
        iss >> accession;
        return accession;
    } else if (mode == MODE_PROTEIN) {
        size_t firstPipe = header.find('|');
        size_t secondPipe = header.find('|', firstPipe + 1);
        if (firstPipe != std::string::npos && secondPipe != std::string::npos) {
            return header.substr(firstPipe + 1, secondPipe - firstPipe - 1);
        } else {
            // Fallback to the first word if pipes are not found
            return header;
        }
    }
    // If neither mode matches, return the header as is
    return header;
}

/**
 * @brief Extract the gene symbol from a FASTA header.
 *
 * For DNA headers, looks for the first pair of parentheses “(GENE)” and returns GENE.
 * For protein headers, takes the part after the second ‘|’ up to the first underscore.
 *
 * @param header The FASTA header line (without the leading '>').
 * @param mode   MODE_DNA or MODE_PROTEIN.
 * @return       The gene symbol, or an empty string on failure.
 */
std::string getGeneSymbol(const std::string& header, ScoreMode mode) {
    if (mode == MODE_DNA) {

        size_t open  = header.find('(');
        size_t close = (open != std::string::npos) ? header.find(')', open + 1) : std::string::npos;
        if (open != std::string::npos && close != std::string::npos && close > open + 1) {
            return header.substr(open + 1, close - open - 1);
        }
        // fallback: no parentheses found
        return "";
    }
    else if (mode == MODE_PROTEIN) {

        size_t firstPipe  = header.find('|');
        size_t secondPipe = (firstPipe != std::string::npos)
                            ? header.find('|', firstPipe + 1)
                            : std::string::npos;
        if (secondPipe != std::string::npos) {
            size_t underscore = header.find('_', secondPipe + 1);
            if (underscore != std::string::npos && underscore > secondPipe + 1) {
                return header.substr(secondPipe + 1, underscore - secondPipe - 1);
            }
        }
        // fallback: try GN= field
        size_t gnPos = header.find("GN=");
        if (gnPos != std::string::npos) {
            size_t start = gnPos + 3;
            size_t end   = header.find_first_of(" ;", start);
            if (end == std::string::npos) end = header.size();
            if (end > start) return header.substr(start, end - start);
        }
        return "";
    }
    // unknown mode
    return "";
}


/**
 * @brief Read the first record from a FASTA file.
 *
 * Extracts the very first header (minus the leading `>`), and concatenates
 * all subsequent lines into a single sequence string. If a target chromosome
 * is specified, only that record is processed.
 *
 * @param filename  Path to the FASTA file.
 * @param[out] header   On return, the header line without leading `>`.
 * @param[out] sequence On return, the full sequence (no newlines).
 * @param targetChr  Optional chromosome ID to filter by (empty means no filter).
 * @throws std::runtime_error if the file cannot be opened.
 */
void processFasta(const std::string &filename,
                  std::string &header,
                  std::string &sequence,
                  const std::string &targetChr = "")
{
    std::ifstream file(filename);
    if (!file) throw std::runtime_error("Error: Unable to open " + filename);
    std::string line;
    header.clear();
    sequence.clear();
    bool inRec = targetChr.empty();   // if no filter, grab first record
    while (std::getline(file, line)) {
        if (line.empty()) continue;
        if (line[0] == '>') {
            // strip leading '>' and up to first space
            std::string h = line.substr(1);
            std::string id = h.substr(0, h.find_first_of(" \t"));
            if (!targetChr.empty()) {
                if (id == targetChr) {
                    inRec = true;
                    header = h;
                } else if (inRec) {
                    // we just left our target record → done
                    break;
                } else {
                    inRec = false;
                }
            } else if (!inRec) {
                // first header ever
                header = h;
                inRec = true;
            }
            continue;
        }
        if (inRec) sequence += line;
    }
}


/* * @brief Read the first record from a gzipped FASTA file.
 *
 * Extracts the very first header (minus the leading `>`), and concatenates
 * all subsequent lines into a single sequence string. If a target chromosome
 * is specified, only that record is processed.
 *
 * @param filename  Path to the gzipped FASTA file.
 * @param[out] header   On return, the header line without leading `>`.
 * @param[out] sequence On return, the full sequence (no newlines).
 * @param targetChr  Optional chromosome ID to filter by (empty means no filter).
 * @throws std::runtime_error if the file cannot be opened.
 */
void processFastaGz(const std::string &filename,
                    std::string &header,
                    std::string &sequence,
                    const std::string &targetChr = "")
{
    gzFile file = gzopen(filename.c_str(), "rb");
    if (!file) throw std::runtime_error("Error opening " + filename);
    char buf[8192];
    header.clear();
    sequence.clear();
    bool inRec = targetChr.empty();
    while (gzgets(file, buf, sizeof(buf))) {
        std::string line(buf);
        if (!line.empty() && line.back()=='\n') line.pop_back();
        if (line.empty()) continue;
        if (line[0]=='>') {
            std::string h = line.substr(1);
            std::string id = h.substr(0, h.find_first_of(" \t"));
            if (!targetChr.empty()) {
                if (id==targetChr) {
                    inRec = true; header = h;
                } else if (inRec) {
                    break;
                } else {
                    inRec = false;
                }
            } else if (!inRec) {
                header = h; inRec = true;
            }
            continue;
        }
        if (inRec) sequence += line;
    }
    gzclose(file);
}


/**
 * @brief Write two aligned sequences in FASTA format.
 *
 * @param header1    Identifier for sequence 1 (no leading '>').
 * @param header2    Identifier for sequence 2 (no leading '>').
 * @param aligned1   Aligned sequence 1 (with gaps).
 * @param aligned2   Aligned sequence 2 (with gaps).
 * @param os         Output stream (e.g. ofstream).
 */
void savePlainAlignment(const std::string &header1,
                        const std::string &header2,
                        const std::string &aligned1,
                        const std::string &aligned2,
                        std::ostream   &os)
{
    // Write seq1
    os << '>' << header1 << '\n'
       << aligned1  << '\n'
       // Write seq2
       << '>' << header2 << '\n'
       << aligned2  << '\n';
}

/**
 * @brief Write the LCS (Longest Common Subsequence) in FASTA format.
 *
 * @param id       Identifier for the LCS (no leading '>').
 * @param lcs_str  The LCS string.
 * @param os       Output stream (e.g. ofstream).
 */
void saveLCS(const std::string &id,
             const std::string &lcs_str,
             std::ostream   &os)
{
    os << '>' << id << "_LCS_len=" << lcs_str.size() << "\n";
    // wrap at 80 chars per line:
    for (size_t i = 0; i < lcs_str.size(); i += 80)
        os << lcs_str.substr(i, 80) << "\n";
}


/**
 * @brief Print two aligned sequences side-by-side with colors.
 *
 * Matches are shown in green, gaps in red, mismatches in cyan.  Long
 * alignments are chunked into blocks of LINE_WIDTH for readability.
 *
 * @param seq1 First aligned sequence (may contain ‘-’).
 * @param seq2 Second aligned sequence (may contain ‘-’).
 * @param os   Output stream to use (defaults to std::cout).
 */
void printColoredAlignment(const string &seq1, const string &seq2, ostream &os = cout) {
    size_t length1 = seq1.size();
    size_t length2 = seq2.size();
    size_t maxLength = max(length1, length2);

    string aligned1 = seq1;
    string aligned2 = seq2;
    if (length1 < maxLength) aligned1.append(maxLength - length1, '-');
    if (length2 < maxLength) aligned2.append(maxLength - length2, '-');

    for (size_t i = 0; i < maxLength; i += LINE_WIDTH) {
        size_t end = min(i + LINE_WIDTH, maxLength);
        os << "\nPosition: " << i + 1 << " - " << end << "\n";

        for (size_t j = i; j < end; j++) {
            if (aligned1[j] == aligned2[j]) os << GREEN << aligned1[j] << RESET;
            else if (aligned1[j] == '-' || aligned2[j] == '-') os << RED << aligned1[j] << RESET;
            else os << CYAN << aligned1[j] << RESET;
        }
        os << "\n";

        for (size_t j = i; j < end; j++) {
            if (aligned1[j] == aligned2[j]) os << GREEN << aligned2[j] << RESET;
            else if (aligned1[j] == '-' || aligned2[j] == '-') os << RED << aligned2[j] << RESET;
            else os << CYAN << aligned2[j] << RESET;
        }
        os << "\n";
    }
}

/**
 * @brief Write the DP matrix to a file.
 *
 * @param dp       The DP matrix (2D vector).
 * @param filename The output filename.
 */
void writeRawDPMatrix(const std::vector<std::vector<int>>& dp, const std::string& filename) {
    std::ofstream out(filename);
    if (!out) {
        std::cerr << "Error: Cannot write DP matrix to " << filename << "\n";
        return;
    }

    for (const auto& row : dp) {
        for (size_t j = 0; j < row.size(); ++j) {
            out << std::setw(5) << row[j];
            if (j != row.size() - 1) out << " ";
        }
        out << "\n";
    }

    out.close();
}

/**
 * @brief Write the DP matrix to a binary file.
 *
 * @param dp       The DP matrix (2D vector).
 * @param filename The output filename.
 */
void writeDPMatrix(const std::vector<std::vector<int>>& dp, const std::string& filename) {
    std::ofstream out(filename, std::ios::binary);
    if (!out) {
        std::cerr << "Error: Cannot write binary DP matrix to " << filename << "\n";
        return;
    }

    int32_t rows = dp.size();
    int32_t cols = dp[0].size();

    // Write matrix shape as 2 int32 headers
    out.write(reinterpret_cast<const char*>(&rows), sizeof(int32_t));
    out.write(reinterpret_cast<const char*>(&cols), sizeof(int32_t));

    for (const auto& row : dp) {
        out.write(reinterpret_cast<const char*>(row.data()), cols * sizeof(int32_t));
    }

    out.close();
}

/**
 * @brief Write a character matrix to a file.
 *
 * @param mat      The character matrix (2D vector).
 * @param filename The output filename.
 */
void writeRawCharMatrix(const std::vector<std::vector<char>>& mat, const std::string& filename) {
    std::ofstream out(filename);
    if (!out) {
        std::cerr << "Error: Cannot open " << filename << "\n";
        return;
    }

    size_t max_cols = 0;
    for (const auto& row : mat) max_cols = std::max(max_cols, row.size());

    for (const auto& row : mat) {
        for (size_t i = 0; i < max_cols; ++i) {
            char ch = (i < row.size()) ? row[i] : ' ';
            out << ch;
            if (i + 1 < max_cols) out << ' ';
        }
        out << '\n';
    }
}

  /**
 * @brief Write a character matrix to a binary file.
 *
 * @param mat      The character matrix (2D vector).
 * @param filename The output filename.
 */
void writeCharMatrix(const std::vector<std::vector<char>>& mat, const std::string& filename) {
    std::ofstream out(filename, std::ios::binary);
    if (!out) {
        std::cerr << "Error: Cannot open " << filename << " for writing.\n";
        return;
    }

    int32_t rows = mat.size();
    int32_t cols = 0;
    for (const auto& row : mat) cols = std::max<int32_t>(cols, row.size());

    // Write dimensions
    out.write(reinterpret_cast<const char*>(&rows), sizeof(int32_t));
    out.write(reinterpret_cast<const char*>(&cols), sizeof(int32_t));

    // Write character data row by row, pad with ' ' if shorter
    for (const auto& row : mat) {
        for (int32_t i = 0; i < cols; ++i) {
            char ch = (i < static_cast<int32_t>(row.size())) ? row[i] : ' ';
            out.write(&ch, sizeof(char));
        }
    }

    out.close();
}

/* * @brief Initialize the DP row and gap arrays for affine gap scoring.
 *
 * @param n         Length of the second sequence (Y).
 * @param prev_row  Row i-1, match/mismatch scores.
 * @param prev_gapX Row i-1, gap scores in X.
 * @param prev_gapY Row i-1, gap scores in Y.
 * @param isGlobal  Flag to indicate if this is a global alignment.
 */
void initAffineDP(int n,
                  vector<int>& prev_row,
                  vector<int>& prev_gapX,
                  vector<int>& prev_gapY,
                  bool isGlobal)
{
    prev_row .assign(n+1, isGlobal ? INT_MIN/2 : 0);
    prev_gapX.assign(n+1, INT_MIN/2);
    prev_gapY.assign(n+1, isGlobal ? INT_MIN/2 : 0);

    if (isGlobal) {
        // Global: j gaps from the very start
        prev_row[0] = 0;
        for (int j = 1; j <= n; ++j) {
            // opening + (j-1)×extension
            prev_gapY[j] = GAP_OPEN + (j-1)*GAP_EXTEND;
            prev_row [j] = prev_gapY[j];
        }
    }
}


/* * @brief Compute a single row of the DP matrix for affine gap scoring.
 *
 * @param i         Current row index (1-based).
 * @param x         First sequence (X).
 * @param y         Second sequence (Y).
 * @param mode      Scoring mode (MODE_DNA or MODE_PROTEIN).
 * @param prev_row  Row i-1, match/mismatch scores.
 * @param prev_gapX Row i-1, gap scores in X.
 * @param prev_gapY Row i-1, gap scores in Y.
 * @param curr_row  Row i, match/mismatch scores.
 * @param curr_gapX Row i, gap scores in X.
 * @param curr_gapY Row i, gap scores in Y.
 */
void computeAffineDPRow(int i,
                        const string& x,
                        const string& y,
                        vector<int>& prev_row,
                        vector<int>& prev_gapX,
                        vector<int>& prev_gapY,
                        vector<int>& curr_row,
                        vector<int>& curr_gapX,
                        vector<int>& curr_gapY,
                        ScoreFn score_fn)
{
    (void)prev_gapY;
    int n = y.size();

    curr_row.resize(n+1);
    curr_gapX.resize(n+1);
    curr_gapY.resize(n+1);

    std::fill_n(curr_row.data(), n+1, INT_MIN/2);
    std::fill_n(curr_gapX.data(), n+1, INT_MIN/2);
    std::fill_n(curr_gapY.data(), n+1, INT_MIN/2);

    // j=0 column: gap in X down to (i,0)
    curr_gapX[0] = max(prev_row[0]  + GAP_OPEN + GAP_EXTEND,
                       prev_gapX[0] + GAP_EXTEND);
    curr_row[0] = curr_gapX[0];

    for (int j = 1; j <= n; ++j) {
        // open or extend a gap in X (vertical)
        curr_gapX[j] = max(
          prev_row[j]   + GAP_OPEN + GAP_EXTEND,
          prev_gapX[j]  + GAP_EXTEND
        );

        // open or extend a gap in Y (horizontal)
        curr_gapY[j] = max(
          curr_row[j-1] + GAP_OPEN + GAP_EXTEND,
          curr_gapY[j-1] + GAP_EXTEND
        );

        // match/mismatch diagonal
        int mscore = score_fn(x[i-1], y[j-1]);
        int diag   = prev_row[j-1] + mscore;

        // choose best of the three
        curr_row[j] = max({ diag, curr_gapX[j], curr_gapY[j] });
    }
}

/**
 * @brief Perform a global (Needleman–Wunsch) alignment of two sequences.
 *
 * Uses linear-space DP (two-row) with simple traceback to reconstruct the
 * full alignment, then prints and saves it.
 *
 * @param x First sequence.
 * @param y Second sequence.
 */
void globalalign(const string &x, const string &y,
                 const string &header1, const string &header2,
                 const std::string &outdir, ScoreMode mode, ScoreFn score_fn)
{
    int m = x.size(), n = y.size();

    // 1) allocate two rows + trace matrix
    allocThreadBuffers(n);
    // trace[i][j] = 'D','U','L'
    vector<vector<char>> trace(m+1, vector<char>(n+1,'0'));

    // initialize row 0
    for(int j=1; j<=n; ++j) {
        dp_prev[j] = j * (score_t)GAP_OPEN;
        trace[0][j] = 'L';
    }
    dp_prev[0] = 0;

    using Clock = chrono::high_resolution_clock;
    auto t_start = Clock::now();

    // 2) fill by antidiagonals
    #pragma omp parallel
    {
        for(int k=1; k<=m+n; ++k) {
            #pragma omp for schedule(dynamic)
            for(int i = max(1, k-n); i <= min(m, k); ++i) {
                int j = k - i;
                if (j < 1 || j > n) continue;
                // scores
                score_t d = dp_prev[j-1] + (score_t)score_fn(x[i-1], y[j-1]);
                score_t u = dp_prev[j]   + (score_t)GAP_OPEN;
                score_t l = dp_curr[j-1] + (score_t)GAP_OPEN;
                // best
                score_t best = max({d,u,l});
                dp_curr[j] = best;
                // record trace
                if      (best == d) trace[i][j] = 'D';
                else if (best == u) trace[i][j] = 'U';
                else                trace[i][j] = 'L';
            }
            #pragma omp single
            swap(dp_prev, dp_curr);
        }
    }

    // 3) traceback from (m,n)
    vector<pair<int,int>> global_path;
    string alignedX, alignedY;
    int i = m, j = n;
    while(i>0 || j>0) {
        global_path.emplace_back(j,i);
        char c = trace[i][j];
        if      (i>0 && j>0 && c=='D') {
            alignedX.push_back(x[i-1]);
            alignedY.push_back(y[j-1]);
            --i; --j;
        }
        else if (i>0 && c=='U') {
            alignedX.push_back(x[i-1]);
            alignedY.push_back('-');
            --i;
        }
        else {  // 'L'
            alignedX.push_back('-');
            alignedY.push_back(y[j-1]);
            --j;
        }
    }
    reverse(alignedX.begin(),  alignedX.end());
    reverse(alignedY.begin(),  alignedY.end());
    reverse(global_path.begin(), global_path.end());

    // 4) timing & stats
    auto t_end = Clock::now();
    long time_ms = chrono::duration_cast<chrono::milliseconds>(t_end - t_start).count();
    size_t total = alignedX.size(), gaps=0, matches=0;
    for(size_t k=0; k<total; ++k) {
        if (alignedX[k]=='-' || alignedY[k]=='-') ++gaps;
        else if (alignedX[k]==alignedY[k])       ++matches;
    }
    double identity = double(matches)/total;
    double coverage = double(total-gaps)/total;
    score_t finalScore = dp_prev[n];

    // 5) write global_path.txt
    {
      ofstream pf(outdir + "/global_path.txt");
      for (auto &p : global_path)
        pf << p.first << " " << p.second << "\n";
    }

    // 6) verbose print
    if (verbose) {
        cout << "\n\nGlobal Alignment Score: " << finalScore << "\n"
             << "Gap Open: " << GAP_OPEN << "\n"
             << "Gap Extend: " << GAP_EXTEND << "\n"
             << "Matches: " << matches << "\n"
             << "Gaps:    " << gaps   << "\n"
             << "Total:   " << total  << "\n"
             << "Identity: " << identity*100.0 << "%\n"
             << "Coverage: " << coverage*100.0 << "%\n"
             << "Time:    " << time_ms << " ms\n"
             << "Query:   " << header1 << "\n"
             << "Target:  " << header2 << "\n";
        printColoredAlignment(alignedX, alignedY);
    }

    // 7) write FASTA
    {
      ofstream outf(outdir + "/global_alignment.fasta");
      savePlainAlignment(header1, header2, alignedX, alignedY, outf);
    }

    // 8) write JSON stats
    {
      ofstream js(outdir + "/global_stats.json");
      js << fixed << setprecision(6)
         << "{\n"
         << "  \"method\":      \"global\",\n"
         << "  \"gap_open\":    " << GAP_OPEN   << ",\n"
         << "  \"gap_extend\":  " << GAP_EXTEND << ",\n"
         << "  \"score\":       " << finalScore << ",\n"
         << "  \"matches\":     " << matches    << ",\n"
         << "  \"gaps\":        " << gaps       << ",\n"
         << "  \"total\":       " << total      << ",\n"
         << "  \"identity\":    " << identity   << ",\n"
         << "  \"coverage\":    " << coverage   << ",\n"
         << "  \"time_ms\":     " << time_ms    << ",\n"
         << "  \"query\":       \"" << getAccession(header1,mode) << "\",\n"
         << "  \"target\":      \"" << getAccession(header2,mode) << "\"\n"
         << "}\n";
    }
}

//----------------------------------------------------------------------
// C) Striped (AVX2) Smith–Waterman, plus band-window traceback
//----------------------------------------------------------------------

void stripedSmithWaterman(const char*     x,   int m,
                          const char*     y,   int n,
                          int             gap_open,
                          int             gap_ext,
                          ScoreFn         score_fn,
                          std::string&    alnX,
                          std::string&    alnY,
                          std::vector<std::pair<int,int>>& path)

{
    // --- 1) the fast “striped” forward pass to find (end_i,end_j,maxScore) ---
    const int W = 16;                           // 16 lanes of int16
    int segLen = (n + W - 1) / W;
    score_t* Hs = (score_t*)aligned_alloc(64, segLen*W*sizeof(score_t));
    score_t* Es = (score_t*)aligned_alloc(64, segLen*W*sizeof(score_t));
    memset(Hs, 0, segLen*W*sizeof(score_t));
    memset(Es, 0, segLen*W*sizeof(score_t));
    __m256i vGapO = _mm256_set1_epi16((int16_t)gap_open);
    __m256i vGapE = _mm256_set1_epi16((int16_t)gap_ext);

    int end_i=0, end_j=0;
    score_t maxScore = 0;

    for(int i=0;i<m;i++){
        __m256i vF = _mm256_setzero_si256();
        __m256i vHprev = _mm256_setzero_si256();
        for(int s=0;s<segLen;s++){
            // load stripe
            __m256i vH = _mm256_load_si256((__m256i*)&Hs[s*W]);
            __m256i vE = _mm256_load_si256((__m256i*)&Es[s*W]);

            // build match vector
            int16_t mat[W];
            for(int k=0;k<W;k++){
                int j = s*W + k;
                mat[k] = (j<n ? score_fn(x[i], y[j]) : 0);
            }
            __m256i vM = _mm256_load_si256((__m256i*)mat);

            // E = max(E - ext, H - open)
            vE = _mm256_max_epi16(
                    _mm256_subs_epi16(vE, vGapE),
                    _mm256_subs_epi16(vH, vGapO)
                 );
            // F = max(F - ext, Hprev - open)
            vF = _mm256_max_epi16(
                    _mm256_subs_epi16(vF, vGapE),
                    _mm256_subs_epi16(vHprev, vGapO)
                 );
            // H' = H + M
            __m256i vHp = _mm256_adds_epi16(vH, vM);
            // H = max(H', E, F, 0)
            vH = _mm256_max_epi16(vHp, vE);
            vH = _mm256_max_epi16(vH, vF);
            vH = _mm256_max_epi16(vH, _mm256_setzero_si256());

            // store back
            _mm256_store_si256((__m256i*)&Hs[s*W], vH);
            _mm256_store_si256((__m256i*)&Es[s*W], vE);
            vHprev = vH;

            // scan for max
            int16_t buf[W];
            _mm256_store_si256((__m256i*)buf, vH);
            for(int k=0;k<W;k++){
                int j = s*W + k;
                if(j<n && buf[k] > maxScore){
                    maxScore = buf[k];
                    end_i = i;
                    end_j = j;
                }
            }
        }
    }

    // free striped arrays
    free(Hs);
    free(Es);

    // --- 2) banded SW in a small window around (end_i,end_j) for full traceback ---
    const int B = 32;                          // half‐window size
    int i0 = std::max(0, end_i - B), i1 = std::min(m, end_i + B);
    int j0 = std::max(0, end_j - B), j1 = std::min(n, end_j + B);
    int R = i1 - i0 + 1, C = j1 - j0 + 1;

    // allocate band‐matrix and traceback
    std::vector<std::vector<score_t>> H(R, std::vector<score_t>(C, 0));
    std::vector<std::vector<char>>    T(R, std::vector<char>(C, 0));

    // fill local SW
    for(int ii = 1; ii < R; ++ii){
        for(int jj = 1; jj < C; ++jj){
            score_t msc = score_fn(x[i0+ii-1], y[j0+jj-1]);
            score_t v = H[ii-1][jj-1] + msc;
            score_t u = H[ii-1][jj]   + gap_open;
            score_t l = H[ii][jj-1]   + gap_open;
            score_t best = std::max({ (score_t)0, v, u, l });
            H[ii][jj] = best;
            if(best == 0)       T[ii][jj] = '0';
            else if(best == v)  T[ii][jj] = 'D';
            else if(best == u)  T[ii][jj] = 'U';
            else                T[ii][jj] = 'L';
        }
    }

    // traceback from local coords (ei,ej)
    int ii = end_i - i0 + 1, jj = end_j - j0 + 1;
    while(ii > 0 && jj > 0 && H[ii][jj] > 0){
        path.emplace_back(j0 + jj, i0 + ii);
        char c = T[ii][jj];
        if(c == 'D'){
            alnX.push_back(x[i0 + ii - 1]);
            alnY.push_back(y[j0 + jj - 1]);
            --ii; --jj;
        }
        else if(c == 'U'){
            alnX.push_back(x[i0 + ii - 1]);
            alnY.push_back('-');
            --ii;
        }
        else { // 'L'
            alnX.push_back('-');
            alnY.push_back(y[j0 + jj - 1]);
            --jj;
        }
    }
    std::reverse(alnX.begin(), alnX.end());
    std::reverse(alnY.begin(), alnY.end());
    std::reverse(path.begin(), path.end());
}

/**
 * @brief Perform a local (Smith-Waterman) alignment of two sequences.
 *
 * Uses a full DP matrix for the local alignment, and performs traceback
 * to reconstruct the best local alignment.
 *
 * @param x First sequence.
 * @param y Second sequence.
 */
void localalign(const std::string &x,
                const std::string &y,
                const std::string &header1,
                const std::string &header2,
                const std::string &outdir,
                ScoreMode mode,
                ScoreFn score_fn)
{
    (void)mode;
    using Clock = std::chrono::high_resolution_clock;
    auto t_start = Clock::now();

    // 1) Striped‐SW forward pass + banded‐traceback
    std::string alignedX, alignedY;
    std::vector<std::pair<int,int>> local_path;
    stripedSmithWaterman(
        x.data(), (int)x.size(),
        y.data(), (int)y.size(),
        (int)GAP_OPEN, (int)GAP_OPEN,
        score_fn,
        alignedX, alignedY, local_path
    );

    // 2) Timing
    auto t_end   = Clock::now();
    auto time_ms = std::chrono::duration_cast<std::chrono::milliseconds>(t_end - t_start).count();

    // 3) Compute stats
    size_t total = alignedX.size(), gaps = 0, matches = 0;
    for (size_t i = 0; i < total; ++i) {
        if (alignedX[i]=='-' || alignedY[i]=='-')      ++gaps;
        else if (alignedX[i]==alignedY[i])            ++matches;
    }
    double identity = double(matches) / total;
    double coverage = double(total - gaps) / total;
    int bestScore = matches - gaps * int(-GAP_OPEN);
      // (or store maxScore from striped pass if you capture it)

    // 4) Verbose console output
    if (verbose) {
        std::cout << "\n\nLocal Alignment Score: " << bestScore << "\n"
                  << "Gap Open: "   << GAP_OPEN   << "\n"
                  << "Matches: "    << matches    << "\n"
                  << "Gaps:    "    << gaps       << "\n"
                  << "Total:   "    << total      << "\n"
                  << "Identity: "   << identity*100.0 << "%\n"
                  << "Coverage: "   << coverage*100.0 << "%\n"
                  << "Time:    "    << time_ms    << " ms\n"
                  << "Query:   "    << header1    << "\n"
                  << "Target:  "    << header2    << "\n";
        printColoredAlignment(alignedX, alignedY);
    }

    // 5) Write local_path.txt
    {
      std::ofstream pf(outdir + "/local_path.txt");
      for (auto &p : local_path)
        pf << p.first << " " << p.second << "\n";
    }

    // 6) Write FASTA alignment
    {
      std::ofstream outf(outdir + "/local_alignment.fasta");
      savePlainAlignment(header1, header2, alignedX, alignedY, outf);
    }

    // 7) Write JSON stats
    {
      std::ofstream js(outdir + "/local_stats.json");
      js << std::fixed << std::setprecision(6)
         << "{\n"
         << "  \"method\":   \"local\",\n"
         << "  \"gap_open\": "   << GAP_OPEN   << ",\n"
         << "  \"score\":    "   << bestScore  << ",\n"
         << "  \"matches\":  "   << matches    << ",\n"
         << "  \"gaps\":     "   << gaps       << ",\n"
         << "  \"total\":    "   << total      << ",\n"
         << "  \"identity\": "   << identity   << ",\n"
         << "  \"coverage\": "   << coverage   << ",\n"
         << "  \"time_ms\":  "   << time_ms    << ",\n"
         << "  \"query\":    \"" << header1    << "\",\n"
         << "  \"target\":   \"" << header2    << "\"\n"
         << "}\n";
    }
}


/**
 * @brief Compute the Longest Common Subsequence (LCS) of two strings.
 *
 * Uses a classic DP with blocking via AVX2 for byte‐wise comparisons
 * (32‐byte lanes) and OpenMP to parallelize across rows.
 *
 * @param x First sequence.
 * @param y Second sequence.
 */
void lcs(const string &x, const string &y,
         const string &header1, const string &header2,
         const std::string &outdir, ScoreMode mode) {
    const int m = x.size(), n = y.size();
    vector<int> prev(n+1), curr(n+1);
    vector<vector<char>> b(m+1, vector<char>(n+1,' '));
    vector<pair<int,int>> lcs_path;

    // bytes per AVX2 register
    const int B = 32;
    const int j_limit = n - (n % B);

    #pragma omp parallel
    {
        for (int i = 1; i <= m; ++i) {
            __m256i vx = _mm256_set1_epi8(x[i-1]);
            // process full 32-byte blocks
            #pragma omp for
            for (int j = 1; j <= j_limit; j += B) {
                __m256i vy   = _mm256_loadu_si256((__m256i const*)(y.data() + j - 1));
                __m256i eq   = _mm256_cmpeq_epi8(vx, vy);
                uint32_t mask = _mm256_movemask_epi8(eq);
                for (int k = 0; k < B; ++k) {
                    int jj = j + k;

                    // mask bit k == 1 means x[i-1]==y[jj-1]
                    if ((mask >> k) & 1) {
                        curr[jj] = prev[jj-1] + 1;
                        b[i][jj] = 'D';
                    } else if (prev[jj] >= curr[jj-1]) {
                        curr[jj] = prev[jj];
                        b[i][jj] = 'U';
                    } else {
                        curr[jj] = curr[jj-1];
                        b[i][jj] = 'L';
                    }
                }
            }

            // fallback for the tail
            #pragma omp for
            for (int j = j_limit+1; j <= n; ++j) {
                if (x[i-1] == y[j-1]) {
                    curr[j] = prev[j-1] + 1;
                    b[i][j] = 'D';
                } else if (prev[j] >= curr[j-1]) {
                    curr[j] = prev[j];
                    b[i][j] = 'U';
                } else {
                    curr[j] = curr[j-1];
                    b[i][j] = 'L';
                }
            }

            #pragma omp single
            {
                swap(prev, curr);
                if (verbose) {if (i % 1000 == 0 || i == m) showProgressBar(i, m);
                }
            }
        }
    }

    // traceback
    int i = m, j = n;
    string lcs_str;
    while (i > 0 && j > 0) {
        if (b[i][j] == 'D') {
        lcs_path.emplace_back(j, i);
        lcs_str += x[i-1]; --i; --j; }
        else if (b[i][j] == 'U')  { --i; }
        else                      { --j; }
    }
    reverse(lcs_str.begin(), lcs_str.end());
    reverse(lcs_path.begin(), lcs_path.end());

    {
      ofstream pf(outdir + "/lcs_path.txt");
      for (auto &p : lcs_path)
        pf << p.first << " " << p.second << "\n";
    }

    std::ofstream outfile(outdir +"/lcs.fasta");

    if (binary) {
        writeCharMatrix(b, outdir +"/lcs_traceback.bin");
    } else if (txt) {
        writeRawCharMatrix(b, outdir +"/lcs_traceback.txt");
    } else {
        ;
    }

    // write to fasta
    if (!outfile) throw runtime_error("Cannot open lcs.fasta for writing");
    saveLCS(getAccession(header1,mode) + "_" + getAccession(header2,mode),
            lcs_str,
            outfile);
    outfile.close();

    // print to console in chunks:
    if (verbose) {
      std::cout << "\nLCS length: " << lcs_str.size() << "\n\n";
      for (size_t i = 0; i < lcs_str.size(); i += LINE_WIDTH) {
        std::cout << lcs_str.substr(i, LINE_WIDTH) << "\n";
      }
    }
}


/**
* @brief Program entry point: read arguments, load FASTA, dispatch chosen method.
*
* @param argc Number of command-line arguments.
* @param argv Command-line arguments.
* @return 0 on success, non-zero on error.
*/
int main(int argc, char** argv) {
    MPI_Init(&argc, &argv);
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    try {
        std::string file1, file2, outdir = ".";
        int choice = -1;
        ScoreMode mode = MODE_DNA;
        bool useGz = false;
        std::string targetChr;

        // Parse arguments
        for (int i = 1; i < argc; ++i) {
            std::string arg = argv[i];
            if (arg == "--query" && i + 1 < argc) {
                file1 = argv[++i];
            } else if (arg == "--target" && i + 1 < argc) {
                file2 = argv[++i];
            } else if (arg == "--chrom" && i+1<argc) {
                targetChr = argv[++i];
            } else if (arg == "--gzipped" || arg == "--gz") {
                useGz = true;
            } else if (arg == "--choice" && i + 1 < argc) {
                choice = std::stoi(argv[++i]);
            } else if (arg == "--mode" && i + 1 < argc) {
                std::string val = argv[++i];
                if (val == "dna") mode = MODE_DNA;
                else if (val == "protein") mode = MODE_PROTEIN;
                else {
                    if (rank == 0) std::cerr << "Unknown mode: " << val << "\n";
                    MPI_Finalize(); return 1;
                }
            } else if (arg == "--outdir" && i + 1 < argc) {
                outdir = argv[++i];
            } else if (arg == "--verbose") {
                verbose = true;
            } else if (arg == "--binary") {
                binary = true;
            } else if (arg == "--txt") {
                txt = true;
            } else if (arg == "--help") {
                if (rank == 0) {
                    std::cout <<
                      "Usage: ./aligner --query <file1> --target <file2> --choice <1|2|3|4> [--mode dna|protein] [--outdir DIR] [--verbose]\n"
                      "  --choice: 1=global, 2=local, 3=LCS, 4=all\n";
                }
                MPI_Finalize(); return 0;
            } else if (arg == "--gap_open" && i + 1 < argc) {
                 GAP_OPEN = std::stod(argv[++i]);
            } else if (arg == "--gap_extend" && i + 1 < argc) {
                 GAP_EXTEND = std::stod(argv[++i]);
            } else {
                if (rank == 0) std::cerr << "Unknown option: " << arg << "\n";
                MPI_Finalize(); return 1;
            }
        }

        ScoreFn score_fn = nullptr;
        if (mode == MODE_DNA) {
            score_fn = &edna_score;
        } else {
            score_fn = &blosum62_score;
        }

        if (file1.empty() || file2.empty() || choice == -1) {
            if (rank == 0)
                std::cerr << "Missing required arguments: --query, --target, --choice\n";
            MPI_Finalize(); return 1;
        }

        if (rank == 0) {
            std::filesystem::create_directories(outdir);
        }

        std::string seq1, seq2, header1, header2;
        if (rank==0) {
            if (useGz) {
                processFastaGz(file1, header1, seq1);
                processFastaGz(file2, header2, seq2, targetChr);
            } else {
                processFasta(file1, header1, seq1);
                processFasta(file2, header2, seq2, targetChr);
            }
        }

        int len1 = seq1.size(), len2 = seq2.size();
        MPI_Bcast(&len1, 1, MPI_INT, 0, MPI_COMM_WORLD);
        MPI_Bcast(&len2, 1, MPI_INT, 0, MPI_COMM_WORLD);

        if (rank != 0) {
            seq1.resize(len1);
            seq2.resize(len2);
        }

        MPI_Bcast(seq1.data(), len1, MPI_CHAR, 0, MPI_COMM_WORLD);
        MPI_Bcast(seq2.data(), len2, MPI_CHAR, 0, MPI_COMM_WORLD);

        if (choice == 1) globalalign(seq1, seq2, header1, header2, outdir, mode, score_fn);
        else if (choice == 2) localalign(seq1, seq2, header1, header2, outdir, mode, score_fn);
        else if (choice == 3) lcs(seq1, seq2, header1, header2, outdir, mode);
        else if (choice == 4) {
            globalalign(seq1, seq2, header1, header2, outdir, mode, score_fn);
            MPI_Barrier(MPI_COMM_WORLD);
            localalign(seq1, seq2, header1, header2, outdir, mode, score_fn);
            MPI_Barrier(MPI_COMM_WORLD);
            lcs(seq1, seq2, header1, header2, outdir, mode);
            MPI_Barrier(MPI_COMM_WORLD);
        } else if (rank == 0) {
            std::cerr << "Invalid method. Use --choice 1/2/3/4.\n";
        }

    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << "\n";
    }

    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Finalize();
    return 0;
}
