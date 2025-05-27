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
enum ScoreMode { MODE_DNA, MODE_PROTEIN };
using ScoreFn = int(*)(char,char);

/// Gap penalty
double GAP_OPEN   = -5.0;  // penalty to open a gap
double GAP_EXTEND = -1.0;  // penalty to extend an existing gap in global alignment
/// Number of columns per line when printing alignments
static const int LINE_WIDTH = 80;

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
 * all subsequent lines into a single sequence string.
 *
 * @param filename  Path to the FASTA file.
 * @param[out] header   On return, the header line without leading `>`.
 * @param[out] sequence On return, the full sequence (no newlines).
 * @throws std::runtime_error if the file cannot be opened.
 */
void processFasta(const string &filename, string &header, string &sequence) {
    ifstream file(filename);
    if (!file) throw runtime_error("Error: Unable to open " + filename);
    string line;
    header = "";
    sequence = "";
    bool headerSet = false;
    while (getline(file, line)) {
        if (line.empty()) continue;
        if (line[0] == '>') {
            if (!headerSet) {
                header = line.substr(1);
                headerSet = true;
            }
            continue;
        }
        sequence += line;
    }
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
 * @brief Print the aligned sequences with color coding.
 *
 * This function prints two aligned sequences side by side, with color coding:
 * - Green for matches
 * - Red for gaps
 * - Cyan for mismatches
 *
 * @param seq1 The first aligned sequence.
 * @param seq2 The second aligned sequence.
 * @param os   Output stream (default is std::cout).
 */
void printColoredAlignment(const string &seq1, const string &seq2, ostream &os = cout) {
    size_t aln_length = seq1.size(); // Assuming seq1 and seq2 are already aligned and have same length
    if (aln_length == 0) {
        os << "No alignment to print.\n";
        return;
    }

    size_t pos1_ungapped_count = 0; // Running count of non-gap characters from seq1
    size_t pos2_ungapped_count = 0; // Running count of non-gap characters from seq2

    for (size_t i = 0; i < aln_length; i += LINE_WIDTH) {
        size_t end_block = min(i + LINE_WIDTH, aln_length); // End of the current block in the alignment

        // Store the ungapped start positions for this block
        size_t block_start_pos1 = pos1_ungapped_count + 1;
        size_t block_start_pos2 = pos2_ungapped_count + 1;

        // Temporary counters for ungapped characters within the current block,
        // to correctly display end positions for the block lines.
        size_t current_block_end_pos1 = pos1_ungapped_count;
        size_t current_block_end_pos2 = pos2_ungapped_count;

        // Print Sequence 1 part of the block
        os << std::setw(6) << block_start_pos1 << " ";
        for (size_t j = i; j < end_block; ++j) {
            if (seq1[j] == seq2[j]) os << GREEN << seq1[j] << RESET;
            else if (seq1[j] == '-' || seq2[j] == '-') os << RED << seq1[j] << RESET;
            else os << CYAN << seq1[j] << RESET;
            if (seq1[j] != '-') {
                current_block_end_pos1++;
            }
        }
        os << " " << current_block_end_pos1 << "\n";

        // Print Sequence 2 part of the block
        os << std::setw(6) << block_start_pos2 << " ";
        for (size_t j = i; j < end_block; ++j) {
            if (seq1[j] == seq2[j]) os << GREEN << seq2[j] << RESET;
            else if (seq1[j] == '-' || seq2[j] == '-') os << RED << seq2[j] << RESET; // Note: color based on seq1[j] vs seq2[j]
            else os << CYAN << seq2[j] << RESET;
            if (seq2[j] != '-') {
                current_block_end_pos2++;
            }
        }
        os << " " << current_block_end_pos2 << "\n";

        os << "\n"; // Add an extra newline between blocks

        // Update the master ungapped counts
        pos1_ungapped_count = current_block_end_pos1;
        pos2_ungapped_count = current_block_end_pos2;
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
**/
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

/**
 * @brief Initialize the DP structures for affine gap scoring.
 *
 * @param n         Length of the second sequence (Y).
 * @param prev_row  Row i-1, match/mismatch scores.
 * @param prev_gapX Row i-1, gap scores in X.
 * @param prev_gapY Row i-1, gap scores in Y.
 * @param isGlobal  Whether this is a global alignment (affects initialization).
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


/**
 * @brief Compute a single row of the affine DP matrix.
 *
 * This function computes the i-th row of the DP matrix for affine gap scoring,
 * including match/mismatch scores and gap scores in both directions.
 *
 * @param i         The current row index (1-based).
 * @param x         The first sequence (string).
 * @param y         The second sequence (string).
 * @param prev_row  Previous row S_{i-1}.
 * @param prev_gapX Previous gap row F_{i-1} (gap in Y).
 * @param prev_gapY Previous gap row E_{i-1} (gap in X).
 * @param curr_row  Current row S_i to be filled.
 * @param curr_gapX Current gap row F_i to be filled.
 * @param curr_gapY Current gap row E_i to be filled.
 * @param curr_trace_row Pointer/Trace for S_i to be filled.
 * @param score_fn Scoring function for match/mismatch.
 */
void computeAffineDPRow(int i,
                        const string& x,
                        const string& y,
                        vector<int>& prev_row,    // S_{i-1}
                        vector<int>& prev_gapX,   // F_{i-1} (Gap in Y, vertical)
                        vector<int>& prev_gapY,   // E_{i-1} (Gap in X, horizontal)
                        vector<int>& curr_row,    // S_i
                        vector<int>& curr_gapX,   // F_i
                        vector<int>& curr_gapY,   // E_i
                        vector<char>& curr_trace_row, // Pointer/Trace for S_i
                        ScoreFn score_fn)
{
    // (void)prev_gapY; // prev_gapY is actually used for curr_gapY calculation if it's E_{i-1}
    int n = y.size();

    curr_row.resize(n+1);
    curr_gapX.resize(n+1);
    curr_gapY.resize(n+1);
    curr_trace_row.resize(n+1); // Resize trace row

    // Initialize first elements of current row/gaps/trace
    // curr_row[0] (S_i,0)
    // curr_gapX[0] (F_i,0 - Gap in Y)
    // curr_gapY[0] (E_i,0 - Gap in X) - typically NEG_INF for global

    // For S_i,0 (all x_k vs gaps)
    // This depends on global F_{i,0} which is S_{i-1,0} + GAP_OPEN or F_{i-1,0} + GAP_EXTEND
    if (i > 0) { // Only for rows i > 0
        // F_{i,0} comes from S_{i-1,0} (prev_row[0]) or F_{i-1,0} (prev_gapX[0])
        int f_open = prev_row[0] + GAP_OPEN;
        int f_extend = prev_gapX[0] + GAP_EXTEND;
        curr_gapX[0] = max(f_open, f_extend);
        curr_row[0] = curr_gapX[0]; // S_{i,0} is F_{i,0} for global

        // Determine pointer for (i,0) based on F_open vs F_extend
        if (curr_gapX[0] == f_open && curr_gapX[0] >= f_extend) { // Prefer opening if scores are equal
             curr_trace_row[0] = 'F'; // Indicates F state, from S_{i-1,0} (vertical move)
        } else {
             curr_trace_row[0] = 'f'; // Indicates F state, from F_{i-1,0} (vertical extend)
        }
        // E_{i,0} is usually -infinity for global alignment as you can't have a gap in X before Y starts.
        curr_gapY[0] = INT_MIN / 2;
    }


    for (int j = 1; j <= n; ++j) {
        // Calculate F_{i,j} (gap in Y, score ending with x_i vs '-')
        // F_{i,j} = max(S_{i-1,j} + GAP_OPEN, F_{i-1,j} + GAP_EXTEND)
        int f_open_score = prev_row[j] + GAP_OPEN;
        int f_extend_score = prev_gapX[j] + GAP_EXTEND;
        curr_gapX[j] = max(f_open_score, f_extend_score);

        // Calculate E_{i,j} (gap in X, score ending with y_j vs '-')
        // E_{i,j} = max(S_{i,j-1} + GAP_OPEN, E_{i,j-1} + GAP_EXTEND)
        int e_open_score = curr_row[j-1] + GAP_OPEN; // curr_row[j-1] is S_{i,j-1}
        int e_extend_score = curr_gapY[j-1] + GAP_EXTEND; // curr_gapY[j-1] is E_{i,j-1}
        curr_gapY[j] = max(e_open_score, e_extend_score);

        // Calculate M_{i,j} (match/mismatch score)
        // M_{i,j} = max(S_{i-1,j-1}, E_{i-1,j-1}, F_{i-1,j-1}) + score_fn(x[i-1],y[j-1])
        int mscore = score_fn(x[i-1], y[j-1]);
        int diag_predecessor_score = max({prev_row[j-1], prev_gapY[j-1], prev_gapX[j-1]}); // Max of S, E, F from top-left
        int match_val = diag_predecessor_score + mscore;

        // Determine S_{i,j} and the pointer
        // Preference: Match > Gap in X (horizontal E) > Gap in Y (vertical F) - adjust as needed
        if (match_val >= curr_gapY[j] && match_val >= curr_gapX[j]) {
            curr_row[j] = match_val;
            curr_trace_row[j] = 'M'; // Match/Mismatch
        } else if (curr_gapY[j] >= curr_gapX[j]) { // Prefer E over F if match_val is not the max
            curr_row[j] = curr_gapY[j];
            // Distinguish if E came from open or extend
            if (curr_gapY[j] == e_open_score && curr_gapY[j] >= e_extend_score) { // Prefer open
                 curr_trace_row[j] = 'E'; // Gap in X, opened from S_{i,j-1}
            } else {
                 curr_trace_row[j] = 'e'; // Gap in X, extended from E_{i,j-1}
            }
        } else {
            curr_row[j] = curr_gapX[j];
            // Distinguish if F came from open or extend
            if (curr_gapX[j] == f_open_score && curr_gapX[j] >= f_extend_score) { // Prefer open
                curr_trace_row[j] = 'F'; // Gap in Y, opened from S_{i-1,j}
            } else {
                curr_trace_row[j] = 'f'; // Gap in Y, extended from F_{i-1,j}
            }
        }
    }
}


/**
 * @brief Perform global sequence alignment using the Needleman-Wunsch algorithm.
 *
 * This function computes the global alignment of two sequences using affine gap scoring.
 * It initializes the DP matrix, fills it row by row, and performs traceback to get the aligned sequences.
 *
 * @param x         The first sequence (string).
 * @param y         The second sequence (string).
 * @param header1   Header for the first sequence.
 * @param header2   Header for the second sequence.
 * @param outdir    Output directory for results.
 * @param mode      Scoring mode (MODE_DNA or MODE_PROTEIN).
 * @param score_fn  Scoring function for match/mismatch.
 */
void globalalign(const string &x, const string &y,
                 const string &header1, const string &header2,
                 const std::string &outdir, ScoreMode mode, ScoreFn score_fn) {
    int m = x.size(), n = y.size();

    vector<int> prev_row, prev_gapX, prev_gapY;
    vector<int> curr_row, curr_gapX, curr_gapY;
    vector<char> curr_trace_row; // For pointers of the current row

    // Initialize DP structures for row 0
    // S_0,j, E_0,j, F_0,j
    prev_row.assign(n + 1, 0);
    prev_gapX.assign(n + 1, INT_MIN / 2); // F_0,j (gap in Y)
    prev_gapY.assign(n + 1, INT_MIN / 2); // E_0,j (gap in X)

    vector<vector<char>> fullTraceMatrix(m + 1, vector<char>(n + 1)); // To store all pointers

    // Initialization for S[0][0], E[0][0], F[0][0]
    prev_row[0] = 0;
    // E[0][0] and F[0][0] are usually NEG_INF or handled by context.
    // Let's ensure they are deeply negative if not 0.
    prev_gapX[0] = INT_MIN / 2;
    prev_gapY[0] = INT_MIN / 2;
    fullTraceMatrix[0][0] = 'S'; // Start or Stop

    // Initialize first row (i=0) S_0,j, E_0,j, F_0,j
    // For global alignment, S_0,j must be a gap in X (sequence Y vs gaps)
    for (int j = 1; j <= n; ++j) {
        prev_gapY[j] = (j == 1) ? (prev_row[j-1] + GAP_OPEN) : (prev_gapY[j-1] + GAP_EXTEND); // E_0,j from S_0,j-1 or E_0,j-1
                                                                                         // S_0,j-1 is prev_row[j-1]
        prev_row[j] = prev_gapY[j]; // S_0,j = E_0,j
        prev_gapX[j] = INT_MIN / 2; // F_0,j is not typically possible
        fullTraceMatrix[0][j] = (j==1 && prev_gapY[j] == prev_row[j-1] + GAP_OPEN) ? 'E' : 'e';
    }


    // For storing fullDP if needed for output (optional if only alignment needed)
    vector<vector<int>> fullDP(m + 1, vector<int>(n + 1));
    for(int j=0; j<=n; ++j) fullDP[0][j] = prev_row[j];


    using Clock = std::chrono::high_resolution_clock;
    auto t_start = Clock::now();

    for (int i = 1; i <= m; ++i) {
        // computeAffineDPRow now fills curr_trace_row as well
        computeAffineDPRow(i, x, y,
                           prev_row, prev_gapX, prev_gapY,
                           curr_row, curr_gapX, curr_gapY,
                           curr_trace_row, // Pass the trace row to be filled
                           score_fn);

        // Store the trace row
        for(int j=0; j<=n; ++j) fullTraceMatrix[i][j] = curr_trace_row[j];
        // Store the S scores if you need the fullDP matrix later
        for(int j=0; j<=n; ++j) fullDP[i][j] = curr_row[j];


        // Swap for next iteration
        prev_row.swap(curr_row);
        prev_gapX.swap(curr_gapX);
        prev_gapY.swap(curr_gapY);
        // No need to swap curr_trace_row as it's freshly computed each time based on S,E,F of previous.

        if (verbose) { if (i % 100 == 0 || i == m) { // Adjusted progress update frequency
            showProgressBar(i, m);
        }}
    }

    // Score is in the last element of the S matrix (which is now in prev_row after the loop)
    int final_score = prev_row[n];


    // TRACEBACK using fullTraceMatrix
    string alignedX, alignedY;
    vector<pair<int,int>> global_path; // Declare/use global_path variable
    int cur_i = m, cur_j = n;

    // Add the starting point of the traceback (end of alignment)
    // The path is usually traced from (m,n) back to (0,0) or where alignment starts
    global_path.emplace_back(cur_j, cur_i);

    while (cur_i > 0 || cur_j > 0) {
        char trace_char = fullTraceMatrix[cur_i][cur_j];

        if (trace_char == 'M') { // Match/Mismatch
            alignedX += x[cur_i - 1];
            alignedY += y[cur_j - 1];
            cur_i--;
            cur_j--;
        } else if (trace_char == 'F' || trace_char == 'f') { // Gap in Y (x_i vs '-'), came from F state
            alignedX += x[cur_i - 1];
            alignedY += '-';
            cur_i--;
        } else if (trace_char == 'E' || trace_char == 'e') { // Gap in X (y_j vs '-'), came from E state
            alignedX += '-';
            alignedY += y[cur_j - 1];
            cur_j--;
        } else { // Boundary conditions for i=0 or j=0
            if (cur_i == 0 && cur_j > 0) { // Reached top row, only gaps in X possible
                 alignedX += '-';
                 alignedY += y[cur_j - 1];
                 cur_j--;
            } else if (cur_j == 0 && cur_i > 0) { // Reached first col, only gaps in Y possible
                 alignedX += x[cur_i - 1];
                 alignedY += '-';
                 cur_i--;
            } else {
                // Should not happen if matrix is well-formed and traceback starts from m,n
                // Or if fullTraceMatrix[0][0] was 'S' and we hit it.
                if (cur_i == 0 && cur_j == 0 && fullTraceMatrix[0][0] == 'S') break;
                // cerr << "Traceback error: unexpected pointer " << trace_char << " at " << cur_i << "," << cur_j << endl;
                break;
            }
        }
        // Add the current position to the global path
        global_path.emplace_back(cur_j, cur_i);
    }

    auto t_end = Clock::now();
    auto time_ms = std::chrono::duration_cast<std::chrono::milliseconds>(t_end - t_start).count();

    if (binary) {
        writeDPMatrix(fullDP, outdir +"/global_dp_matrix.bin");
    } else if (txt) {
        writeRawDPMatrix(fullDP, outdir +"/global_dp_matrix.txt");
    } else {
        ;
    }

    {
      ofstream pf(outdir + "/global_path.txt");
      for (auto &p : global_path)
        pf << p.first << " " << p.second << "\n";
    }

    size_t total   = alignedX.size();
    size_t gaps    = 0, matches = 0;
    for (size_t i = 0; i < total; ++i) {
        if (alignedX[i]=='-' || alignedY[i]=='-') ++gaps;
        else if (alignedX[i]==alignedY[i])       ++matches;
    }
    double identity = double(matches) / total;
    double coverage = double(total - gaps) / total;

    std::string accession1 = getAccession(header1, mode);
    std::string accession2 = getAccession(header2, mode);

    std::string gene1 = getGeneSymbol(header1, mode);
    std::string gene2 = getGeneSymbol(header2, mode);

    reverse(alignedX.begin(), alignedX.end());
    reverse(alignedY.begin(), alignedY.end());

    if (verbose) {
        cout << "\n\nGlobal Alignment Score: " << final_score << "\n";
        cout << "Gap Open: " << GAP_OPEN << "\n";
        cout << "Gap Extend: " << GAP_EXTEND << "\n";
        cout << "Matches: " << matches << "\n";
        cout << "Gaps:    " << gaps << "\n";
        cout << "Total:   " << total << "\n";
        cout << "Identity: " << identity * 100.0f << "%\n";
        cout << "Coverage: " << coverage * 100.0f << "%\n";
        cout << "Time:    " << time_ms << " ms\n";
        cout << "Query:   " << accession1 << "\n";
        cout << "Target:  " << accession2 << "\n";
        cout << "QueryID:  " << gene1 << "\n";
        cout << "TargetID:  " << gene2 << "\n\n\n";
        printColoredAlignment(alignedX, alignedY);
    }

    std::ofstream outfile(outdir +"/global_alignment.fasta");
    if (outfile) {
        savePlainAlignment(getAccession(header1,mode),
                           getAccession(header2,mode),
                           alignedX, alignedY,
                           outfile);
            outfile.close();
    } else {
        cerr << "Error: Unable to open output file global_alignment.fasta\n";
    }

    ofstream js(outdir +"/global_stats.json");
    if (js) {
      js << fixed << setprecision(6)
         << "{\n"
         << "  \"method\":      \"global\",\n"
         << " \"gap_open\":   " << GAP_OPEN << ",\n"
         << "  \"gap_extend\": " << GAP_EXTEND << ",\n"
         << "  \"score\":       " << prev_row[n]  << ",\n"
         << "  \"matches\":     " << matches << ",\n"
         << "  \"gaps\":        " << gaps << ",\n"
         << "  \"total\":       " << total << ",\n"
         << "  \"identity\":    " << identity << ",\n"
         << "  \"coverage\":    " << coverage << ",\n"
         << "  \"time_ms\":     " << time_ms << ",\n"
         << "  \"query\":       \"" << accession1 << "\",\n"
         << "  \"target\":      \"" << accession2 << "\",\n"
         << "  \"queryid\":       \"" << gene1 << "\",\n"
         << "  \"targetid\":       \"" << gene2 << "\",\n"
         << "  \"query_length\": " << m << ",\n"
         << "  \"target_length\": " << n << "\n"
         << "}\n";
      js.close();
    } else {
      cerr << "Error: cannot open global_stats.json\n";
    }

}


// Struct to hold scores and pointers for a single DP cell in affine local alignment
struct AffineDPScores {
    int s_val = 0;
    int e_val = 0; // Score ending with gap in query (seqX), horizontal move
    int f_val = 0; // Score ending with gap in target (seqY), vertical move
    char ptr = 'X'; // Pointer for S_val traceback: 'M', 'E'/'e' (open/extend E), 'F'/'f' (open/extend F), 'X' (stop/zero)
};

/**
 * @brief Compute a single cell in the local affine DP matrix.
 *
 * This function computes the scores for a single cell in the local alignment DP matrix
 * using affine gap penalties. It returns an AffineDPScores struct containing the scores
 * and the pointer for traceback.
 *
 * @param s_diag_prev  S value from the diagonal cell (i-1, j-1)
 * @param e_diag_prev  E value from the diagonal cell (i-1, j-1)
 * @param f_diag_prev  F value from the diagonal cell (i-1, j-1)
 * @param s_left       S value from the left cell (i, j-1)
 * @param e_left       E value from the left cell (i, j-1)
 * @param s_up         S value from the upper cell (i-1, j)
 * @param f_up         F value from the upper cell (i-1, j)
 * @param char_x       Character from sequence X at position i
 * @param char_y       Character from sequence Y at position j
 * @param score_fn     Scoring function for match/mismatch
 * @return AffineDPScores containing computed scores and pointer
 */
AffineDPScores compute_local_affine_cell(
    int s_diag_prev, int e_diag_prev, int f_diag_prev, // S, E, F of (i-1, j-1)
    int s_left, int e_left,                             // S, E of (i, j-1)
    int s_up, int f_up,                               // S, F of (i-1, j)
    char char_x, char char_y, ScoreFn score_fn) {

    AffineDPScores cell_scores;

    // Gap penalties (ensure these are negative scores from command line)
    // double local_gap_open = GAP_OPEN;
    // double local_gap_extend = GAP_EXTEND;
    // In C++, GAP_OPEN and GAP_EXTEND are already global and negative.

    int match_mismatch_score = score_fn(char_x, char_y);

    // M_val: score from match/mismatch
    int m_val_pred = std::max({s_diag_prev, e_diag_prev, f_diag_prev});
    int m_val = m_val_pred + match_mismatch_score;

    // E_val: score ending with gap in query X (horizontal move)
    int e_open_score = s_left + GAP_OPEN;
    int e_extend_score = e_left + GAP_EXTEND;
    int current_E_val = std::max(e_open_score, e_extend_score);

    // F_val: score ending with gap in target Y (vertical move)
    int f_open_score = s_up + GAP_OPEN;
    int f_extend_score = f_up + GAP_EXTEND;
    int current_F_val = std::max(f_open_score, f_extend_score);

    // Apply floor of 0 for Smith-Waterman
    m_val = std::max(0, m_val);
    current_E_val = std::max(0, current_E_val);
    current_F_val = std::max(0, current_F_val);

    cell_scores.e_val = current_E_val;
    cell_scores.f_val = current_F_val;

    // Determine S_val and pointer, preference M > E > F
    if (m_val >= current_E_val && m_val >= current_F_val) {
        cell_scores.s_val = m_val;
        cell_scores.ptr = (m_val > 0) ? 'M' : 'X';
    } else if (current_E_val >= current_F_val) {
        cell_scores.s_val = current_E_val;
        if (current_E_val > 0) {
             cell_scores.ptr = (current_E_val == e_open_score && current_E_val >= e_extend_score) ? 'E' : 'e';
        } else {
            cell_scores.ptr = 'X';
        }
    } else {
        cell_scores.s_val = current_F_val;
        if (current_F_val > 0) {
            cell_scores.ptr = (current_F_val == f_open_score && current_F_val >= f_extend_score) ? 'F' : 'f';
        } else {
            cell_scores.ptr = 'X';
        }
    }
    // Final S score must also be >= 0 (implicitly handled if M,E,F are already floored and S is max of them)
    // If all M, E, F were 0, then S is 0.
    if (cell_scores.s_val == 0 && (m_val > 0 || current_E_val > 0 || current_F_val > 0) ) {
        // This case should not happen if M, E, F are already max(0, val) and S = max(M,E,F)
        // However, if S is strictly max(0, M, E, F) as one combined step:
        // cell_scores.s_val = std::max({0, m_val_unfloored, current_E_val_unfloored, current_F_val_unfloored});
        // For clarity, ensure S is max of already-floored M, E, F, and also >=0
    }
     if (cell_scores.s_val == 0) cell_scores.ptr = 'X';


    return cell_scores;
}

// Struct to hold the best local alignment score and its position
struct Loc {
    int score;
    int i; // Represents global row index (1-indexed for matrix cell)
    int j; // Represents global column index (1-indexed for matrix cell)
};

/**
 * @brief Perform local sequence alignment using the Smith-Waterman algorithm.
 *
 * This function computes the local alignment of two sequences using affine gap scoring.
 * It initializes the DP matrix, fills it row by row, and performs traceback to get the aligned sequences.
 *
 * @param x         The first sequence (string).
 * @param y         The second sequence (string).
 * @param header1   Header for the first sequence.
 * @param header2   Header for the second sequence.
 * @param outdir    Output directory for results.
 * @param mode      Scoring mode (MODE_DNA or MODE_PROTEIN).
 * @param score_fn  Scoring function for match/mismatch.
 */
void localalign(const std::string &x,
                const std::string &y,
                const std::string &header1,
                const std::string &header2,
                const std::string &outdir,
                ScoreMode mode,
                ScoreFn score_fn)
{
    int m = x.size(), n = y.size();
    int rank, num_procs;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &num_procs);


    // Each rank will now need to manage S, E, F rows.
    // For simplicity, let's assume each rank computes its S, E, F rows independently for its block,
    // receiving the necessary S, E, F values from the previous rank's last computed row.

    int chunkSize = (m + num_procs - 1) / num_procs;
    int start_row_global = rank * chunkSize; // 0-indexed global start row for this rank
    int end_row_global = std::min(start_row_global + chunkSize, m); // 0-indexed global end row (exclusive)
    int local_num_rows = end_row_global - start_row_global;

    // DP matrices for the current rank's block (S, E, F)
    // +1 for the initial row received from prev rank or initialized to 0s

    std::vector<std::vector<int>> s_block(local_num_rows + 1, std::vector<int>(n + 1, 0));
    std::vector<std::vector<int>> e_block(local_num_rows + 1, std::vector<int>(n + 1, 0));
    std::vector<std::vector<int>> f_block(local_num_rows + 1, std::vector<int>(n + 1, 0));

    // Pointers are only strictly needed for the final traceback on bestRank,
    // but could be stored per block if memory allows and helps.

    // Initialization for local alignment: S[0,j]=0, E[0,j]=0, F[0,j]=0. S[i,0]=0, E[i,0]=0, F[i,0]=0.
    // S[0,0]=0. E[0,0]=NEG_INF, F[0,0]=NEG_INF (conceptually for transitions into first cell).
    // We initialize S,E,F edges to 0. Let's stick to that for simplicity.

    if (rank > 0) {

        // Receive S, E, F values of the last row from rank-1
        // These become the 'previous' row (index 0 in s_block, e_block, f_block)

        MPI_Recv(s_block[0].data(), n + 1, MPI_INT, rank - 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Recv(e_block[0].data(), n + 1, MPI_INT, rank - 1, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Recv(f_block[0].data(), n + 1, MPI_INT, rank - 1, 2, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    } else {

        // Rank 0 initializes its 'previous' row (s_block[0], e_block[0], f_block[0]) to all zeros.
        // s_block[0] is already all zeros by vector initialization.
        // e_block[0] and f_block[0] should also be all zeros for local affine.
    }

    using Clock = std::chrono::high_resolution_clock;
    auto t_start = Clock::now();

    // score, global_i (0-indexed), global_j (0-indexed)
    Loc localBest{0, 0, 0};

    // Fill DP block using affine gap penalties
    for (int ii = 1; ii <= local_num_rows; ++ii) { // ii is 1-indexed local row
        int global_i_idx = start_row_global + ii - 1; // 0-indexed global character from x

        // S[i,0], E[i,0], F[i,0] are 0 for local
        s_block[ii][0] = 0;
        e_block[ii][0] = 0;
        f_block[ii][0] = 0;

        for (int j = 1; j <= n; ++j) { // j is 1-indexed column
            int global_j_idx = j - 1; // 0-indexed global character from y

            AffineDPScores cell = compute_local_affine_cell(
                s_block[ii-1][j-1], e_block[ii-1][j-1], f_block[ii-1][j-1], // S,E,F from diag_prev
                s_block[ii][j-1],   e_block[ii][j-1],                     // S,E from left
                s_block[ii-1][j],                     f_block[ii-1][j],   // S,F from up
                x[global_i_idx], y[global_j_idx], score_fn
            );
            s_block[ii][j] = cell.s_val;
            e_block[ii][j] = cell.e_val;
            f_block[ii][j] = cell.f_val;

            // Pointers not stored globally here, only during final traceback recomputation

            if (cell.s_val > localBest.score) {
                localBest = {cell.s_val, global_i_idx + 1, j}; // Store 1-indexed matrix cell coords
            }
        }

        if (verbose && rank == 0) {

             // Progress bar only by rank 0 for combined progress

             if ((start_row_global + ii) % 100 == 0 || (start_row_global + ii) == m) {
                showProgressBar(start_row_global + ii, m);
             }
        }
    }

    // Send last computed S, E, F rows to the next rank
    if (rank + 1 < num_procs) {
        MPI_Send(s_block[local_num_rows].data(), n + 1, MPI_INT, rank + 1, 0, MPI_COMM_WORLD);
        MPI_Send(e_block[local_num_rows].data(), n + 1, MPI_INT, rank + 1, 1, MPI_COMM_WORLD);
        MPI_Send(f_block[local_num_rows].data(), n + 1, MPI_INT, rank + 1, 2, MPI_COMM_WORLD);
    }

    // Gather best hits (similar to before, localBest.i and localBest.j are 1-indexed matrix cell coords)
    int mine[4] = {localBest.score, rank, localBest.i, localBest.j};
    std::vector<int> all_best_hits_info;
    if (rank == 0) all_best_hits_info.resize(4 * num_procs);

    MPI_Gather(mine, 4, MPI_INT,
               rank == 0 ? all_best_hits_info.data() : nullptr,
               4, MPI_INT, 0, MPI_COMM_WORLD);

    int bestScore_global = 0, bestRank_global = 0, bestI_global = 0, bestJ_global = 0; // 1-indexed cell coords
    if (rank == 0) {
        for (int r_idx = 0; r_idx < num_procs; ++r_idx) {
            int s = all_best_hits_info[4 * r_idx];
            if (s > bestScore_global) {
                bestScore_global = s;
                bestRank_global = all_best_hits_info[4 * r_idx + 1];
                bestI_global = all_best_hits_info[4 * r_idx + 2]; // This is 1-indexed global row
                bestJ_global = all_best_hits_info[4 * r_idx + 3]; // This is 1-indexed global col
            }
        }
    }

    MPI_Bcast(&bestRank_global, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&bestScore_global, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&bestI_global, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&bestJ_global, 1, MPI_INT, 0, MPI_COMM_WORLD);

    std::string alignedX, alignedY;
    std::vector<std::pair<int,int>> local_path; // To store the path coordinates

    if (rank == bestRank_global) {

        // Recompute DP matrix and Pointers for traceback on bestRank_global

        // This will be up to (bestI_global, bestJ_global)
        // using the characters x[0...bestI_global-1] and y[0...bestJ_global-1]

        // Effectively, perform small_align_core logic here for the relevant sub-sequences.

        int tb_m = bestI_global; // Length of x subseq for traceback (1-indexed matrix size)
        int tb_n = bestJ_global; // Length of y subseq for traceback

        if (tb_m > 0 && tb_n > 0) { // Only proceed if there's a region to traceback
            std::vector<std::vector<int>> s_trace(tb_m + 1, std::vector<int>(tb_n + 1, 0));
            std::vector<std::vector<int>> e_trace(tb_m + 1, std::vector<int>(tb_n + 1, 0));
            std::vector<std::vector<int>> f_trace(tb_m + 1, std::vector<int>(tb_n + 1, 0));
            std::vector<std::vector<char>> p_trace(tb_m + 1, std::vector<char>(tb_n + 1, 'X'));

            // Initialize 0th row and col for S, E, F matrices to 0 for local affine

            // S[0][0] is 0. E[0][0], F[0][0] are effectively NEG_INF for transitions but stored as 0.

            for (int i = 1; i <= tb_m; ++i) {
                for (int j = 1; j <= tb_n; ++j) {
                    AffineDPScores cell = compute_local_affine_cell(
                        s_trace[i-1][j-1], e_trace[i-1][j-1], f_trace[i-1][j-1],
                        s_trace[i][j-1],   e_trace[i][j-1],
                        s_trace[i-1][j],   f_trace[i-1][j],
                        x[i-1], y[j-1], score_fn // Use original x, y characters
                    );
                    s_trace[i][j] = cell.s_val;
                    e_trace[i][j] = cell.e_val;
                    f_trace[i][j] = cell.f_val;
                    p_trace[i][j] = cell.ptr;
                }
            }

            // Actual traceback using p_trace from (bestI_global, bestJ_global)
            int cur_i = bestI_global;
            int cur_j = bestJ_global;

            local_path.emplace_back(cur_j, cur_i); // Add starting point of traceback

            while (cur_i > 0 || cur_j > 0) {
                if (s_trace[cur_i][cur_j] == 0 && p_trace[cur_i][cur_j] == 'X') break; // End of local alignment
                if (p_trace[cur_i][cur_j] == 'X' && (cur_i > 0 || cur_j > 0) ) break; // Safety break if X is hit prematurely

                char ptr_char = p_trace[cur_i][cur_j];
                if (ptr_char == 'M') {
                    alignedX.push_back(x[cur_i - 1]);
                    alignedY.push_back(y[cur_j - 1]);
                    --cur_i; --cur_j;
                } else if (ptr_char == 'F' || ptr_char == 'f') { // Gap in Y (target), x[cur_i-1] vs '-'
                    alignedX.push_back(x[cur_i - 1]);
                    alignedY.push_back('-');
                    --cur_i;
                } else if (ptr_char == 'E' || ptr_char == 'e') { // Gap in X (query), y[cur_j-1] vs '-'
                    alignedX.push_back('-');
                    alignedY.push_back(y[cur_j - 1]);
                    --cur_j;
                } else { // Reached 0 or an edge case not covered by M,E,F pointers
                    break;
                }
                local_path.emplace_back(cur_j, cur_i); // Add current point after moving
            }
            std::reverse(alignedX.begin(), alignedX.end());
            std::reverse(alignedY.begin(), alignedY.end());
            std::reverse(local_path.begin(), local_path.end()); // Path is built backwards
        }
    }

    // Timing and statistics (similar to before, using bestScore_global)
    auto t_end = Clock::now();
    auto time_ms = std::chrono::duration_cast<std::chrono::milliseconds>(t_end - t_start).count();

    // --- Gather full S-matrix (affine scores) to rank 0 ---
    // Each rank sends its s_block (excluding the received overlap row s_block[0])
    // This is similar to your original Gatherv logic for 'dp' matrix

    std::vector<std::vector<int>> fullDP_S_affine; // This will be the gathered S-matrix on Rank 0

    if (rank == 0) {
        fullDP_S_affine.assign(m + 1, std::vector<int>(n + 1, 0));

        // Rank 0 already has its first block (s_block[0]...s_block[local_num_rows])
        // Copy its own contribution (s_block[1] to s_block[local_num_rows])

        for (int i = 0; i < local_num_rows; ++i) {
             if (start_row_global + i < m) { // Ensure not to write out of bounds
                fullDP_S_affine[start_row_global + 1 + i] = s_block[i + 1];
             }
        }
    }

    // Each rank (except 0) sends its s_block[1]...s_block[local_num_rows] data to rank 0
    // Rank 0 receives them. This requires a more complex Gatherv or Send/Recv loop.
    // For simplicity, let's assume a Gatherv approach for the S-matrix part:

    // 1. Each rank prepares its S-block data to send (local_num_rows * (n+1) integers)

    std::vector<int> send_s_buffer;
    if (local_num_rows > 0) {
        send_s_buffer.resize(local_num_rows * (n + 1));
        for (int i = 0; i < local_num_rows; ++i) {
            std::memcpy(&send_s_buffer[i * (n + 1)], s_block[i + 1].data(), (n + 1) * sizeof(int));
        }
    }

    std::vector<int> recv_s_counts, s_displs;
    std::vector<int> toàn_s_recvbuf;

    if (rank == 0) {
        recv_s_counts.resize(num_procs);
        s_displs.resize(num_procs);
        int total_s_elements_to_recv = 0;
        for (int r = 0; r < num_procs; ++r) {
            int r_start = r * chunkSize;
            int r_end = std::min(r_start + chunkSize, m);
            int r_local_rows = r_end - r_start;
            recv_s_counts[r] = r_local_rows * (n + 1);
            s_displs[r] = (r == 0) ? 0 : s_displs[r-1] + recv_s_counts[r-1];
            total_s_elements_to_recv += recv_s_counts[r];
        }
        if (total_s_elements_to_recv > 0) {
             toàn_s_recvbuf.resize(total_s_elements_to_recv);
        }
    }

    MPI_Gatherv(
        send_s_buffer.data(),
        send_s_buffer.size(), MPI_INT,

        // send_s_buffer.size() is local_num_rows * (n+1)

        rank == 0 && !toàn_s_recvbuf.empty() ? toàn_s_recvbuf.data() : nullptr,
        rank == 0 ? recv_s_counts.data() : nullptr,
        rank == 0 ? s_displs.data() : nullptr,
        MPI_INT,
        0, MPI_COMM_WORLD
    );

    if (rank == 0) {

        // Reconstruct fullDP_S_affine from toàn_s_recvbuf
        // The 0-th row of fullDP_S_affine is already 0s.

        for (int r = 0; r < num_procs; ++r) {
            int r_start_global_row = r * chunkSize;

            // 0-indexed global start row

            int r_local_rows = recv_s_counts[r] / (n + 1);
            int offset_in_recvbuf = s_displs[r];
            for (int i = 0; i < r_local_rows; ++i) {
                if (r_start_global_row + i < m) {

                    // Ensure writing within bounds of m original rows

                    std::memcpy(
                        fullDP_S_affine[r_start_global_row + 1 + i].data(), // Write to global_row + 1
                        &toàn_s_recvbuf[offset_in_recvbuf + i * (n + 1)],
                        (n + 1) * sizeof(int)
                    );
                }
            }
        }
        // Now rank 0 has fullDP_S_affine (the S-matrix from the affine calculation)
        // Save fullDP_S_affine if binary or txt flags are set

        if (binary && !fullDP_S_affine.empty()) {
            writeDPMatrix(fullDP_S_affine, outdir + "/local_dp_matrix.bin");
        } else if (txt && !fullDP_S_affine.empty()) {
            writeRawDPMatrix(fullDP_S_affine, outdir + "/local_dp_matrix.txt");
        }
    }
    // END of Gathering full S-matrix

    // Rank 0 handles output (after potentially receiving alignedX, alignedY from bestRank_global)

    if (rank == 0) {
        if (bestRank_global != 0 && bestScore_global > 0) { /* ... MPI_Recv alignedX, alignedY ... */
            int len_ax, len_ay; MPI_Status status;
            MPI_Recv(&len_ax, 1, MPI_INT, bestRank_global, 10, MPI_COMM_WORLD, &status);
            alignedX.resize(len_ax); if (len_ax > 0) MPI_Recv(&alignedX[0], len_ax, MPI_CHAR, bestRank_global, 11, MPI_COMM_WORLD, &status);
            MPI_Recv(&len_ay, 1, MPI_INT, bestRank_global, 12, MPI_COMM_WORLD, &status);
            alignedY.resize(len_ay); if (len_ay > 0) MPI_Recv(&alignedY[0], len_ay, MPI_CHAR, bestRank_global, 13, MPI_COMM_WORLD, &status);
        }

        // Save local_path (which was populated on bestRank_global and sent to rank 0 if needed, or just written by bestRank_global)
        // For simplicity, let's assume bestRank_global writes its own path file if it computed one.
        // Or, if rank 0 is always bestRank_global (if not using MPI effectively or small data), it has local_path.
        // If bestRank_global is not 0, local_path needs to be sent to rank 0 or written by bestRank_global.
        // Let's assume rank 0 will write it if it has it (i.e., if rank 0 was bestRank_global)

        if (rank == bestRank_global && bestScore_global > 0 && !local_path.empty()) {
             std::ofstream pf_local(outdir + "/local_path.txt");
             if (pf_local) {
                for (const auto &p : local_path) {
                    pf_local << p.first << " " << p.second << "\n";
                }
                pf_local.close();
             } else {
                cerr << "Error: Cannot open " << outdir + "/local_path.txt" << " for writing.\n";
             }
        } else if (rank == 0 && bestRank_global != 0 && bestScore_global > 0) {
            // If rank 0 needs to write the path but wasn't the bestRank, it needs to receive local_path
            // This part is complex to add here; for now, assume bestRank handles its path file.
            // Or, only save path if rank 0 *is* bestRank.
        }
    }

    // Rank 0 does the printing and file writing

    if (rank == 0) {

        // (Broadcast alignedX and alignedY from bestRank_global to rank 0 if needed,
        // or have bestRank_global do the printing of alignment.
        // For simplicity, let's assume bestRank_global (if not 0) sends its alignedX, alignedY to rank 0)

        if (bestRank_global != 0 && bestScore_global > 0) {
            int len_ax = alignedX.size();
            int len_ay = alignedY.size();
            MPI_Send(&len_ax, 1, MPI_INT, 0, 10, MPI_COMM_WORLD);
            if (len_ax > 0) MPI_Send(alignedX.data(), len_ax, MPI_CHAR, 0, 11, MPI_COMM_WORLD);
            MPI_Send(&len_ay, 1, MPI_INT, 0, 12, MPI_COMM_WORLD);
            if (len_ay > 0) MPI_Send(alignedY.data(), len_ay, MPI_CHAR, 0, 13, MPI_COMM_WORLD);
        }
    } else if (rank == bestRank_global && bestScore_global > 0) {

        // bestRank sends its result to rank 0
        int len_ax = alignedX.size();
        int len_ay = alignedY.size();
        MPI_Send(&len_ax, 1, MPI_INT, 0, 10, MPI_COMM_WORLD);
        if (len_ax > 0) MPI_Send(alignedX.data(), len_ax, MPI_CHAR, 0, 11, MPI_COMM_WORLD);
        MPI_Send(&len_ay, 1, MPI_INT, 0, 12, MPI_COMM_WORLD);
        if (len_ay > 0) MPI_Send(alignedY.data(), len_ay, MPI_CHAR, 0, 13, MPI_COMM_WORLD);
    }

    if (rank == 0) {
        if (bestRank_global != 0 && bestScore_global > 0) {
            int len_ax, len_ay;
            MPI_Status status;
            MPI_Recv(&len_ax, 1, MPI_INT, bestRank_global, 10, MPI_COMM_WORLD, &status);
            alignedX.resize(len_ax);
            if (len_ax > 0) MPI_Recv(&alignedX[0], len_ax, MPI_CHAR, bestRank_global, 11, MPI_COMM_WORLD, &status);

            MPI_Recv(&len_ay, 1, MPI_INT, bestRank_global, 12, MPI_COMM_WORLD, &status);
            alignedY.resize(len_ay);
            if (len_ay > 0) MPI_Recv(&alignedY[0], len_ay, MPI_CHAR, bestRank_global, 13, MPI_COMM_WORLD, &status);
        }

        // Now rank 0 has bestScore_global, alignedX, alignedY
        size_t total = alignedX.size(), gaps = 0, matches = 0;
        if (total > 0) { // Calculate stats only if alignment exists
            for (size_t k = 0; k < total; ++k) {
                if (alignedX[k]=='-' || alignedY[k]=='-') ++gaps;
                else if (alignedX[k]==alignedY[k]) ++matches;
            }
        }
        double identity = (total > 0) ? (double(matches) / total) : 0.0;
        double coverage = (total > 0) ? (double(total - gaps) / total) : 0.0; // This is coverage of aligned region

        std::string acc1 = getAccession(header1, mode);
        std::string acc2 = getAccession(header2, mode);
        std::string gene1 = getGeneSymbol(header1, mode);
        std::string gene2 = getGeneSymbol(header2, mode);

        if (verbose) {
            std::cout << "\n\nLocal Alignment Score: " << bestScore_global << "\n"
                      << "Gap Open: " << GAP_OPEN << "\n"
                      << "Gap Extend: " << GAP_EXTEND << "\n"; // Now GAP_EXTEND is conceptually used
            if (total > 0) {
                 std::cout << "Matches: "  << matches << "\n"
                           << "Gaps:    "  << gaps    << "\n"
                           << "Total Aligned Length: " << total << "\n" // Renamed for clarity
                           << "Identity (of aligned region): " << identity*100.0 << "%\n"
                           << "Coverage (of aligned region): " << coverage*100.0 << "%\n";
            }
            std::cout << "Time:    "  << time_ms << " ms\n"
                      << "Query:   "  << acc1   << "\n"
                      << "Target:  "  << acc2   << "\n"
                      << "QueryID: "  << gene1  << "\n"
                      << "TargetID: " << gene2  << "\n";
            printColoredAlignment(alignedX, alignedY);
        }

        std::ofstream outf(outdir +"/local_alignment.fasta");
        if (outf) {
            savePlainAlignment(acc1, acc2, alignedX, alignedY, outf);
            outf.close();
        }
        // ... (JSON stats saving, similar to global, using bestScore_global etc.)
         std::ofstream js(outdir +"/local_stats.json");
        if (js) {
            js << std::fixed << std::setprecision(6)
               << "{\n"
               << "  \"method\":   \"local\",\n"
               << "  \"gap_open\": " << GAP_OPEN << ",\n"
               << "  \"gap_extend\": " << GAP_EXTEND << ",\n" // Include extend in report
               << "  \"score\":    " << bestScore_global << ",\n";
            if (total > 0) {
               js<< "  \"matches\":  " << matches   << ",\n"
                 << "  \"gaps\":     " << gaps      << ",\n"
                 << "  \"aligned_length\":    " << total     << ",\n"
                 << "  \"identity\": " << identity  << ",\n"
                 << "  \"coverage_aligned\": " << coverage  << ",\n";
            }
               js<< "  \"time_ms\":  " << time_ms   << ",\n"
                 << "  \"query\":    \"" << acc1  << "\",\n"
                 << "  \"target\":   \"" << acc2  << "\",\n"
                 << "  \"queryid\":  \"" << gene1 << "\",\n"
                 << "  \"targetid\": \"" << gene2 << "\",\n"
                 << "  \"query_length_original\": " << m << ",\n" // Original lengths
                 << "  \"target_length_original\": " << n << "\n"
                 << "}\n";
            js.close();
        }
    }
    // Ensure all ranks reach this point before finalizing
    MPI_Barrier(MPI_COMM_WORLD);
}

/**
* @brief Compute the Longest Common Subsequence (LCS) of two sequences.
*
* This function computes the LCS using dynamic programming and saves the results.
*
* @param x         The first sequence (string).
* @param y         The second sequence (string).
* @param header1   Header for the first sequence.
* @param header2   Header for the second sequence.
* @param outdir    Output directory for results.
* @param mode      Scoring mode (MODE_DNA or MODE_PROTEIN).
**/
void lcs(const string &x, const string &y,
         const string &header1, const string &header2,
         const std::string &outdir, ScoreMode mode) { // ScoreMode might not be relevant for LCS scoring itself
    const int m = x.size(), n = y.size();
    vector<vector<int>> L(m + 1, vector<int>(n + 1, 0)); // DP table for lengths
    vector<vector<char>> b(m + 1, vector<char>(n + 1, ' ')); // Traceback pointers
    vector<pair<int, int>> lcs_path_coords; // To store (j,i) matrix coordinates of the LCS path

    // --- Step 1: Fill DP table (L) and pointer table (b) ---
    // The 'b' matrix stores 'D' (diagonal), 'U' (up), 'L' (left).
    for (int i = 1; i <= m; ++i) {
        for (int j = 1; j <= n; ++j) {
            if (x[i - 1] == y[j - 1]) {
                L[i][j] = L[i - 1][j - 1] + 1;
                b[i][j] = 'D'; // Diagonal
            } else if (L[i - 1][j] >= L[i][j - 1]) { // Tie-breaking: prefer up over left
                L[i][j] = L[i - 1][j];
                b[i][j] = 'U'; // Up
            } else {
                L[i][j] = L[i][j - 1];
                b[i][j] = 'L'; // Left
            }
        }
        if (verbose) { // Simplified progress for LCS DP fill
            if (i % 100 == 0 || i == m) showProgressBar(i, m);
        }
    }
    if (verbose && m > 0) std::cout << std::endl; // Newline after progress bar

    // --- Step 2: Traceback to reconstruct LCS string AND the gapped alignment ---
    string lcs_str = "";
    string alignedX = "";
    string alignedY = "";
    int cur_i = m, cur_j = n;

    while (cur_i > 0 && cur_j > 0) {
        lcs_path_coords.emplace_back(cur_j, cur_i); // Store current cell coordinates
        if (b[cur_i][cur_j] == 'D') {
            lcs_str += x[cur_i - 1];
            alignedX += x[cur_i - 1];
            alignedY += y[cur_j - 1];
            cur_i--;
            cur_j--;
        } else if (b[cur_i][cur_j] == 'U') {
            alignedX += x[cur_i - 1];
            alignedY += '-';
            cur_i--;
        } else { // Must be 'L'
            alignedX += '-';
            alignedY += y[cur_j - 1];
            cur_j--;
        }
    }

    // Handle remaining characters if one sequence is longer
    while (cur_i > 0) {
        lcs_path_coords.emplace_back(cur_j, cur_i);
        alignedX += x[cur_i - 1];
        alignedY += '-';
        cur_i--;
    }
    while (cur_j > 0) {
        lcs_path_coords.emplace_back(cur_j, cur_i);
        alignedX += '-';
        alignedY += y[cur_j - 1];
        cur_j--;
    }
    lcs_path_coords.emplace_back(0,0); // Add origin

    reverse(lcs_str.begin(), lcs_str.end());
    reverse(alignedX.begin(), alignedX.end());
    reverse(alignedY.begin(), alignedY.end());
    reverse(lcs_path_coords.begin(), lcs_path_coords.end());

    // --- Step 3: Saving outputs ---
    // Save the LCS string itself
    std::ofstream outfile_lcs(outdir + "/lcs.fasta");
    if (!outfile_lcs) {
        cerr << "Error: Cannot open " << outdir + "/lcs.fasta for writing\n";
    } else {
        saveLCS(getAccession(header1, mode) + "_" + getAccession(header2, mode),
                lcs_str,
                outfile_lcs);
        outfile_lcs.close();
    }

    // Save the alignment showing the LCS
    std::ofstream outfile_alignment(outdir + "/lcs_alignment.fasta");
     if (!outfile_alignment) {
        cerr << "Error: Cannot open " << outdir + "/lcs_alignment.fasta for writing\n";
    } else {
        savePlainAlignment(getAccession(header1, mode) + "_LCS_aligned",
                           getAccession(header2, mode) + "_LCS_aligned",
                           alignedX, alignedY,
                           outfile_alignment);
        outfile_alignment.close();
    }

    // Save lcs_path_coords (matrix cell coordinates of the path)
    {
        ofstream pf(outdir + "/lcs_path.txt");
        if (!pf) {
            cerr << "Error: Cannot open " << outdir + "/lcs_path.txt" << " for writing.\n";
        } else {
            for (const auto &p : lcs_path_coords) {
                pf << p.first << " " << p.second << "\n";
            }
            pf.close();
        }
    }

    // Save traceback pointer matrix 'b' if flags are set
    if (binary) {
        writeCharMatrix(b, outdir + "/lcs_traceback_pointers.bin");
    } else if (txt) {
        writeRawCharMatrix(b, outdir + "/lcs_traceback_pointers.txt");
    }

    // Save DP length matrix 'L' if flags are set
    if (binary) {
        writeDPMatrix(L, outdir + "/lcs_dp_lengths.bin");
    } else if (txt) {
        writeRawDPMatrix(L, outdir + "/lcs_dp_lengths.txt");
    }


    // --- Step 4: Displaying the alignment ---
    if (verbose) {

        std::cout << "LCS Length: " << lcs_str.size() << "\n";
        std::cout << "\n\nLCS Based Alignment:\n";

        // This is true LCS length (no explicit gap penalties, just +1 for match).
        // Display the alignment constructed
        printColoredAlignment(alignedX, alignedY);

        std::cout << "\nLongest Common Subsequence String:\n";
        for (size_t k = 0; k < lcs_str.size(); k += LINE_WIDTH) {
            std::cout << lcs_str.substr(k, LINE_WIDTH) << "\n";
        }
        std::cout << "\n";
    }
}

/**
 * @brief Main function to run the MPI-based sequence alignment tool.
 *
 * This function initializes MPI, parses command line arguments, reads input files,
 * and performs sequence alignments based on user choices.
 *
 * @param argc Number of command line arguments.
 * @param argv Array of command line arguments.
 * @return Exit status code.
 */
int main(int argc, char** argv) {
    MPI_Init(&argc, &argv);
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    try {
        std::string file1, file2, outdir = ".";
        int choice = -1;
        ScoreMode mode = MODE_DNA;

        // Parse arguments
        for (int i = 1; i < argc; ++i) {
            std::string arg = argv[i];
            if (arg == "--query" && i + 1 < argc) {
                file1 = argv[++i];
            } else if (arg == "--target" && i + 1 < argc) {
                file2 = argv[++i];
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
        if (rank == 0) {
            processFasta(file1, header1, seq1);
            processFasta(file2, header2, seq2);
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
