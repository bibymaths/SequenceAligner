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

#include <iostream>
#include <chrono>
#include <iomanip>
#include <fstream>
#include <vector>
#include <bitset>
#include <algorithm>
#include <stdexcept>
#include <omp.h>
#include <immintrin.h>
#include <mpi.h>
#include <cstdint>
#include <sstream>
#include "EDNAFULL.h"
#include "EBLOSUM62.h"
#include <unordered_map>
#include <filesystem>
#include <climits>

bool verbose = false;
using namespace std;
enum ScoreMode { MODE_DNA, MODE_PROTEIN };

/// Gap penalty
static const int GAP      = -2;    // penalty for a gap in local alignment
static const int GAP_OPEN   = -1;  // penalty to open a gap in global alignment
static const double GAP_EXTEND = -0.5;  // penalty to extend an existing gap in global alignment
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


// The EDNAFULL matrix is defined as a 2D C array of size N×N:
static const std::unordered_map<char,int> edna_index = {
  {'A',0},{'C',1},{'G',2},{'T',3},
  {'U',3}, // T=U
  {'R',4},{'Y',5},{'S',6},{'W',7},
  {'K',8},{'M',9},{'B',10},{'D',11},
  {'H',12},{'V',13},{'N',14},{'X',14}
};

// The BLOSUM62 matrix is defined as a 2D C array of size 24×24:
static const std::unordered_map<char, int> blosum62_index = {
  {'A', 0}, {'R', 1}, {'N', 2}, {'D', 3}, {'C', 4}, {'Q', 5}, {'E', 6}, {'G', 7},
  {'H', 8}, {'I', 9}, {'L',10}, {'K',11}, {'M',12}, {'F',13}, {'P',14}, {'S',15},
  {'T',16}, {'W',17}, {'Y',18}, {'V',19}, {'B',20}, {'Z',21}, {'X',22}, {'*',23}
};


/**
 * @brief Lookup the EDNAFULL score between two bases (incl. ambiguous codes).
 * @param x First base (IUPAC code).
 * @param y Second base (IUPAC code).
**/
inline int edna_score(char x, char y) {
    auto ix = edna_index.find(x);
    auto iy = edna_index.find(y);
    if (ix==edna_index.end() || iy==edna_index.end())
        throw std::runtime_error(std::string("Invalid base: ")+x+","+y);
    return static_cast<int>(EDNAFULL_matrix[ix->second][iy->second]);
}

/**
 * @brief Lookup the BLOSUM62 score between two amino acids.
 * @param x First amino acid (1-letter code).
 * @param y Second amino acid (1-letter code).
**/
inline int blosum62_score(char x, char y) {
    auto ix = blosum62_index.find(x);
    auto iy = blosum62_index.find(y);
    if (ix == blosum62_index.end() || iy == blosum62_index.end())
        throw std::runtime_error(std::string("Invalid protein code: ") + x + "," + y);
    return static_cast<int>(EBLOSUM62_matrix[ix->second][iy->second]);
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
 * @brief Dump an already-aligned pair of sequences (no colors) into a stream,
 *        in the same LINE_WIDTH blocks as printColoredAlignment.
 *
 * @param seq1 First aligned sequence (may contain ‘-’).
 * @param seq2 Second aligned sequence (may contain ‘-’).
 * @param os   Output stream to use (e.g. a std::ofstream).
 */
void savePlainAlignment(const string &seq1,
                        const string &seq2,
                        ostream &os)
{
    size_t len1 = seq1.size();
    size_t len2 = seq2.size();
    size_t maxLength = max(len1, len2);

    string a1 = seq1, a2 = seq2;
    if (len1 < maxLength) a1.append(maxLength - len1, '-');
    if (len2 < maxLength) a2.append(maxLength - len2, '-');

    for (size_t i = 0; i < maxLength; i += LINE_WIDTH) {
        size_t end = min(i + LINE_WIDTH, maxLength);
        os << "\nPosition: " << (i+1) << " - " << end << "\n";
        os << a1.substr(i, end - i) << "\n";
        os << a2.substr(i, end - i) << "\n";
    }
    os << "\n";
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
                        ScoreMode mode,
                        vector<int>& prev_row,
                        vector<int>& prev_gapX,
                        vector<int>& prev_gapY,
                        vector<int>& curr_row,
                        vector<int>& curr_gapX,
                        vector<int>& curr_gapY)
{
    (void)prev_gapY;
    int n = y.size();
    curr_row .assign(n+1, INT_MIN/2);
    curr_gapX.assign(n+1, INT_MIN/2);
    curr_gapY.assign(n+1, INT_MIN/2);

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
        int mscore = score(x[i-1], y[j-1], mode);
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
                 const std::string &outdir, ScoreMode mode) {
    int m = x.size(), n = y.size();

    vector<int> prev_row, prev_gapX, prev_gapY;
    initAffineDP(n, prev_row, prev_gapX, prev_gapY, true);
    vector<int> curr_row, curr_gapX, curr_gapY;

    vector<char> prev_trace(n + 1, '0');
    vector<char> curr_trace(n + 1, '0');

    using Clock = std::chrono::high_resolution_clock;
    auto t_start = Clock::now();

     for (int i = 1; i <= m; ++i) {
          computeAffineDPRow(i, x, y, mode,
                             prev_row, prev_gapX, prev_gapY,
                             curr_row, curr_gapX, curr_gapY);
          prev_row.swap(curr_row);
          prev_gapX.swap(curr_gapX);
          prev_gapY.swap(curr_gapY);

        if (verbose) { if (i % 1000 == 0 || i == m) {
            showProgressBar(i, m);
        }}
    }

    // Traceback
    string alignedX, alignedY;
    int i = m, j = n;

    while (i > 0 || j > 0) {
        if (i == 0) {
            alignedX += '-';
            alignedY += y[--j];
            continue;
        }
        if (j == 0) {
            alignedX += x[--i];
            alignedY += '-';
            continue;
        }

        int matchScore = score(x[i - 1], y[j - 1], mode);
        int diag = prev_row[j - 1] + matchScore;
        int up = prev_row[j] + GAP;
        int left = curr_row[j - 1] + GAP;

        int score = max({diag, up, left});
        curr_row[j] = score;

        char move;
        if (score == diag) move = 'd';
        else if (score == up) move = 'u';
        else move = 'l';

        if (move == 'd') {
            alignedX += x[i - 1];
            alignedY += y[j - 1];
            i--; j--;
        } else if (move == 'u') {
            alignedX += x[i - 1];
            alignedY += '-';
            i--;
        } else {
            alignedX += '-';
            alignedY += y[j - 1];
            j--;
        }
        swap(prev_row, curr_row);
    }

    auto t_end = Clock::now();
    auto time_ms = std::chrono::duration_cast<std::chrono::milliseconds>(t_end - t_start).count();

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

    std::string modeDir = (mode == MODE_DNA ? "dna" : "protein");

    if (verbose) {
        cout << "\n\nGlobal Alignment Score: " << prev_row[n] << "\n";
        cout << "Matches: " << matches << "\n";
        cout << "Gaps:    " << gaps << "\n";
        cout << "Total:   " << total << "\n";
        cout << "Identity: " << identity * 100.0f << "%\n";
        cout << "Coverage: " << coverage * 100.0f << "%\n";
        cout << "Time:    " << time_ms << " ms\n";
        cout << "Query:   " << accession1 << "\n";
        cout << "Target:  " << accession2 << "\n";
        cout << "QueryID:  " << gene1 << "\n";
        cout << "TargetID:  " << gene2 << "\n";
        printColoredAlignment(alignedX, alignedY);
    }

    std::ofstream outfile(outdir + "/" + modeDir + "/global_alignment.txt");
    if (outfile) {
        outfile << "\nSequence 1: " << header1;
        outfile << "\nSequence 2: " << header2;
        outfile << "\n\nGlobal Alignment Score: " << prev_row[n] << "\n\n";
        savePlainAlignment(alignedX, alignedY, outfile);
        outfile.close();
    } else {
        cerr << "Error: Unable to open output file global_alignment.txt\n";
    }

    ofstream js(outdir + "/" + modeDir + "/global_stats.json");
    if (js) {
      js << fixed << setprecision(6)
         << "{\n"
         << "  \"method\":      \"global\",\n"
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
         << "  \"targetid\":       \"" << gene2 << "\"\n"
         << "}\n";
      js.close();
    } else {
      cerr << "Error: cannot open global_stats.json\n";
    }

}

/**
 * @brief Perform a parallel, MPI‐distributed local (Smith–Waterman) alignment.
 *
 * Splits X across ranks, does a two‐row DP on each chunk (with OpenMP + SIMD),
 * finds the local max, gathers to pick the global max, then traceback on the
 * winning rank.
 *
 * @param x Full first sequence.
 * @param y Full second sequence.
 */
void localalign(const std::string &x, const std::string &y,
                const std::string &header1, const std::string &header2,
                const std::string &outdir, ScoreMode mode) {
    int m = x.size(), n = y.size();
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    int chunkSize = (m + size - 1) / size;
    int start    = rank * chunkSize;
    int end      = std::min(start + chunkSize, m);
    int localRows = end - start;

    // two-row DP per rank
    std::vector<int> prev(n+1, 0), curr(n+1, 0);
    if (rank > 0) {
        // receive the “row 0” from rank-1
        MPI_Recv(prev.data(), n+1, MPI_INT, rank-1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }

    using Clock = std::chrono::high_resolution_clock;
    auto t_start = Clock::now();

    struct Loc { int score, i, j; };
    Loc localBest{0,0,0};

//    #pragma omp parallel
    {
        Loc thrBest{0,0,0};
        #pragma omp for
        for (int ii = 1; ii <= localRows; ++ii) {
            int gi = start + ii - 1;
            curr[0] = 0;
            #pragma omp simd
            for (int j = 1; j <= n; ++j) {
                int ms = score(x[gi], y[j - 1], mode);
                int v  = prev[j-1] + ms;
                int u  = prev[j]   + GAP;
                int l  = curr[j-1] + GAP;
                int s  = v;
                if (u > s) s = u;
                if (l > s) s = l;
                if (s < 0) s = 0;
                curr[j] = s;
                if (s > thrBest.score) thrBest = { s, gi, j };
            }
            if (verbose) {
                if (ii % 1000 == 0 || ii == localRows) {
                    showProgressBar(ii, localRows);
                }
            }
            std::swap(prev, curr);
        }
        #pragma omp critical
        if (thrBest.score > localBest.score)
            localBest = thrBest;
    }

    // pass our final row down
    if (rank + 1 < size) {
        MPI_Send(prev.data(), n+1, MPI_INT, rank+1, 0, MPI_COMM_WORLD);
    }

    // gather (score, rank, I, J) on rank 0
    int mine[4] = { localBest.score, rank, localBest.i, localBest.j };
    std::vector<int> all;
    if (rank == 0) all.resize(4*size);

    MPI_Gather(
      mine, 4, MPI_INT,
      rank==0 ? all.data() : nullptr, 4, MPI_INT,
      0, MPI_COMM_WORLD
    );

    int bestScore=0, bestRank=0, bestI=0, bestJ=0;
    if (rank == 0) {
      for (int r = 0; r < size; ++r) {
        int s = all[4*r+0];
        if (s > bestScore) {
          bestScore = s;
          bestRank  = all[4*r+1];
          bestI     = all[4*r+2];
          bestJ     = all[4*r+3];
        }
      }
    }
    // broadcast winner info
    MPI_Bcast(&bestRank, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&bestScore, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&bestI,     1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&bestJ,     1, MPI_INT, 0, MPI_COMM_WORLD);

    // traceback on the winning rank
    std::string alignedX, alignedY;
    if (rank == bestRank) {
        int I = bestI, J = bestJ;
        std::vector<std::vector<int>> dp(I+1, std::vector<int>(J+1,0));
        for (int i = 1; i <= I; ++i) {
            for (int j = 1; j <= J; ++j) {
                int ms = score(x[i - 1], y[j - 1], mode);
                int v  = dp[i-1][j-1] + ms;
                int u  = dp[i-1][j]   + GAP;
                int l  = dp[i][j-1]   + GAP;
                int s  = v;
                if (u > s) s = u;
                if (l > s) s = l;
                if (s < 0) s = 0;
                dp[i][j] = s;
            }
        }
        int i = I, j = J;
        while (i > 0 && j > 0 && dp[i][j] > 0) {
            int cur = dp[i][j];
            int ms  = score(x[i - 1], y[j - 1], mode);
            if (cur == dp[i-1][j-1] + ms) {
                alignedX.push_back(x[i-1]);
                alignedY.push_back(y[j-1]);
                --i; --j;
            }
            else if (cur == dp[i-1][j] + GAP) {
                alignedX.push_back(x[i-1]);
                alignedY.push_back('-');
                --i;
            }
            else {
                alignedX.push_back('-');
                alignedY.push_back(y[j-1]);
                --j;
            }
        }
        std::reverse(alignedX.begin(), alignedX.end());
        std::reverse(alignedY.begin(), alignedY.end());
    }

    // ship back to rank 0
    if (rank == bestRank && rank != 0) {
        int a = alignedX.size(), b = alignedY.size();
        MPI_Send(&a, 1, MPI_INT, 0, 1, MPI_COMM_WORLD);
        MPI_Send(alignedX.data(), a, MPI_CHAR, 0, 2, MPI_COMM_WORLD);
        MPI_Send(&b, 1, MPI_INT, 0, 3, MPI_COMM_WORLD);
        MPI_Send(alignedY.data(), b, MPI_CHAR, 0, 4, MPI_COMM_WORLD);
    }
    else if (rank == 0) {
        if (bestRank != 0) {
            int a,b;
            MPI_Recv(&a, 1, MPI_INT, bestRank, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            alignedX.resize(a);
            MPI_Recv(&alignedX[0], a, MPI_CHAR, bestRank, 2, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            MPI_Recv(&b, 1, MPI_INT, bestRank, 3, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            alignedY.resize(b);
            MPI_Recv(&alignedY[0], b, MPI_CHAR, bestRank, 4, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }

        auto t_end = Clock::now();
        auto time_ms = std::chrono::duration_cast<std::chrono::milliseconds>(t_end - t_start).count();

        size_t total = alignedX.size();
        size_t gaps  = 0, matches = 0;
        for (size_t i = 0; i < total; ++i) {
            if (alignedX[i] == '-' || alignedY[i] == '-')
            ++gaps;
            else if (alignedX[i] == alignedY[i])
            ++matches;
        }
        double identity = double(matches) / double(total);
        double coverage = double(total - gaps) / double(total);

        std::string accession1 = getAccession(header1, mode);
        std::string accession2 = getAccession(header2, mode);
        std::string gene1 = getGeneSymbol(header1, mode);
        std::string gene2 = getGeneSymbol(header2, mode);
        std::string modeDir = (mode == MODE_DNA ? "dna" : "protein");

        if (verbose) {
            std::cout << "\n\nLocal Alignment Score: " << bestScore << "\n";
            std::cout << "Matches: " << matches << "\n";
            std::cout << "Gaps:    " << gaps << "\n";
            std::cout << "Total:   " << total << "\n";
            std::cout << "Identity: " << identity * 100.0f << "%\n";
            std::cout << "Coverage: " << coverage * 100.0f << "%\n";
            cout << "Time:    " << time_ms << " ms\n";
            cout << "Query:   " << accession1 << "\n";
            cout << "Target:  " << accession2 << "\n";
            cout << "QueryID:  " << gene1 << "\n";
            cout << "TargetID:  " << gene2 << "\n";
            printColoredAlignment(alignedX, alignedY);
        }
        std::ofstream outfile(outdir + "/" + modeDir + "/local_alignment.txt");
        if (outfile) {
            outfile << "\nSequence 1: " << header1;
            outfile << "\nSequence 2: " << header2;
            outfile << "\n\nLocal Alignment Score: " << bestScore << "\n\n";
            savePlainAlignment(alignedX, alignedY, outfile);
            outfile.close();
        } else {
            cerr << "Error: Unable to open output file local_alignment.txt\n";
        }

        ofstream js(outdir + "/" + modeDir + "/local_stats.json");
        if (js) {
          js << fixed << setprecision(6)
             << "{\n"
             << "  \"method\":      \"local\",\n"
             << "  \"score\":       " << bestScore << ",\n"
             << "  \"matches\":     " << matches << ",\n"
             << "  \"gaps\":        " << gaps << ",\n"
             << "  \"total\":       " << total << ",\n"
             << "  \"identity\":    " << identity << ",\n"
             << "  \"coverage\":    " << coverage << ",\n"
             << "  \"time_ms\":     " << time_ms << ",\n"
             << "  \"query\":       \"" << accession1 << "\",\n"
             << "  \"target\":      \"" << accession2 << "\",\n"
             << "  \"queryid\":       \"" << gene1 << "\",\n"
             << "  \"targetid\":       \"" << gene2 << "\"\n"
             << "}\n";
          js.close();
        } else {
          cerr << "Error: cannot open local_stats.json\n";
        }

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
        if (b[i][j] == 'D')       { lcs_str += x[i-1]; --i; --j; }
        else if (b[i][j] == 'U')  { --i; }
        else                      { --j; }
    }
    reverse(lcs_str.begin(), lcs_str.end());
    std::string modeDir = (mode == MODE_DNA ? "dna" : "protein");
    std::ofstream outfile(outdir + "/" + modeDir + "/lcs.txt");
    if (verbose) {
        std::cout << "\n\nLCS length: " << prev[n] << "\n\nLCS: \n";
    }
    outfile << "\nSequence 1: " << header1
            << "\nSequence 2: " << header2;
    outfile << "\n\nLCS length: " << prev[n];
    outfile << "\n\nLCS: \n";
    for (size_t i = 0; i < lcs_str.length(); i += LINE_WIDTH) {
        if (verbose) {
            std::cout << lcs_str.substr(i, LINE_WIDTH) << "\n";
        }
        outfile << lcs_str.substr(i, LINE_WIDTH) << "\n";
    }
    outfile.close();
}


/**
 * @brief Program entry point: read arguments, load FASTA, dispatch chosen method.
 *
 * Usage: `aligner <fasta1> <fasta2> <method>`
 *  - method=1  → global alignment
 *  - method=2  → local alignment
 *  - method=3  → LCS
 *  - --outdir <output_directory> (optional)
 *  - --mode <dna|protein> (optional, default: dna)
 * @param argc Argument count.
 * @param argv Argument vector.
 * @return     0 on success, nonzero on failure.
 */
int main(int argc, char** argv) {

    MPI_Init(&argc, &argv);
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    try {
        if (argc < 4) {
            if (rank == 0) {
                std::cout <<
                    "Usage: " << argv[0] << " <fasta_file1> <fasta_file2> <method>\n"
                    "  method: 1=global, 2=local, 3=LCS\n"
                    "Options:\n"
                    "  --mode <dna|protein>     Scoring mode (default: dna)\n"
                    "  --outdir <directory>     Output directory (default: .)\n"
                    "  --verbose                Print progress and debug info\n"
                    "  --help                   Show this help message\n";
            }
            MPI_Finalize();
            return 1;
        }

        // Now it's safe to use argv[1]..argv[3]
        std::string file1 = argv[1];
        std::string file2 = argv[2];
        int choice = std::atoi(argv[3]);

        // Optional args
        std::string outdir = ".";
        ScoreMode mode = MODE_DNA;

        for (int i = 4; i < argc; ++i) {
            std::string arg = argv[i];
            if (arg == "--mode" && i + 1 < argc) {
                std::string val = argv[++i];
                if (val == "dna") mode = MODE_DNA;
                else if (val == "protein") mode = MODE_PROTEIN;
                else {
                    if (rank == 0) std::cerr << "Unknown mode: " << val << "\n";
                    MPI_Finalize();
                    return 1;
                }
            } else if (arg == "--outdir" && i + 1 < argc) {
                outdir = argv[++i];
            } else if (arg == "--verbose") {
                verbose = true;
            } else {
                if (rank == 0) std::cerr << "Unknown option: " << arg << "\n";
                MPI_Finalize();
                return 1;
            }
        }


        // Create output directory if needed
        if (rank == 0) {
            std::string modeDir = (mode == MODE_DNA ? "dna" : "protein");
            std::filesystem::create_directories(outdir + "/" + modeDir);
        }

        // Read FASTA and broadcast
        std::string seq1, seq2, header1, header2;
        if (rank == 0) {
            processFasta(file1, header1, seq1);
            processFasta(file2, header2, seq2);
        }
        if (verbose && rank == 0) {
            std::cout << "Sequence 1 Header: " << header1 << "\n";
            std::cout << "Sequence 1 Length: " << seq1.size() << " bases\n";
            std::cout << "Sequence 2 Header: " << header2 << "\n";
            std::cout << "Sequence 2 Length: " << seq2.size() << " bases\n";
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
        MPI_Bcast(&choice, 1, MPI_INT, 0, MPI_COMM_WORLD);

        // Dispatch
        switch (choice) {
            case 1:
                if (rank == 0)
                    globalalign(seq1, seq2, header1, header2, outdir, mode);
                break;
            case 2:
                localalign(seq1, seq2, header1, header2, outdir, mode);
                break;
            case 3:
                if (rank == 0)
                    lcs(seq1, seq2, header1, header2, outdir, mode);
                break;
            default:
                if (rank == 0)
                    std::cerr << "Invalid method. Use 1=global, 2=local, 3=LCS.\n";
        }

    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << "\n";
    }

    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Finalize();
    return 0;
}
