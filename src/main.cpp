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
#include <unordered_map>
using namespace std;

/**
 * @file aligner.cpp
 * @brief Sequence aligner supporting global, local and LCS methods using MPI, OpenMP, and AVX2.
 */

constexpr int EDNA_SIZE = 15;
/// Score awarded for a character match
static const int MATCH    =  1;
/// Score penalty for a character mismatch
static const int MISMATCH = -1;
/// Gap penalty
static const int GAP      = -1;
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
static constexpr int N = EDNA_SIZE;
static const double (*EDNA_mat)[N] = EDNAFULL_matrix;

// Map each IUPAC base code to its row/col in EDNAFULL
static const std::unordered_map<char,int> edna_index = {
  {'A',0},{'C',1},{'G',2},{'T',3},
  {'U',3}, // T=U
  {'R',4},{'Y',5},{'S',6},{'W',7},
  {'K',8},{'M',9},{'B',10},{'D',11},
  {'H',12},{'V',13},{'N',14},{'X',14}
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
    return static_cast<int>(EDNA_mat[ix->second][iy->second]);
}

/**
 * @brief Display a textual progress bar on stdout.
 * @param progress Number of units completed so far.
 * @param total    Total number of units to complete.
 */
void showProgressBar(int progress, int total) {
    int barWidth = 50;
    float percentage = (float)progress / total;
    int pos = barWidth * percentage;

    cout << "[";
    for (int i = 0; i < barWidth; ++i) {
        if (i < pos) cout << "=";
        else if (i == pos) cout << ">";
        else cout << " ";
    }
    cout << "] " << int(percentage * 100.0) << "%\r";
    cout.flush();
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

/**
 * @brief Perform a global (Needleman–Wunsch) alignment of two sequences.
 *
 * Uses linear-space DP (two-row) with simple traceback to reconstruct the
 * full alignment, then prints and saves it.
 *
 * @param x First sequence.
 * @param y Second sequence.
 */
void globalalign(const string &x, const string &y) {
    int m = x.size(), n = y.size();

    vector<int> prev_row(n + 1, 0);
    vector<int> curr_row(n + 1, 0);
    vector<char> prev_trace(n + 1, '0');
    vector<char> curr_trace(n + 1, '0');

    using Clock = std::chrono::high_resolution_clock;
    auto t_start = Clock::now();

    for (int j = 0; j <= n; ++j) {
        prev_row[j] = j * GAP;
        prev_trace[j] = 'l';
    }

    for (int i = 1; i <= m; i++) {
        curr_row[0] = i * GAP;
        curr_trace[0] = 'u';

        for (int j = 1; j <= n; j++) {
//            int matchScore = (x[i - 1] == y[j - 1]) ? MATCH : MISMATCH;
            int matchScore = edna_score(x[i-1], y[j-1]);
            int diag = prev_row[j - 1] + matchScore;
            int up = prev_row[j] + GAP;
            int left = curr_row[j - 1] + GAP;

            int score = max({diag, up, left});
            curr_row[j] = score;

            if (score == diag) curr_trace[j] = 'd';
            else if (score == up) curr_trace[j] = 'u';
            else curr_trace[j] = 'l';
        }

        prev_row.swap(curr_row);
        prev_trace.swap(curr_trace);

        if (i % 1000 == 0 || i == m) {
            showProgressBar(i, m);
        }
    }

    // Traceback
    string alignedX, alignedY;
    int i = m, j = n;
    vector<char> *curr_trace_ptr = &prev_trace, *prev_trace_ptr = &curr_trace;

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

//        int matchScore = (x[i - 1] == y[j - 1]) ? MATCH : MISMATCH;
        int matchScore = edna_score(x[i-1], y[j-1]);

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

    reverse(alignedX.begin(), alignedX.end());
    reverse(alignedY.begin(), alignedY.end());

    cout << "\nGlobal Alignment Score: " << prev_row[n] << endl;
    printColoredAlignment(alignedX, alignedY);

    ofstream outfile("global_alignment.txt");
    if (outfile) {
        outfile << "Global Alignment Score: " << prev_row[n] << "\n\n";
        printColoredAlignment(alignedX, alignedY);
        savePlainAlignment(alignedX, alignedY, outfile);
        outfile.close();
    } else {
        cerr << "Error: Unable to open output file global_alignment.txt\n";
    }

    ofstream js("global_stats.json");
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
         << "  \"time_ms\":     " << time_ms << "\n"
         << "}\n";
      js.close();
    } else {
      cerr << "Error: cannot open global_stats.json\n";
    }

}


// --- Advanced MPI Local Alignment with Recursive Multi-Chunk Traceback ---
// Tags for traceback messages.
const int TRACEBACK_REQ_TAG = 100;
const int TRACEBACK_RESP_TAG = 101;


/**
 * @brief Holds the result of a chunked traceback step in MPI local alignment.
 */
struct TracebackResult {
    string alignedX;   ///< Partial aligned segment from sequence X
    string alignedY;   ///< Partial aligned segment from sequence Y
    int    new_global_i; ///< Global row index where traceback stopped
    int    new_j;        ///< Column index where traceback stopped
};

/**
 * @brief Compute a 1D index for a 2D DP array of size (m+1)×(n+1).
 * @param i Row index (0..m)
 * @param j Column index (0..n)
 * @param n Number of columns in the sequences
 * @return  Flattened index = i*(n+1) + j
 */
inline int idx(int i, int j, int n) {
    return i * (n + 1) + j;
}

/**
 * @brief Recursively traceback through a local DP chunk in MPI local alignment.
 *
 * If the traceback reaches the top of this chunk and still has positive score,
 * it will send a request to the previous MPI rank to continue.
 *
 * @param global_i   Current global row index in X.
 * @param j          Current column index in Y.
 * @param rank       MPI rank of this process.
 * @param dpStart    Row index in the local DP array that corresponds to global start.
 * @param chunkStart Global row index of the first row in this chunk.
 * @param dp         Flattened DP matrix for this chunk.
 * @param trace      Flattened traceback directions for this chunk.
 * @param n          Number of columns in Y.
 * @param x          The full X sequence.
 * @param y          The full Y sequence.
 * @return           A TracebackResult containing aligned fragments and new indices.
 */
TracebackResult recursiveTraceback(int global_i, int j, int rank, int dpStart, int chunkStart,
                                    const vector<int>& dp, const vector<char>& trace, int n,
                                    const string &x, const string &y)
{
    TracebackResult result;
    result.alignedX = "";
    result.alignedY = "";

    // Convert global row to local DP index.
    int local_i = global_i - chunkStart + dpStart;

    // Walk the DP/traceback in this chunk.
    while (local_i > dpStart && j > 0 && dp[idx(local_i, j, n)] > 0) {
        char dir = trace[idx(local_i, j, n)];
        if (dir == 'd') {
            result.alignedX.push_back(x[global_i]);
            result.alignedY.push_back(y[j - 1]);
            local_i--; j--; global_i--;
        } else if (dir == 'u') {
            result.alignedX.push_back(x[global_i]);
            result.alignedY.push_back('-');
            local_i--; global_i--;
        } else if (dir == 'l') {
            result.alignedX.push_back('-');
            result.alignedY.push_back(y[j - 1]);
            j--;
        } else {
            break;
        }
    }
    reverse(result.alignedX.begin(), result.alignedX.end());
    reverse(result.alignedY.begin(), result.alignedY.end());
    result.new_global_i = global_i;
    result.new_j = j;

    // If reached the top boundary and still nonzero, request continuation.
    if (local_i == dpStart && dp[idx(local_i, j, n)] > 0 && rank > 0) {
        int req[2] = { global_i, j };
        MPI_Send(req, 2, MPI_INT, rank - 1, TRACEBACK_REQ_TAG, MPI_COMM_WORLD);

        int lenX, lenY, new_global_i, new_j;
        MPI_Recv(&lenX, 1, MPI_INT, rank - 1, TRACEBACK_RESP_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        char *bufX = new char[lenX + 1];
        MPI_Recv(bufX, lenX, MPI_CHAR, rank - 1, TRACEBACK_RESP_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        bufX[lenX] = '\0';
        MPI_Recv(&lenY, 1, MPI_INT, rank - 1, TRACEBACK_RESP_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        char *bufY = new char[lenY + 1];
        MPI_Recv(bufY, lenY, MPI_CHAR, rank - 1, TRACEBACK_RESP_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        bufY[lenY] = '\0';
        MPI_Recv(&new_global_i, 1, MPI_INT, rank - 1, TRACEBACK_RESP_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Recv(&new_j, 1, MPI_INT, rank - 1, TRACEBACK_RESP_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        string segX(bufX), segY(bufY);
        delete[] bufX; delete[] bufY;

        result.alignedX = segX + result.alignedX;
        result.alignedY = segY + result.alignedY;
        result.new_global_i = new_global_i;
        result.new_j = new_j;
    }
    return result;
}

/**
 * @brief Service an incoming MPI traceback request from a later rank.
 *
 * Receives the (i,j) position, performs recursiveTraceback, and
 * sends back the aligned segments and updated indices.
 *
 * @param chunkStart Global index of the first row in this chunk.
 * @param dpStart    Row index in the local DP array that corresponds to chunkStart.
 * @param dp         Flattened DP matrix for this chunk.
 * @param trace      Flattened traceback directions for this chunk.
 * @param n          Number of columns in Y.
 * @param x          The full X sequence.
 * @param y          The full Y sequence.
 * @param rank       MPI rank of this process.
 */
void handleTracebackRequest(int chunkStart, int dpStart,
                            const vector<int>& dp, const vector<char>& trace, int n,
                            const string &x, const string &y, int rank)
{
    int req[2];
    MPI_Status status;
    MPI_Recv(req, 2, MPI_INT, MPI_ANY_SOURCE, TRACEBACK_REQ_TAG, MPI_COMM_WORLD, &status);
    int requester = status.MPI_SOURCE;
    int global_i = req[0];
    int j = req[1];
    TracebackResult res = recursiveTraceback(global_i, j, rank, dpStart, chunkStart, dp, trace, n, x, y);
    int lenX = res.alignedX.size(), lenY = res.alignedY.size();
    MPI_Send(&lenX, 1, MPI_INT, requester, TRACEBACK_RESP_TAG, MPI_COMM_WORLD);
    MPI_Send(res.alignedX.c_str(), lenX, MPI_CHAR, requester, TRACEBACK_RESP_TAG, MPI_COMM_WORLD);
    MPI_Send(&lenY, 1, MPI_INT, requester, TRACEBACK_RESP_TAG, MPI_COMM_WORLD);
    MPI_Send(res.alignedY.c_str(), lenY, MPI_CHAR, requester, TRACEBACK_RESP_TAG, MPI_COMM_WORLD);
    MPI_Send(&res.new_global_i, 1, MPI_INT, requester, TRACEBACK_RESP_TAG, MPI_COMM_WORLD);
    MPI_Send(&res.new_j, 1, MPI_INT, requester, TRACEBACK_RESP_TAG, MPI_COMM_WORLD);
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
void localalign(const std::string &x, const std::string &y) {
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
            int gi = start + ii - 1; // global index in x
            curr[0] = 0;
            #pragma omp simd
            for (int j = 1; j <= n; ++j) {
//                int ms = (x[gi] == y[j-1] ? MATCH : MISMATCH);
                int ms = edna_score(x[gi], y[j-1]);
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
            if (rank == 0 && (ii % 1000 == 0 || ii == localRows))
                showProgressBar(ii, localRows);
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
//                int ms = (x[i-1]==y[j-1]?MATCH:MISMATCH);
                int ms = edna_score(x[i-1], y[j-1]);
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
//            int ms  = (x[i-1]==y[j-1]?MATCH:MISMATCH);
            int ms  = edna_score(x[i-1], y[j-1]);
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

        std::cout << "\nLocal Alignment Score: " << bestScore << "\n";
        printColoredAlignment(alignedX, alignedY);
        ofstream outfile("local_alignment.txt");
        if (outfile) {
            outfile << "Local Alignment Score: " << bestScore << "\n\n";
            savePlainAlignment(alignedX, alignedY, outfile);
            outfile.close();
        } else {
            cerr << "Error: Unable to open output file local_alignment.txt\n";
        }

        ofstream js("local_stats.json");
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
             << "  \"time_ms\":     " << time_ms << "\n"
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
void lcs(const string &x, const string &y) {
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
                if (i % 1000 == 0 || i == m) showProgressBar(i, m);
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

    std::ofstream outfile ("lcs.txt");
    cout << "\n\nLCS length: " << prev[n]
         << "\n\nLCS: "        << lcs_str << "\n";
    outfile << "\nLCS length: " << prev[n]
         << "\n\nLCS: "        << lcs_str << "\n";
    outfile.close();
}


/**
 * @brief Program entry point: read arguments, load FASTA, dispatch chosen method.
 *
 * Usage: `aligner <fasta1> <fasta2> <method>`
 *  - method=1  → global alignment
 *  - method=2  → local alignment
 *  - method=3  → LCS
 *
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
            if (rank == 0)
                cout << "Usage: " << argv[0] << " <fasta_file1> <fasta_file2> <method: 1=global, 2=local, 3=LCS>" << endl;
            MPI_Finalize();
            return 1;
        }


        string seq1, seq2, header1, header2;
        if (rank == 0) {
            processFasta(argv[1], header1, seq1);
            processFasta(argv[2], header2, seq2);
            cout << "Sequence 1 Header: " << header1 << "\n";
            cout << "Sequence 1 Length: " << seq1.size() << " bases\n";
            cout << "Sequence 2 Header: " << header2 << "\n";
            cout << "Sequence 2 Length: " << seq2.size() << " bases\n";
        }

        int seq1_len = seq1.size(), seq2_len = seq2.size();
        MPI_Bcast(&seq1_len, 1, MPI_INT, 0, MPI_COMM_WORLD);
        MPI_Bcast(&seq2_len, 1, MPI_INT, 0, MPI_COMM_WORLD);

        if (rank != 0) {
            seq1.resize(seq1_len);
            seq2.resize(seq2_len);
        }

        MPI_Bcast(&seq1[0], seq1_len, MPI_CHAR, 0, MPI_COMM_WORLD);
        MPI_Bcast(&seq2[0], seq2_len, MPI_CHAR, 0, MPI_COMM_WORLD);

        int choice = atoi(argv[3]);
        MPI_Bcast(&choice, 1, MPI_INT, 0, MPI_COMM_WORLD);

        switch (choice) {
            case 1:
                if (rank == 0)
                    globalalign(seq1, seq2);
                break;
            case 2:
                localalign(seq1, seq2);
                break;
            case 3:
                lcs(seq1, seq2);
                break;
            default:
                if (rank == 0)
                    cout << "Invalid choice! Use 1 for global, 2 for local." << endl;
        }

    } catch (const exception &e) {
        cerr << e.what() << endl;
    }


    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Finalize();
    return 0;
}