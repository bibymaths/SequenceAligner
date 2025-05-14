#include <iostream>
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

using namespace std;

// Constants
const int MATCH = 1;
const int MISMATCH = -1;
const int GAP = -1;
const int LINE_WIDTH = 150;

// ANSI color codes
#define RESET "\033[0m"
#define GREEN "\033[32m"
#define RED "\033[31m"
#define CYAN "\033[36m"

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

// Process a FASTA file, extracting the first header (without '>') and the sequence.
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
                header = line.substr(1); // remove '>' from header
                headerSet = true;
            }
            continue; // skip header lines thereafter
        }
        sequence += line;
    }
}

// Overloaded printColoredAlignment accepts an output stream (default = cout)
void printColoredAlignment(const string &seq1, const string &seq2, ostream &os = cout) {
    size_t length1 = seq1.size();
    size_t length2 = seq2.size();
    size_t maxLength = max(length1, length2);  // Ensure full length is displayed

    // Pad shorter sequence with gaps ('-') to match length
    string aligned1 = seq1;
    string aligned2 = seq2;
    if (length1 < maxLength) aligned1.append(maxLength - length1, '-');
    if (length2 < maxLength) aligned2.append(maxLength - length2, '-');

    // Display alignment in chunks of LINE_WIDTH
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

void globalalign(const string &x, const string &y) {
    int m = x.size(), n = y.size();

    vector<int> prev_row(n + 1, 0);
    vector<int> curr_row(n + 1, 0);
    vector<char> prev_trace(n + 1, '0');
    vector<char> curr_trace(n + 1, '0');

    for (int j = 0; j <= n; ++j) {
        prev_row[j] = j * GAP;
        prev_trace[j] = 'l';
    }

    for (int i = 1; i <= m; i++) {
        curr_row[0] = i * GAP;
        curr_trace[0] = 'u';

        for (int j = 1; j <= n; j++) {
            int matchScore = (x[i - 1] == y[j - 1]) ? MATCH : MISMATCH;

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

    // Recompute rows during traceback
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

        int matchScore = (x[i - 1] == y[j - 1]) ? MATCH : MISMATCH;

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

        // recompute rows for the next traceback step
        swap(prev_row, curr_row);
    }

    reverse(alignedX.begin(), alignedX.end());
    reverse(alignedY.begin(), alignedY.end());

    cout << "\nGlobal Alignment Score: " << prev_row[n] << endl;
    printColoredAlignment(alignedX, alignedY);

    ofstream outfile("alignment.txt");
    if (outfile) {
        outfile << "Global Alignment Score: " << prev_row[n] << "\n\n";
        printColoredAlignment(alignedX, alignedY, outfile);
        outfile.close();
        cout << "\nAlignment saved to alignment.txt\n";
    } else {
        cerr << "Error: Unable to open output file alignment.txt\n";
    }
}


// --- Advanced MPI Local Alignment with Recursive Multi-Chunk Traceback ---
// Tags for traceback messages.
const int TRACEBACK_REQ_TAG = 100;
const int TRACEBACK_RESP_TAG = 101;

// Structure to hold traceback results.
struct TracebackResult {
    string alignedX;
    string alignedY;
    int new_global_i;  // Global row where traceback ended in this chunk.
    int new_j;         // Column where traceback ended.
};

// Helper to compute 2D index into 1D array.
inline int idx(int i, int j, int n) {
    return i * (n + 1) + j;
}

/*
  Recursive traceback function.
  Performs traceback in the local DP matrix. If the traceback reaches the top boundary
  (local row == dpStart) with a nonzero cell and the process is not rank 0, it sends a
  traceback request to the previous process.
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

/*
  Handle an incoming traceback request.
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

/*
  Advanced MPI Local Alignment with Recursive Multi-Chunk Traceback.
  (The DP and traceback arrays are computed with an extra overlap region.)
*/
void localalign(const string &x, const string &y) {
    int m = x.size(), n = y.size();
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    // split rows of x across ranks
    int chunkSize = (m + size - 1) / size;
    int start    = rank * chunkSize;
    int end      = min(start + chunkSize, m);
    int localRows = end - start;

    // allocate local DP table (rows 0..localRows) × (cols 0..n)
    vector<vector<int>> dp(localRows+1, vector<int>(n+1, 0));

    // receive the “row 0” from previous rank if any
    if (rank > 0) {
        MPI_Recv(dp[0].data(), n+1, MPI_INT, rank-1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }

    // fill local table and track the best cell
    int local_maxScore = 0, local_maxI = 0, local_maxJ = 0;
    for (int i = 1; i <= localRows; ++i) {
        int gi = start + i - 1;              // global index in x
        for (int j = 1; j <= n; ++j) {
            int ms = (x[gi] == y[j-1]) ? MATCH : MISMATCH;
            dp[i][j] = max({ 0,
                             dp[i-1][j-1] + ms,
                             dp[i-1][j]   + GAP,
                             dp[i][j-1]   + GAP });
            if (dp[i][j] > local_maxScore) {
                local_maxScore = dp[i][j];
                local_maxI     = gi;   // store global row
                local_maxJ     = j;
            }
        }
    }

    // pass our last row down to the next rank
    if (rank+1 < size) {
        MPI_Send(dp[localRows].data(), n+1, MPI_INT, rank+1, 0, MPI_COMM_WORLD);
    }

    // gather (score, rank, I, J) on rank 0
    struct Loc { int score, rank, I, J; };
    Loc mine{ local_maxScore, rank, local_maxI, local_maxJ };
    vector<Loc> all;
    if (rank == 0) all.resize(size);
    MPI_Gather(&mine, sizeof(Loc)/sizeof(int), MPI_INT,
               all.data(), sizeof(Loc)/sizeof(int), MPI_INT,
               0, MPI_COMM_WORLD);

    // pick the true best on rank 0
    int bestScore=0, bestRank=0, bestI=0, bestJ=0;
    if (rank == 0) {
        for (auto &c : all) {
            if (c.score > bestScore) {
                bestScore = c.score;
                bestRank  = c.rank;
                bestI     = c.I;
                bestJ     = c.J;
            }
        }
    }

    // broadcast the winner info to everyone
    MPI_Bcast(&bestScore, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&bestRank,  1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&bestI,     1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&bestJ,     1, MPI_INT, 0, MPI_COMM_WORLD);

    // now only the winning rank does the traceback
    string alignedX, alignedY;
    if (rank == bestRank) {
        // map global bestI back into our local table row index
        int local_i = bestI - start + 1;
        int j       = bestJ;
        while (local_i > 0 && j > 0 && dp[local_i][j] > 0) {
            int cur = dp[local_i][j];
            int gi  = start + local_i - 1;
            int ms  = (x[gi] == y[j-1]) ? MATCH : MISMATCH;
            if (local_i>0 && j>0 && cur == dp[local_i-1][j-1] + ms) {
                alignedX.push_back(x[gi]);
                alignedY.push_back(y[j-1]);
                --local_i; --j;
            } else if (local_i>0 && cur == dp[local_i-1][j] + GAP) {
                alignedX.push_back(x[gi]);
                alignedY.push_back('-');
                --local_i;
            } else {
                alignedX.push_back('-');
                alignedY.push_back(y[j-1]);
                --j;
            }
        }
        reverse(alignedX.begin(), alignedX.end());
        reverse(alignedY.begin(), alignedY.end());
    }

    // send the two strings back to rank 0 if needed
    if (rank == bestRank && rank != 0) {
        int lenX = alignedX.size(), lenY = alignedY.size();
        MPI_Send(&lenX, 1, MPI_INT, 0, 1, MPI_COMM_WORLD);
        MPI_Send(alignedX.data(), lenX, MPI_CHAR, 0, 2, MPI_COMM_WORLD);
        MPI_Send(&lenY, 1, MPI_INT, 0, 3, MPI_COMM_WORLD);
        MPI_Send(alignedY.data(), lenY, MPI_CHAR, 0, 4, MPI_COMM_WORLD);
    } else if (rank == 0) {
        // rank 0 may need to receive from bestRank
        if (bestRank != 0) {
            int lenX, lenY;
            MPI_Recv(&lenX, 1, MPI_INT, bestRank, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            alignedX.resize(lenX);
            MPI_Recv(alignedX.data(), lenX, MPI_CHAR, bestRank, 2, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            MPI_Recv(&lenY, 1, MPI_INT, bestRank, 3, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            alignedY.resize(lenY);
            MPI_Recv(alignedY.data(), lenY, MPI_CHAR, bestRank, 4, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }
        // finally, print
        cout << "\nLocal Alignment Score: " << bestScore << "\n";
        printColoredAlignment(alignedX, alignedY);
    }
}


void lcs(const string &x, const string &y) {
    const int m = x.size(), n = y.size();
    vector<int> prev(n + 1, 0), curr(n + 1, 0);
    vector<vector<char>> b(m + 1, vector<char>(n + 1, ' '));

    const int blockSize = 32;  // For AVX2 (256 bits = 32 bytes)
    const int j_limit = n - (n % blockSize);

    #pragma omp parallel
    {
        vector<int> tmp(blockSize);

        for (int i = 1; i <= m; ++i) {
            #pragma omp for schedule(static)
            for (int j = 1; j <= j_limit; j += blockSize) {
                __m256i vx = _mm256_set1_epi8(x[i - 1]);
                __m256i vy = _mm256_loadu_si256((__m256i const*)(y.data() + j - 1));
                __m256i eq = _mm256_cmpeq_epi8(vx, vy);
                _mm256_storeu_si256((__m256i*)tmp.data(), eq);

                for (int k = 0; k < blockSize; ++k) {
                    int jj = j + k;
                    if (jj > n) break;

                    if (tmp[k] == -1) {  // true from AVX2 comparison (all bits set)
                        curr[jj] = prev[jj - 1] + 1;
                        b[i][jj] = 'D';
                    } else if (prev[jj] >= curr[jj - 1]) {
                        curr[jj] = prev[jj];
                        b[i][jj] = 'U';
                    } else {
                        curr[jj] = curr[jj - 1];
                        b[i][jj] = 'L';
                    }
                }
            }

            // Fallback for remainder (non-AVX2-aligned)
            #pragma omp for schedule(static)
            for (int j = j_limit + 1; j <= n; ++j) {
                if (x[i - 1] == y[j - 1]) {
                    curr[j] = prev[j - 1] + 1;
                    b[i][j] = 'D';
                } else if (prev[j] >= curr[j - 1]) {
                    curr[j] = prev[j];
                    b[i][j] = 'U';
                } else {
                    curr[j] = curr[j - 1];
                    b[i][j] = 'L';
                }
            }

            #pragma omp single
            {
                swap(prev, curr);
                if (i % 1000 == 0 || i == m)
                    showProgressBar(i, m);
            }
        }
    }

    // Traceback
    int i = m, j = n;
    string lcs_str;
    while (i > 0 && j > 0) {
        if (b[i][j] == 'D') {
            lcs_str += x[i - 1];
            --i; --j;
        } else if (b[i][j] == 'U') {
            --i;
        } else {
            --j;
        }
    }
    reverse(lcs_str.begin(), lcs_str.end());
    cout << "\n\nLCS length: " << prev[n] << "\n\nLCS: " << lcs_str << endl;
}

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

    MPI_Finalize();
    return 0;
}


