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

void simdGlobalAlign(const string &x, const string &y) {
    int m = x.size(), n = y.size();
    vector<vector<int>> dp(m + 1, vector<int>(n + 1, 0));
    vector<vector<char>> traceback(m + 1, vector<char>(n + 1, ' '));

    // Initialize first row and column (gap penalties)
    for (int i = 0; i <= m; i++) dp[i][0] = i * GAP;
    for (int j = 0; j <= n; j++) dp[0][j] = j * GAP;

    // Fill DP table using AVX2 vectorized operations
    for (int i = 1; i <= m; i++) {
        #pragma omp parallel for
        for (int j = 1; j <= n; j += 8) {
            __m256i v_left = _mm256_loadu_si256((__m256i*)&dp[i][j - 1]);
            __m256i v_up = _mm256_loadu_si256((__m256i*)&dp[i - 1][j]);
            __m256i v_diag = _mm256_loadu_si256((__m256i*)&dp[i - 1][j - 1]);

            __m256i v_match = _mm256_set1_epi32(x[i - 1] == y[j - 1] ? MATCH : MISMATCH);
            __m256i v_gap = _mm256_set1_epi32(GAP);

            __m256i v_matchResult = _mm256_add_epi32(v_diag, v_match);
            __m256i v_gapLeft = _mm256_add_epi32(v_left, v_gap);
            __m256i v_gapUp = _mm256_add_epi32(v_up, v_gap);

            __m256i v_result = _mm256_max_epi32(v_matchResult, v_gapLeft);
            v_result = _mm256_max_epi32(v_result, v_gapUp);

            _mm256_storeu_si256((__m256i*)&dp[i][j], v_result);
        }

        // Update progress bar every 1000 rows or on the last row.
        if (i % 1000 == 0 || i == m) {
            showProgressBar(i, m);
        }
    }

    // Traceback to reconstruct alignment.
    string alignedX, alignedY;
    int i = m, j = n;
    while (i > 0 || j > 0) {
        if (i > 0 && j > 0 && dp[i][j] == dp[i - 1][j - 1] + (x[i - 1] == y[j - 1] ? MATCH : MISMATCH)) {
            alignedX += x[i - 1];
            alignedY += y[j - 1];
            i--; j--;
        } else if (i > 0 && dp[i][j] == dp[i - 1][j] + GAP) {
            alignedX += x[i - 1];
            alignedY += '-';
            i--;
        } else {
            alignedX += '-';
            alignedY += y[j - 1];
            j--;
        }
    }
    reverse(alignedX.begin(), alignedX.end());
    reverse(alignedY.begin(), alignedY.end());

    cout << "\nGlobal Alignment Score: " << dp[m][n] << endl;
    printColoredAlignment(alignedX, alignedY);

    // Save colored alignment to a file.
    ofstream outfile("alignment.txt");
    if (!outfile) {
        cerr << "Error: Unable to open output file alignment.txt" << endl;
    } else {
        outfile << "\nGlobal Alignment Score: " << dp[m][n] << "\n\n";
        printColoredAlignment(alignedX, alignedY, outfile);
        outfile.close();
        cout << "\nAlignment saved to alignment.txt" << endl;
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
void mpiLocalAlign(const string &x, const string &y) {
    int m = x.size(), n = y.size();
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    const int overlap = 500; // Adjust based on expected alignment length

    int chunkSize = (m + size - 1) / size;
    int chunkStart = rank * chunkSize;
    int chunkEnd = min(chunkStart + chunkSize, m);
    int localRows = chunkEnd - chunkStart;

    int effectiveRows = localRows;
    if (rank > 0) effectiveRows += overlap;

    vector<int> dp((effectiveRows + 1) * (n + 1), 0);
    vector<char> trace((effectiveRows + 1) * (n + 1), '0');

    auto local_idx = [n](int i, int j) -> int {
        return i * (n + 1) + j;
    };

    if (rank > 0) {
        int rowsToRecv = min(overlap, chunkSize);
        MPI_Recv(&dp[local_idx(0, 0)], (rowsToRecv + 1) * (n + 1), MPI_INT, rank - 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        for (int i = 0; i <= rowsToRecv; i++) {
            for (int j = 0; j <= n; j++) {
                trace[local_idx(i, j)] = 'x';
            }
        }
    } else {
        for (int j = 0; j <= n; j++) {
            dp[local_idx(0, j)] = 0;
            trace[local_idx(0, j)] = '0';
        }
    }

    int dpStart = (rank > 0) ? overlap + 1 : 1;
    int dpRows = dpStart - 1 + localRows;

    int local_maxScore = 0, local_maxRow = 0, local_maxCol = 0;
    for (int i = dpStart; i <= dpRows; i++) {
        int global_i = chunkStart + i - dpStart;
        dp[local_idx(i, 0)] = 0;
        trace[local_idx(i, 0)] = '0';
        for (int j = 1; j <= n; j++) {
            int matchScore = (x[global_i] == y[j - 1]) ? MATCH : MISMATCH;
            int score_diag = dp[local_idx(i - 1, j - 1)] + matchScore;
            int score_up   = dp[local_idx(i - 1, j)] + GAP;
            int score_left = dp[local_idx(i, j - 1)] + GAP;
            int score = 0;
            char direction = '0';
            if (score_diag > 0 || score_up > 0 || score_left > 0) {
                score = score_diag;
                direction = 'd';
                if (score_up > score)   { score = score_up;   direction = 'u'; }
                if (score_left > score) { score = score_left; direction = 'l'; }
                if (score < 0) { score = 0; direction = '0'; }
            }
            dp[local_idx(i, j)] = score;
            trace[local_idx(i, j)] = direction;
            if (score > local_maxScore) {
                local_maxScore = score;
                local_maxRow = i;
                local_maxCol = j;
            }
        }
        // (Optional) update progress bar on rank 0.
    }

    if (rank < size - 1) {
        MPI_Send(&dp[local_idx(dpRows, 0)], n + 1, MPI_INT, rank + 1, 0, MPI_COMM_WORLD);
    }

    int local_data[3] = { local_maxScore, chunkStart + local_maxRow - dpStart, local_maxCol };
    int *gathered = nullptr;
    if (rank == 0)
        gathered = new int[3 * size];
    MPI_Gather(local_data, 3, MPI_INT, gathered, 3, MPI_INT, 0, MPI_COMM_WORLD);

    int global_maxScore = 0, winnerRank = 0, winner_global_i = 0, winner_j = 0;
    if (rank == 0) {
        for (int r = 0; r < size; r++) {
            int score = gathered[3 * r];
            if (score > global_maxScore) {
                global_maxScore = score;
                winnerRank = r;
                winner_global_i = gathered[3 * r + 1];
                winner_j = gathered[3 * r + 2];
            }
        }
        delete[] gathered;
    }
    MPI_Bcast(&winnerRank, 1, MPI_INT, 0, MPI_COMM_WORLD);

    string alignedX, alignedY;
    if (rank == winnerRank) {
        TracebackResult res = recursiveTraceback(winner_global_i, winner_j, rank, dpStart, chunkStart, dp, trace, n, x, y);
        alignedX = res.alignedX;
        alignedY = res.alignedY;
    }

    double startTime = MPI_Wtime();
    while (MPI_Wtime() - startTime < 1.0) {
        int flag = 0;
        MPI_Iprobe(MPI_ANY_SOURCE, TRACEBACK_REQ_TAG, MPI_COMM_WORLD, &flag, MPI_STATUS_IGNORE);
        if (flag) {
            handleTracebackRequest(chunkStart, dpStart, dp, trace, n, x, y, rank);
        }
    }

    if (rank == winnerRank && rank != 0) {
        int lenX = alignedX.size(), lenY = alignedY.size();
        MPI_Send(&lenX, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
        MPI_Send(alignedX.c_str(), lenX, MPI_CHAR, 0, 0, MPI_COMM_WORLD);
        MPI_Send(&lenY, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
        MPI_Send(alignedY.c_str(), lenY, MPI_CHAR, 0, 0, MPI_COMM_WORLD);
    } else if (rank == 0 && winnerRank != 0) {
        int lenX, lenY;
        MPI_Recv(&lenX, 1, MPI_INT, winnerRank, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        char *bufX = new char[lenX + 1];
        MPI_Recv(bufX, lenX, MPI_CHAR, winnerRank, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        bufX[lenX] = '\0';
        MPI_Recv(&lenY, 1, MPI_INT, winnerRank, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        char *bufY = new char[lenY + 1];
        MPI_Recv(bufY, lenY, MPI_CHAR, winnerRank, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        bufY[lenY] = '\0';
        alignedX = string(bufX);
        alignedY = string(bufY);
        delete[] bufX; delete[] bufY;
    }

    if (rank == 0) {
        // Save colored alignment to a file.
        ofstream outfile("alignment.txt");
        if (!outfile) {
            cerr << "Error: Unable to open output file alignment.txt" << endl;
        }
        outfile << "Local Alignment Score: " << global_maxScore << "\n\n";
        // Optionally, you can also output FASTA header information here if needed.
        outfile << "Aligned X: " << alignedX << "\n\n";
        outfile << "Aligned Y: " << alignedY << "\n\n";
        // Print the colored alignment (with ANSI codes) into the file.
        printColoredAlignment(alignedX, alignedY, outfile);
        outfile.close();

        cout << "\nLocal Alignment Score: " << global_maxScore << endl;
        if (alignedX.empty() || alignedY.empty()) {
            cout << "Warning: Alignment traceback failed or spans too many chunks.\n";
        } else {
            printColoredAlignment(alignedX, alignedY);
            cout << "\nAlignment saved to alignment.txt" << endl;
        }
    }
}

int main(int argc, char** argv) {
    MPI_Init(&argc, &argv);

    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    try {
        if (argc < 3) {
            if (rank == 0)
                cout << "Usage: " << argv[0] << " <fasta_file1> <fasta_file2>" << endl;
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

        int choice = -1;
        if (rank == 0) {
            cout << "Select Alignment Method:\n1. Global\n2. Local \nChoice: ";
            cin >> choice;
        }
        MPI_Bcast(&choice, 1, MPI_INT, 0, MPI_COMM_WORLD);

        switch (choice) {
            case 1:
                if (rank == 0)
                    simdGlobalAlign(seq1, seq2); // Only rank 0 runs global alignment
                break;
            case 2:
                // All ranks participate in local alignment.
                mpiLocalAlign(seq1, seq2);
                break;
            default:
                if (rank == 0)
                    cout << "Invalid choice!" << endl;
        }
    } catch (const exception &e) {
        cerr << e.what() << endl;
    }

    MPI_Finalize();
    return 0;
}

