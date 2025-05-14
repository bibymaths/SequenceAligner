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

void processFasta(const string &filename, string &sequence) {
    ifstream file(filename);
    if (!file) throw runtime_error("Error: Unable to open " + filename);
    string line;
    sequence.clear();
    while (getline(file, line)) {
        if (line.empty() || line[0] == '>') continue;
        sequence += line;
    }
}

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

void printColoredAlignment(const string &seq1, const string &seq2) {
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
        cout << "\nPosition: " << i + 1 << " - " << end << "\n";

        for (size_t j = i; j < end; j++) {
            if (aligned1[j] == aligned2[j]) cout << GREEN << aligned1[j] << RESET;
            else if (aligned1[j] == '-' || aligned2[j] == '-') cout << RED << aligned1[j] << RESET;
            else cout << CYAN << aligned1[j] << RESET;
        }
        cout << "\n";

        for (size_t j = i; j < end; j++) {
            if (aligned1[j] == aligned2[j]) cout << GREEN << aligned2[j] << RESET;
            else if (aligned1[j] == '-' || aligned2[j] == '-') cout << RED << aligned2[j] << RESET;
            else cout << CYAN << aligned2[j] << RESET;
        }
        cout << "\n";
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

        // Insert progress bar update every 1000 rows or at the last row
        if (i % 1000 == 0 || i == m) {
            showProgressBar(i, m);
        }
    }

    // Traceback to reconstruct alignment
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
}

void mpiLocalAlign(const string &x, const string &y) {
    int m = x.size(), n = y.size();

    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    int numThreads = omp_get_max_threads();

    // Divide the rows of x among ranks.
    int chunkSize = (m + size - 1) / size;
    int start = rank * chunkSize;
    int end = min(start + chunkSize, m);
    int localRows = end - start;

    // Allocate DP table with (localRows+1) rows.
    // For rank 0, dp[0] is all zeros.
    // For rank > 0, dp[0] will be received from rank-1.
    vector<vector<int>> dp(localRows + 1, vector<int>(n + 1, 0));

    // For ranks > 0, receive the last row from the previous process
    if (rank > 0) {
        MPI_Recv(dp[0].data(), n + 1, MPI_INT, rank - 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }

    int local_maxScore = 0, local_maxI = 0, local_maxJ = 0;

    // Fill DP table for rows 1 to localRows.
    for (int i = 1; i <= localRows; i++) {
        int global_i = start + i - 1; // Map local row to global index in x.
        for (int j = 1; j <= n; j++) {
            int matchScore = (x[global_i] == y[j - 1]) ? MATCH : MISMATCH;
            dp[i][j] = max({0, dp[i - 1][j - 1] + matchScore,
                              dp[i - 1][j] + GAP,
                              dp[i][j - 1] + GAP});
            if (dp[i][j] > local_maxScore) {
                local_maxScore = dp[i][j];
                local_maxI = global_i;  // store global index
                local_maxJ = j;
            }
        }
        // Optionally, update progress bar here.
    }

    // For ranks that are not the last, send the last computed row to the next rank.
    if (rank < size - 1) {
        MPI_Send(dp[localRows].data(), n + 1, MPI_INT, rank + 1, 0, MPI_COMM_WORLD);
    }

    // Now, use MPI_Reduce (or MPI_Gather) to find the overall best score and corresponding indices.
    int local_data[3] = {local_maxScore, local_maxI, local_maxJ};
    int global_data[3] = {0, 0, 0};
    MPI_Reduce(local_data, global_data, 3, MPI_INT, MPI_MAX, 0, MPI_COMM_WORLD);

    // The following traceback is simplified and works correctly only if the optimal
    // alignment is entirely contained within one process’s computed DP table.
    // For a cross-chunk alignment, a full traceback across processes requires a more
    // elaborate mechanism (for example, a wavefront approach).
    if (rank == 0) {
        cout << "\nLocal Alignment Score: " << global_data[0] << endl;
        // If the winning score is from rank 0, perform traceback using dp.
        // Otherwise, you’d need to receive the alignment from the winning rank.
        // Here we assume for simplicity that the best alignment is in rank 0.
        int gi = global_data[1]; // global index in x where max was found
        int gj = global_data[2];
        // For rank 0, start == 0, so local index equals global index + 1.
        int li = gi - start + 1;
        string alignedX, alignedY;
        while (li > 0 && gj > 0 && dp[li][gj] > 0) {
            int current = dp[li][gj];
            int diag = (li - 1 >= 0 && gj - 1 >= 0) ? dp[li - 1][gj - 1] : -1000000;
            int up   = (li - 1 >= 0) ? dp[li - 1][gj] : -1000000;
            int left = (gj - 1 >= 0) ? dp[li][gj - 1] : -1000000;
            int matchScore = (x[gi - 1] == y[gj - 1]) ? MATCH : MISMATCH;
            if (current == diag + matchScore) {
                alignedX.push_back(x[gi - 1]);
                alignedY.push_back(y[gj - 1]);
                gi--; li--; gj--;
            } else if (current == up + GAP) {
                alignedX.push_back(x[gi - 1]);
                alignedY.push_back('-');
                gi--; li--;
            } else if (current == left + GAP) {
                alignedX.push_back('-');
                alignedY.push_back(y[gj - 1]);
                gj--;
            } else {
                break;
            }
        }
        reverse(alignedX.begin(), alignedX.end());
        reverse(alignedY.begin(), alignedY.end());
        printColoredAlignment(alignedX, alignedY);
    }
}


int main(int argc, char** argv) {
    MPI_Init(&argc, &argv);

    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    try {
        string seq1, seq2;
        if (rank == 0) {
            processFasta("seq1.fasta", seq1);
            processFasta("seq2.fasta", seq2);
            cout << "Sequence 1 Length: " << seq1.size() << " bases\n";
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
                    simdGlobalAlign(seq1, seq2); // Only rank 0 runs global
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

