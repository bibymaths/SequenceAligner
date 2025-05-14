/**
 * @file main.cpp
 * @brief Implements LCS, Global, and Local alignment for DNA sequences with colored output.
 * @author Abhinav Mishra
 * @date 2025
 */

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <algorithm>
#include <stdexcept>
#include <cstdlib>
using namespace std;

/// Match score for alignment
const int MATCH = 1;
/// Mismatch penalty for alignment
const int MISMATCH = -1;
/// Gap penalty for alignment
const int GAP = -1;
/// Width of output lines for formatted alignment printing
const int LINE_WIDTH = 150;

/// ANSI escape codes for colored output
#define RESET "\033[0m"
#define GREEN "\033[32m"
#define RED "\033[31m"
#define CYAN "\033[36m"


/**
 * @brief Print aligned sequences with color-coded differences.
 *
 * Matching characters are printed in green, mismatches in cyan,
 * and gaps in red. Outputs the sequences in aligned block format.
 *
 * @param seq1 First aligned sequence (may include gaps).
 * @param seq2 Second aligned sequence (may include gaps).
 */
void printColoredAlignment(const string &seq1, const string &seq2) {
    size_t length = seq1.size();
    for (size_t i = 0; i < length; i += LINE_WIDTH) {
        size_t end = min(i + LINE_WIDTH, length);

        cout << i + 1 << " - " << end << "\n";

        for (size_t j = i; j < end; j++) {
            if (seq1[j] == seq2[j]) {
                cout << GREEN << seq1[j] << RESET;
            } else if (seq1[j] == '-' || seq2[j] == '-') {
                cout << RED << seq1[j] << RESET;
            } else {
                cout << CYAN << seq1[j] << RESET;
            }
        }
        cout << "\n";

        for (size_t j = i; j < end; j++) {
            if (seq1[j] == seq2[j]) {
                cout << GREEN << seq2[j] << RESET;
            } else if (seq1[j] == '-' || seq2[j] == '-') {
                cout << RED << seq2[j] << RESET;
            } else {
                cout << CYAN << seq2[j] << RESET;
            }
        }
        cout << "\n\n";
    }
}

/**
 * @brief Reads a sequence from a FASTA file, skipping headers.
 *
 * @param filename Path to the input FASTA file.
 * @return string Concatenated sequence (no newlines or headers).
 * @throws runtime_error if the file cannot be opened.
 */
string readSequence(const string &filename) {
    ifstream file(filename);
    if (!file) {
        throw runtime_error("Error: Unable to open " + filename);
    }

    string line, sequence;
    while (getline(file, line)) {
        if (line.empty() || line[0] == '>') continue; // Skip FASTA headers
        sequence += line;
    }

    // Remove carriage returns in case of Windows-formatted files
    sequence.erase(remove(sequence.begin(), sequence.end(), '\r'), sequence.end());

    return sequence;
}


/**
 * @brief Perform Longest Common Subsequence (LCS) alignment.
 *
 * Prints the LCS length and the LCS string.
 *
 * @param x First input sequence.
 * @param y Second input sequence.
 */
void lcs(const string &x, const string &y) {
    int m = x.size(), n = y.size();
    vector<vector<int>> c(m + 1, vector<int>(n + 1, 0));
    vector<vector<char>> b(m + 1, vector<char>(n + 1, ' '));

    for (int i = 1; i <= m; i++) {
        for (int j = 1; j <= n; j++) {
            if (x[i - 1] == y[j - 1]) {
                c[i][j] = c[i - 1][j - 1] + 1;
                b[i][j] = 'D';
            } else if (c[i - 1][j] >= c[i][j - 1]) {
                c[i][j] = c[i - 1][j];
                b[i][j] = 'U';
            } else {
                c[i][j] = c[i][j - 1];
                b[i][j] = 'L';
            }
        }
    }

    int i = m, j = n;
    string lcs_str;
    while (i > 0 && j > 0) {
        if (b[i][j] == 'D') {
            lcs_str += x[i - 1];
            i--, j--;
        } else if (b[i][j] == 'U') {
            i--;
        } else {
            j--;
        }
    }
    reverse(lcs_str.begin(), lcs_str.end());
    cout << "LCS length: " << c[m][n] << "\nLCS: " << lcs_str << endl;
}

/**
 * @brief Perform Local Alignment using Smith-Waterman algorithm.
 *
 * Finds the best local alignment and prints colored alignment and score.
 *
 * @param x First input sequence.
 * @param y Second input sequence.
 */
void localAlign(const string &x, const string &y) {
    int m = x.size(), n = y.size();
    vector<vector<int>> dp(m + 1, vector<int>(n + 1, 0));
    int maxScore = 0, maxI = 0, maxJ = 0;

    for (int i = 1; i <= m; i++) {
        for (int j = 1; j <= n; j++) {
            int matchScore = (x[i - 1] == y[j - 1]) ? MATCH : MISMATCH;
            dp[i][j] = max({0, dp[i - 1][j - 1] + matchScore, dp[i - 1][j] + GAP, dp[i][j - 1] + GAP});
            if (dp[i][j] > maxScore) {
                maxScore = dp[i][j];
                maxI = i;
                maxJ = j;
            }
        }
    }

    string alignX, alignY;
    int i = maxI, j = maxJ;
    while (i > 0 && j > 0 && dp[i][j] > 0) {
        if (x[i - 1] == y[j - 1]) {
            alignX += x[i - 1];
            alignY += y[j - 1];
            i--, j--;
        } else if (dp[i][j] == dp[i - 1][j] + GAP) {
            alignX += x[i - 1];
            alignY += '-';
            i--;
        } else {
            alignX += '-';
            alignY += y[j - 1];
            j--;
        }
    }
    reverse(alignX.begin(), alignX.end());
    reverse(alignY.begin(), alignY.end());
    cout << "Local Alignment Score: " << maxScore << endl;
    printColoredAlignment(alignX, alignY);
}


/**
 * @brief Perform Global Alignment using Needleman-Wunsch algorithm.
 *
 * Finds the best global alignment and prints colored alignment and score.
 *
 * @param x First input sequence.
 * @param y Second input sequence.
 */
void globalAlign(const string &x, const string &y) {
    int m = x.size(), n = y.size();
    vector<vector<int>> dp(m + 1, vector<int>(n + 1, 0));
    vector<vector<char>> traceback(m + 1, vector<char>(n + 1, ' '));

    for (int i = 0; i <= m; i++) dp[i][0] = i * GAP;
    for (int j = 0; j <= n; j++) dp[0][j] = j * GAP;

    for (int i = 1; i <= m; i++) {
        for (int j = 1; j <= n; j++) {
            int matchScore = (x[i - 1] == y[j - 1]) ? MATCH : MISMATCH;
            int scores[] = {dp[i - 1][j - 1] + matchScore, dp[i - 1][j] + GAP, dp[i][j - 1] + GAP};
            int maxIndex = max_element(scores, scores + 3) - scores;
            dp[i][j] = scores[maxIndex];
            traceback[i][j] = (maxIndex == 0) ? 'D' : (maxIndex == 1) ? 'U' : 'L';
        }
    }

    int i = m, j = n;
    string alignX, alignY;
    while (i > 0 || j > 0) {
        if (i > 0 && j > 0 && traceback[i][j] == 'D') {
            alignX += x[i - 1];
            alignY += y[j - 1];
            i--, j--;
        } else if (i > 0 && traceback[i][j] == 'U') {
            alignX += x[i - 1];
            alignY += '-';
            i--;
        } else {
            alignX += '-';
            alignY += y[j - 1];
            j--;
        }
    }
    reverse(alignX.begin(), alignX.end());
    reverse(alignY.begin(), alignY.end());
    cout << "Global Alignment Score: " << dp[m][n] << endl;
    printColoredAlignment(alignX, alignY);
}

/**
 * @brief Main entry point. Parses arguments and selects alignment strategy.
 *
 * Usage: `./aligner <seq1_file> <seq2_file> <choice>`
 * Choice values:
 * - `1` = LCS
 * - `2` = Global alignment
 * - `3` = Local alignment
 *
 * @param argc Argument count.
 * @param argv Argument values.
 * @return int Exit status.
 */
int main(int argc, char* argv[]) {
    try {
        if (argc != 4) {
            cerr << "Usage: " << argv[0] << " <seq1_file> <seq2_file> <choice (1=LCS, 2=Global, 3=Local)>" << endl;
            return 1;
        }

        string seq1 = readSequence(argv[1]);
        string seq2 = readSequence(argv[2]);
        int choice = stoi(argv[3]);

        switch (choice) {
            case 1: lcs(seq1, seq2); break;
            case 2: globalAlign(seq1, seq2); break;
            case 3: localAlign(seq1, seq2); break;
            default: cout << "Invalid choice! Use 1 (LCS), 2 (Global), or 3 (Local)." << endl;
        }
    } catch (const exception &e) {
        cerr << e.what() << endl;
        return 1;
    }
    return 0;
}
