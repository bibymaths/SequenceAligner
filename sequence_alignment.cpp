#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
using namespace std;

// Constants for scoring
const int MATCH = 3;
const int MISMATCH = -1;
const int GAP = -2;

// Function to read sequence from a file
string readSequence(const string &filename) {
    ifstream file(filename);
    if (!file) {
        cerr << "Error: Unable to open " << filename << endl;
        return "";
    }
    string sequence;
    getline(file, sequence);
    return sequence;
}

// Function to perform Longest Common Subsequence (LCS) alignment
void lcs(const string &x, const string &y) {
    int m = x.size(), n = y.size();
    vector<vector<int>> c(m + 1, vector<int>(n + 1, 0));

    for (int i = 1; i <= m; i++) {
        for (int j = 1; j <= n; j++) {
            if (x[i - 1] == y[j - 1]) {
                c[i][j] = c[i - 1][j - 1] + 1;
            } else {
                c[i][j] = max(c[i - 1][j], c[i][j - 1]);
            }
        }
    }
    cout << "LCS length: " << c[m][n] << endl;
}

// Function to perform Global Alignment
void globalAlign(const string &x, const string &y) {
    int m = x.size(), n = y.size();
    vector<vector<int>> dp(m + 1, vector<int>(n + 1, 0));

    for (int i = 0; i <= m; i++) dp[i][0] = i * GAP;
    for (int j = 0; j <= n; j++) dp[0][j] = j * GAP;

    for (int i = 1; i <= m; i++) {
        for (int j = 1; j <= n; j++) {
            int matchScore = (x[i - 1] == y[j - 1]) ? MATCH : MISMATCH;
            dp[i][j] = max({dp[i - 1][j - 1] + matchScore, dp[i - 1][j] + GAP, dp[i][j - 1] + GAP});
        }
    }
    cout << "Global Alignment Score: " << dp[m][n] << endl;
}

// Function to perform Local Alignment
void localAlign(const string &x, const string &y) {
    int m = x.size(), n = y.size();
    vector<vector<int>> dp(m + 1, vector<int>(n + 1, 0));
    int maxScore = 0;

    for (int i = 1; i <= m; i++) {
        for (int j = 1; j <= n; j++) {
            int matchScore = (x[i - 1] == y[j - 1]) ? MATCH : MISMATCH;
            dp[i][j] = max({0, dp[i - 1][j - 1] + matchScore, dp[i - 1][j] + GAP, dp[i][j - 1] + GAP});
            maxScore = max(maxScore, dp[i][j]);
        }
    }
    cout << "Local Alignment Score: " << maxScore << endl;
}

int main() {
    string seq1 = readSequence("seq1.txt");
    string seq2 = readSequence("seq2.txt");
    if (seq1.empty() || seq2.empty()) {
        cerr << "Error: Unable to read sequences from files." << endl;
        return 1;
    }

    int choice;
    cout << "Select Alignment Method:\n1. LCS\n2. Global\n3. Local\nChoice: ";
    cin >> choice;

    switch (choice) {
        case 1: lcs(seq1, seq2); break;
        case 2: globalAlign(seq1, seq2); break;
        case 3: localAlign(seq1, seq2); break;
        default: cout << "Invalid choice!" << endl;
    }

    return 0;
}
