#include <iostream>
#include <fstream> 
#include <sstream>
#include <vector>
#include <algorithm>
#include <stdexcept>
using namespace std;

// Constants for scoring
const int MATCH = 1;
const int MISMATCH = -1;
const int GAP = -1;

// Function to read sequence from a file
string readSequence(const string &filename) {
    ifstream file(filename);
    if (!file) {
        throw runtime_error("Error: Unable to open " + filename);
    }
    stringstream buffer;
    buffer << file.rdbuf();
    string sequence = buffer.str();
    sequence.erase(remove(sequence.begin(), sequence.end(), '\n'), sequence.end());
    sequence.erase(remove(sequence.begin(), sequence.end(), '\r'), sequence.end());
    return sequence;
}
 
// Function to perform Longest Common Subsequence (LCS) alignment
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

// Function to perform Local Alignment
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
    cout << "Local Alignment Score: " << maxScore << "\n" << alignX << "\n" << alignY << endl;
} 


// Function to perform Global Alignment
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
    cout << "Global Alignment Score: " << dp[m][n] << "\n" << alignX << "\n" << alignY << endl;
}


int main() {
    try {
        string seq1 = readSequence("seq1.txt");
        string seq2 = readSequence("seq2.txt"); 
         
        cout << "Sequence 1: " << seq1 << endl;
        cout << "Sequence 2: " << seq2 << endl;

        int choice;
        cout << "Select Alignment Method:\n1. LCS\n2. Global\n3. Local\nChoice: ";
        cin >> choice;

        switch (choice) {
            case 1: lcs(seq1, seq2); break;
            case 2: globalAlign(seq1, seq2); break;
            case 3: localAlign(seq1, seq2); break;
            default: cout << "Invalid choice!" << endl;
        }
    } catch (const exception &e) {
        cerr << e.what() << endl;
        return 1;
    }
    return 0;
}
