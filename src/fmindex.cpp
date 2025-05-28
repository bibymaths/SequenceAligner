#include <iostream>
#include <vector>
#include <string>
#include <algorithm>
#include <map>
#include <numeric> // For std::iota
#include <fstream> // For file operations
#include <stdexcept> // For std::runtime_error

// -------- Suffix Array Construction --------
std::vector<int> suffix_array(const std::string& s) {
    int n = s.length();
    std::vector<int> sa(n);
    std::iota(sa.begin(), sa.end(), 0); // sa = [0, 1, ..., n-1]

    std::vector<int> rank(n);
    for (int i = 0; i < n; ++i) {
        rank[i] = static_cast<unsigned char>(s[i]); // Initial ranks from char values
    }

    std::vector<int> tmp_rank(n);
    int k = 1;

    while (true) {
        // Sort suffix array based on (rank[i], rank[i+k])
        std::sort(sa.begin(), sa.end(), [&](int a, int b) {
            if (rank[a] != rank[b]) {
                return rank[a] < rank[b];
            }
            int rank_a_k = (a + k < n) ? rank[a + k] : -1;
            int rank_b_k = (b + k < n) ? rank[b + k] : -1;
            return rank_a_k < rank_b_k;
        });

        // Compute new ranks based on sorted suffix array
        tmp_rank[sa[0]] = 0;
        for (int i = 1; i < n; ++i) {
            int prev = sa[i - 1];
            int curr = sa[i];
            bool new_rank = (rank[curr] != rank[prev]);
            if (!new_rank) {
                 int rank_prev_k = (prev + k < n) ? rank[prev + k] : -1;
                 int rank_curr_k = (curr + k < n) ? rank[curr + k] : -1;
                 if (rank_curr_k != rank_prev_k) {
                     new_rank = true;
                 }
            }
            tmp_rank[curr] = tmp_rank[prev] + (new_rank ? 1 : 0);
        }
        rank.swap(tmp_rank);

        if (rank[sa[n - 1]] == n - 1) { // All ranks are unique
            break;
        }
        k <<= 1; // Double k
    }
    return sa;
}

// -------- FMIndex Class --------
class FMIndex {
public:
    std::string text_with_sentinel;
    std::vector<int> sa;
    std::string bwt;
    char sentinel_char;
    std::map<char, int> C;
    std::map<char, std::vector<int>> Occ;

    FMIndex(const std::string& text, char sentinel = '$')
        : sentinel_char(sentinel) {
        text_with_sentinel = text + sentinel_char;

        // Build and save the suffix array on the full text
        this->sa = suffix_array(text_with_sentinel);

        // Build the BWT directly from that SA
        this->bwt.reserve(text_with_sentinel.length());
        for (int i = 0; i < this->sa.size(); ++i) {
            if (this->sa[i] == 0) {
                this->bwt += text_with_sentinel.back();
            } else {
                this->bwt += text_with_sentinel[this->sa[i] - 1];
            }
        }

        build_c_table();
        build_occ_table();
    }

    // Private methods made public for easier access or could be helper functions
    void build_c_table() {
        // C[c] = total number of chars in text < c
        std::map<char, int> counts;
        for (char ch : this->bwt) {
            counts[ch]++;
        }
        this->C.clear();
        int total = 0;
        // std::map iterates keys in sorted order
        for (const auto& pair : counts) {
            this->C[pair.first] = total;
            total += pair.second;
        }
    }

    void build_occ_table() {
        // Occ[ch][i] = number of occurrences of ch in BWT[0:i-1] (prefix of length i)
        this->Occ.clear();
        // Initialize zero-row for all characters present in C table (and thus in BWT)
        for (const auto& pair : this->C) {
            this->Occ[pair.first].push_back(0);
        }
        // Also initialize for any character not in C yet (though unlikely if C is from BWT)
        for (char ch_in_bwt : this->bwt) {
            if (this->Occ.find(ch_in_bwt) == this->Occ.end()) {
                 this->Occ[ch_in_bwt].push_back(0);
            }
        }


        for (int i = 0; i < this->bwt.length(); ++i) {
            char current_char_in_bwt = this->bwt[i];
            for (auto& pair : this->Occ) {
                char c = pair.first;
                int prev_count = pair.second.back();
                pair.second.push_back(prev_count + (c == current_char_in_bwt ? 1 : 0));
            }
        }
    }

    std::pair<int, int> backward_search(const std::string& pattern) const {
        int l = 0;
        int r = bwt.length();
        for (int i = pattern.length() - 1; i >= 0; --i) {
            char ch = pattern[i];
            auto c_it = C.find(ch);
            auto occ_it = Occ.find(ch);

            if (c_it == C.end() || occ_it == Occ.end()) {
                return {0, 0}; // Character not in alphabet
            }

            // Ensure l and r are valid indices for Occ table
            if (l >= occ_it->second.size() || r >= occ_it->second.size()) {
                 // This case should ideally not happen if Occ is built correctly up to bwt.length()
                 return {0,0}; // Index out of bounds
            }

            l = c_it->second + occ_it->second[l];
            r = c_it->second + occ_it->second[r];

            if (l >= r) {
                return {0, 0}; // Empty interval
            }
        }
        return {l, r};
    }

    // Serialization method (simple binary format)
    void save(std::ostream& os) const {
        // Write text_with_sentinel
        size_t len = text_with_sentinel.length();
        os.write(reinterpret_cast<const char*>(&len), sizeof(len));
        os.write(text_with_sentinel.data(), len);

        // Write sentinel_char
        os.write(reinterpret_cast<const char*>(&sentinel_char), sizeof(sentinel_char));

        // Write sa
        len = sa.size();
        os.write(reinterpret_cast<const char*>(&len), sizeof(len));
        os.write(reinterpret_cast<const char*>(sa.data()), len * sizeof(int));

        // Write bwt
        len = bwt.length();
        os.write(reinterpret_cast<const char*>(&len), sizeof(len));
        os.write(bwt.data(), len);

        // Write C table
        len = C.size();
        os.write(reinterpret_cast<const char*>(&len), sizeof(len));
        for (const auto& pair : C) {
            os.write(reinterpret_cast<const char*>(&pair.first), sizeof(char));
            os.write(reinterpret_cast<const char*>(&pair.second), sizeof(int));
        }

        // Write Occ table
        len = Occ.size();
        os.write(reinterpret_cast<const char*>(&len), sizeof(len));
        for (const auto& pair : Occ) {
            os.write(reinterpret_cast<const char*>(&pair.first), sizeof(char));
            size_t vec_len = pair.second.size();
            os.write(reinterpret_cast<const char*>(&vec_len), sizeof(vec_len));
            os.write(reinterpret_cast<const char*>(pair.second.data()), vec_len * sizeof(int));
        }
    }
};


// -------- FASTA Parsing --------
// Reads the next FASTA sequence from the input stream.
// Returns true if a sequence was read, false otherwise (e.g., EOF).
bool read_fasta_sequence(std::istream& in, std::string& header, std::string& sequence) {
    std::string line;
    header.clear();
    sequence.clear();

    // Find the next header line
    while (std::getline(in, line)) {
        if (!line.empty() && line[0] == '>') {
            header = line.substr(1);
            // Remove potential carriage returns if files are from Windows/Mac
            if (!header.empty() && header.back() == '\r') {
                header.pop_back();
            }
            break;
        }
    }

    if (header.empty() && in.eof()) { // No more sequences
        return false;
    }
    if (header.empty() && !in.eof()){ // Found content before a header or malformed
        // This could be an error, or skip to next actual header.
        // For this port, we'll assume we found a header if loop exited.
        // If header is still empty here and not EOF, it's likely an issue
        // with the FASTA or we hit EOF immediately after a header.
         if(in.eof() && !header.empty()) { /* Fine, last sequence had no seq lines */ }
         else if (in.eof() && header.empty()) return false; // Truly end of file.
         // else: header is empty, but not eof. This is odd.
    }


    // Read sequence lines
    while (std::getline(in, line)) {
        if (line.empty()) continue;
        if (line[0] == '>') {
            // Found the next header, put it back for the next call
            // This requires the stream to support putback.
            // A simple way is to rewind by the length of the line + newline.
            for (size_t i = 0; i <= line.length(); ++i) { // line.length() + 1 for newline
                in.unget(); // or in.seekg(- (line.length() + 1), std::ios_base::cur);
                            // unget is safer for basic_istream
            }
            // Correcting unget loop:
            // Actually, it's simpler if the main loop checks if stream is good and
            // the first line read is already processed or it is the first line of the new sequence
            // For now, we'll assume `std::getline` followed by `in.peek()` or careful `unget`
            // A simpler strategy for this parser:
            // After reading a header, read lines until another '>' or EOF.
            // The current getline for sequence will consume the line.
            // If the next line IS a header, we need to put it back.
            // The python version reads until next '>' and then `yield`.
            // The C++ version is called iteratively.

            // Simpler way: read lines and if it's not a header, append.
            // The next call to read_fasta_sequence will handle the line starting with '>'.
            break; // End of current sequence
        }
        // Remove potential carriage returns
        if (!line.empty() && line.back() == '\r') {
            line.pop_back();
        }
        sequence += line;
    }
    return !header.empty(); // Successfully read if header is not empty
}


// -------- Main Function --------
int main(int argc, char* argv[]) {
    if (argc < 2) {
        std::cerr << "Usage: " << argv[0] << " <fasta_file_or_-> [-s SENTINEL_CHAR]" << std::endl;
        return 1;
    }

    std::string fasta_filepath = argv[1];
    char sentinel_char = '$';

    for (int i = 2; i < argc; ++i) {
        std::string arg = argv[i];
        if (arg == "-s" && i + 1 < argc) {
            if (std::string(argv[i+1]).length() != 1) {
                std::cerr << "Error: Sentinel must be a single character." << std::endl;
                return 1;
            }
            sentinel_char = argv[i+1][0];
            i++; // Consume sentinel value
        } else {
            std::cerr << "Warning: Unknown argument '" << arg << "'" << std::endl;
        }
    }

    std::istream* input_stream;
    std::ifstream file_stream;

    if (fasta_filepath == "-") {
        input_stream = &std::cin;
        std::cerr << "Reading FASTA from stdin..." << std::endl;
    } else {
        file_stream.open(fasta_filepath);
        if (!file_stream.is_open()) {
            std::cerr << "Error: Cannot open FASTA file: " << fasta_filepath << std::endl;
            return 1;
        }
        input_stream = &file_stream;
        std::cerr << "Reading FASTA from " << fasta_filepath << "..." << std::endl;
    }

    std::string current_header;
    std::string current_sequence_lines;
    std::string line;

    bool in_sequence = false;

    while(std::getline(*input_stream, line)) {
        if (line.empty()) continue;

        if (line[0] == '>') {
            if (in_sequence) { // Process previous sequence
                // Sanitize header for filename
                std::string sanitized_header = current_header;
                std::replace_if(sanitized_header.begin(), sanitized_header.end(),
                                [](char c){ return !std::isalnum(c) && c != '_' && c != '-'; },
                                '_');
                if (sanitized_header.empty()) sanitized_header = "sequence";


                FMIndex idx(current_sequence_lines, sentinel_char);
                std::string fname = sanitized_header + ".fmidx";
                std::ofstream outfile(fname, std::ios::binary);
                if (!outfile) {
                    std::cerr << "Error: Could not open " << fname << " for writing." << std::endl;
                } else {
                    idx.save(outfile);
                    outfile.close();
                }

                std::cerr << "[" << current_header << "] BWT length=" << idx.bwt.length() << "  C={";
                bool first_c = true;
                for (const auto& pair : idx.C) {
                    if (!first_c) std::cerr << ", ";
                    std::cerr << "'" << pair.first << "': " << pair.second;
                    first_c = false;
                }
                std::cerr << "}  (saved -> " << fname << ")" << std::endl;
            }
            // Start new sequence
            current_header = line.substr(1);
            // Remove potential carriage returns
            if (!current_header.empty() && current_header.back() == '\r') {
                current_header.pop_back();
            }
            current_sequence_lines.clear();
            in_sequence = true;
        } else if (in_sequence) {
            // Remove potential carriage returns
            if (!line.empty() && line.back() == '\r') {
                line.pop_back();
            }
            current_sequence_lines += line;
        }
    }

    // Process the last sequence in the file (if any)
    if (in_sequence && !current_header.empty()) {
         std::string sanitized_header = current_header;
         std::replace_if(sanitized_header.begin(), sanitized_header.end(),
                        [](char c){ return !std::isalnum(c) && c != '_' && c != '-'; },
                        '_');
        if (sanitized_header.empty()) sanitized_header = "sequence";


        FMIndex idx(current_sequence_lines, sentinel_char);
        std::string fname = sanitized_header + ".fmidx";
        std::ofstream outfile(fname, std::ios::binary);
         if (!outfile) {
            std::cerr << "Error: Could not open " << fname << " for writing." << std::endl;
        } else {
            idx.save(outfile);
            outfile.close();
        }

        std::cerr << "[" << current_header << "] BWT length=" << idx.bwt.length() << "  C={";
        bool first_c = true;
        for (const auto& pair : idx.C) {
            if (!first_c) std::cerr << ", ";
            std::cerr << "'" << pair.first << "': " << pair.second;
            first_c = false;
        }
        std::cerr << "}  (saved -> " << fname << ")" << std::endl;
    }


    if (fasta_filepath != "-") {
        file_stream.close();
    }

    return 0;
}