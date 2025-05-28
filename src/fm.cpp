#include <iostream>
#include <vector>
#include <string>
#include <algorithm> // For std::replace_if
#include <map>       // For std::map in print_c_table_like_info
#include <fstream>   // For file operations
#include <stdexcept> // For std::runtime_error
#include <cstdio>    // For std::remove (temporary file for sdsl construction)
#include <cctype>    // For std::isalnum, isprint
#include <iomanip>   // For std::hex, std::dec in printing non-printable chars
#include <ctime>     // For std::time in temporary file naming

// SDSL-Lite headers
#include <sdsl/csa_wt.hpp>
#include <sdsl/suffix_arrays.hpp>
#include <sdsl/util.hpp>

// -------- FMIndex Class using SDSL-Lite --------
class FMIndexSdsl {
public:
    using csa_sdsl_type = sdsl::csa_wt<sdsl::wt_huff<sdsl::rrr_vector<63>>, 32, 32,
                                     sdsl::sa_order_sa_sampling<>, sdsl::isa_sampling<>,
                                     sdsl::byte_alphabet>; // Using sdsl::byte_alphabet

    csa_sdsl_type fm_index;
    char sentinel_char_used;
    bool initialized;

    FMIndexSdsl() : sentinel_char_used('$'), initialized(false) {}

    FMIndexSdsl(const std::string& text, char sentinel = '$')
        : sentinel_char_used(sentinel), initialized(false) {
        std::string text_with_sentinel = text + sentinel_char_used;

        std::string temp_file_name = "temp_sdsl_input_text_" + std::to_string(std::time(nullptr)) + ".dat";
        {
            std::ofstream temp_out_file(temp_file_name, std::ios::binary);
            if (!temp_out_file) {
                throw std::runtime_error("Failed to create temporary file for SDSL construction: " + temp_file_name);
            }
            temp_out_file.write(text_with_sentinel.data(), text_with_sentinel.length());
        }

        sdsl::cache_config config(true, ".", "temp_sdsl_constructor_cache_");
        try {
            sdsl::construct(fm_index, temp_file_name, config, 1);
        } catch (const std::exception& e) {
            std::remove(temp_file_name.c_str());
            throw std::runtime_error(std::string("SDSL FM-Index construction failed: ") + e.what());
        }
        std::remove(temp_file_name.c_str());
        initialized = true;
    }

    // Implemented backward search manually using C-table and BWT.rank (Occ)
    std::pair<size_t, size_t> backward_search(const std::string& pattern) const {
        if (!initialized) throw std::runtime_error("FMIndexSdsl not initialized for search.");

        size_t l = 0;                     // SA range start (exclusive for rank, inclusive for C)
        size_t r = fm_index.size();       // SA range end (exclusive)

        if (pattern.empty()) {
            return {l, r};
        }

        for (auto it = pattern.rbegin(); it != pattern.rend(); ++it) {
            char original_char = *it;
            // value_type is the type of the compressed character in BWT
            typename csa_sdsl_type::value_type compressed_char =
                fm_index.char2comp[static_cast<unsigned char>(original_char)];

            // Check if character is part of the effective alphabet (sigma)
            // fm_index.C has size sigma. compressed_char should be < sigma.
            // If original_char was not in text, char2comp might map it to a value
            // outside the effective alphabet range [0, sigma-1] or to a default
            // (often 0, which might be the sentinel or another valid char).
            // A robust check would be to see if original_char maps to a comp_char < fm_index.sigma
            // and if that comp_char was actually present.
            // However, if C[comp_char] + rank(r, comp_char) <= C[comp_char] + rank(l, comp_char), range is empty.

            size_t new_l = fm_index.C[compressed_char] + fm_index.bwt.rank(l, compressed_char);
            size_t new_r = fm_index.C[compressed_char] + fm_index.bwt.rank(r, compressed_char);

            l = new_l;
            r = new_r;

            if (l >= r) { // Range is empty
                return {0, 0};
            }
        }
        return {l, r};
    }

    size_t count_occurrences(const std::string& pattern) const {
        if (!initialized) throw std::runtime_error("FMIndexSdsl not initialized for count.");
        if (pattern.empty()) return fm_index.size();
        // Alternatively, use the result of backward_search:
        // auto range = backward_search(pattern);
        // return range.second - range.first;
        return sdsl::count(fm_index, pattern.begin(), pattern.end());
    }

    // Corrected return type to sdsl::int_vector<64>
    sdsl::int_vector<64> locate_occurrences(const std::string& pattern) const {
        if (!initialized) throw std::runtime_error("FMIndexSdsl not initialized for locate.");
        return sdsl::locate(fm_index, pattern.begin(), pattern.end());
    }

    bool save(const std::string& filepath) const {
        if (!initialized) {
            std::cerr << "Warning: Attempting to save an uninitialized FMIndexSdsl." << std::endl;
            return false;
        }
        try {
            return sdsl::store_to_file(fm_index, filepath);
        } catch (const std::exception& e) {
            std::cerr << "Error saving FMIndex to " << filepath << ": " << e.what() << std::endl;
            return false;
        }
    }

    bool load(const std::string& filepath) {
        try {
            sdsl::load_from_file(fm_index, filepath);
            initialized = true;
            return true;
        } catch (const std::exception& e) {
            std::cerr << "Error loading FMIndex from " << filepath << ": " << e.what() << std::endl;
            initialized = false;
            return false;
        }
    }

    size_t get_text_length() const {
        if (!initialized) return 0;
        return fm_index.size();
    }

    // Corrected to use fm_index.sigma
    size_t get_alphabet_size() const {
        if (!initialized) return 0;
        return fm_index.sigma; // sigma is the size of the effective alphabet [cite: 33]
    }

    // Corrected to use fm_index.sigma and fm_index.comp2char
    void print_c_table_like_info(std::ostream& os) const {
        if (!initialized) {
            os << "{FMIndexSdsl not initialized}";
            return;
        }
        os << "C Table (from sdsl::csa_wt): {";
        bool first = true;

        std::map<char, size_t> display_map;
        // fm_index.sigma is the size of the effective alphabet [cite: 33]
        // fm_index.comp2char maps compressed codes to original characters [cite: 33]
        // fm_index.C is indexed by compressed character codes [cite: 33]
        for (size_t i = 0; i < fm_index.sigma; ++i) { // i is the compressed char code
            char original_char = fm_index.comp2char[i];
            display_map[original_char] = fm_index.C[i];
        }

        for (const auto& pair : display_map) {
             if (!first) os << ", ";
             if (isprint(static_cast<unsigned char>(pair.first))) {
                 os << "'" << pair.first << "': " << pair.second;
             } else {
                 os << "'\\x" << std::hex << std::setw(2) << std::setfill('0')
                    << static_cast<int>(static_cast<unsigned char>(pair.first))
                    << std::dec << "': " << pair.second;
             }
             first = false;
        }
        os << "}";
    }
};


// -------- Main Function --------
int main(int argc, char* argv[]) {
    if (argc < 2) {
        std::cerr << "Usage: " << argv[0] << " <fasta_file_or_-> [-s SENTINEL_CHAR]" << std::endl;
        return 1;
    }

    std::string fasta_filepath = argv[1];
    char sentinel_char_arg = '$';

    for (int i = 2; i < argc; ++i) {
        std::string arg = argv[i];
        if (arg == "-s" && i + 1 < argc) {
            if (std::string(argv[i+1]).length() != 1) {
                std::cerr << "Error: Sentinel must be a single character." << std::endl;
                return 1;
            }
            sentinel_char_arg = argv[i+1][0];
            i++;
        } else {
            std::cerr << "Warning: Unknown argument or missing value for -s: '" << arg << "'" << std::endl;
        }
    }

    std::istream* input_stream_ptr;
    std::ifstream file_stream_obj;

    if (fasta_filepath == "-") {
        input_stream_ptr = &std::cin;
        std::cerr << "Reading FASTA from stdin..." << std::endl;
    } else {
        file_stream_obj.open(fasta_filepath);
        if (!file_stream_obj.is_open()) {
            std::cerr << "Error: Cannot open FASTA file: " << fasta_filepath << std::endl;
            return 1;
        }
        input_stream_ptr = &file_stream_obj;
        std::cerr << "Reading FASTA from " << fasta_filepath << "..." << std::endl;
    }
    std::istream& input_stream = *input_stream_ptr;

    std::string current_header;
    std::string current_sequence_lines;
    std::string line;
    bool in_sequence = false;
    const size_t EXPECTED_SEQUENCE_RESERVE = 50000000;

    auto process_sequence = [&](const std::string& header, const std::string& sequence_data) {
        if (sequence_data.empty()) {
            std::cerr << "[" << header << "] Skipping empty sequence." << std::endl;
            return;
        }
        std::cerr << "Processing sequence: " << header
                  << " (length " << sequence_data.length() << ")" << std::endl;

        std::string sanitized_header = header;
        std::replace_if(sanitized_header.begin(), sanitized_header.end(),
                        [](char c){ return !std::isalnum(c) && c != '_' && c != '-'; },
                        '_');
        if (sanitized_header.empty()) sanitized_header = "sequence";

        try {
            FMIndexSdsl idx(sequence_data, sentinel_char_arg);
            std::string fname = sanitized_header + ".sdsl.fmidx";

            if (!idx.save(fname)) {
                 std::cerr << "Error: Could not save FM-Index to " << fname << std::endl;
            } else {
                std::cerr << "[" << header << "] Text length (incl. sentinel)=" << idx.get_text_length();
                std::cerr << ", Alphabet size (sigma)=" << idx.get_alphabet_size() << " ";
                idx.print_c_table_like_info(std::cerr);
                std::cerr << " (saved -> " << fname << ")" << std::endl;
                 // --- Example Usage (Optional: Uncomment to test search functionality) ---
                if (idx.get_text_length() > 10) {
                    std::string test_pattern;
                    if (idx.get_text_length() > 50) {
                        test_pattern = sequence_data.substr(std::min((size_t)20, sequence_data.length()-1),
                                                           std::min((size_t)3, sequence_data.length() - std::min((size_t)20, sequence_data.length()-1)));
                    } else {
                         test_pattern = sequence_data.substr(0, std::min((size_t)3, sequence_data.length()));
                    }

                    if (!test_pattern.empty()) {
                        std::cerr << "  Testing with pattern: \"" << test_pattern << "\"" << std::endl;
                        auto sa_range = idx.backward_search(test_pattern);
                        size_t count_val = idx.count_occurrences(test_pattern); // Renamed from count to avoid conflict
                        std::cerr << "    Pattern '" << test_pattern << "': SA range [" << sa_range.first
                                  << ", " << sa_range.second << "), Count (sdsl::count): " << count_val << std::endl;

                        if (count_val > 0 && count_val < 15) {
                            sdsl::int_vector<64> locs = idx.locate_occurrences(test_pattern);
                            std::cerr << "    Locations (0-indexed): ";
                            for(size_t j=0; j<locs.size(); ++j) {
                                std::cerr << locs[j] << (j == locs.size()-1 ? "" : ", ");
                            }
                            std::cerr << std::endl;
                        } else if (count_val >= 15) {
                            std::cerr << "    Pattern occurs " << count_val << " times (locations not printed)." << std::endl;
                        }
                    }
                }
            }
        } catch (const std::exception& e) {
            std::cerr << "ERROR processing sequence " << header << ": " << e.what() << std::endl;
        }
    };

    while(std::getline(input_stream, line)) {
        if (line.empty()) continue;

        if (!line.empty() && line.back() == '\r') {
            line.pop_back();
        }

        if (line[0] == '>') {
            if (in_sequence) {
                process_sequence(current_header, current_sequence_lines);
            }
            current_header = line.substr(1);
            current_sequence_lines.clear();
            current_sequence_lines.reserve(EXPECTED_SEQUENCE_RESERVE);
            in_sequence = true;
        } else if (in_sequence) {
            current_sequence_lines += line;
        }
    }

    if (in_sequence) {
        process_sequence(current_header, current_sequence_lines);
    }

    if (fasta_filepath != "-") {
        file_stream_obj.close();
    }

    return 0;
}