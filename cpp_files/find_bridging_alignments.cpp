//
// parse a bam file to find reads bridging 2 seperate contigs
//

#include <iostream>
#include <fstream>
#include <algorithm>
#include <sys/stat.h>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>
#include <cstdio>
#include <cctype>

using std::string;
using std::unordered_map;
using std::unordered_set;
using std::vector;
using std::fstream;
using std::iostream;
using std::cout;
using std::cerr;

typedef unsigned int uint;

struct edge_alignment {
    string contig;
    char   direction;
    uint   length;

    // define a constructor
    edge_alignment(string cntg, char d, uint len) :
        contig(cntg), direction(d), length(len) {}
};

int
parse_tabular(string &line, vector<string> &parts){

    int start = 0, end = 0;

    while (end < line.size()){
        if (line[end] == '\t'){
            parts.emplace_back(line.substr(start, end - start));
            start = end + 1;
        }
        end++;
    }

    if (start < line.size()){
        parts.emplace_back(line.substr(start));
    }

    return 0;
}

bool
file_exists(const string &path){
    struct stat buffer;
    return stat(path.c_str(), &buffer) == 0;
}

bool
check_for_index(const string& bamfile){
    //
    // check if the bam file has an index
    // file in the same path
    //
    return file_exists(bamfile + ".bai") || file_exists(bamfile + ".csi");
}

int
at_edge(string &pos, string &cigar, const uint &chrom_length){

    char first = '\0', last = '\0';
    int length = 0, from_start = 0, from_end = 0;
    int align_pos = std::stoi(pos);
    int cumulative = align_pos;

    int p  = 0;
    char c;
    string lstr; // length string

    while (p < cigar.size()){
        if (std::isdigit(cigar[p])){
            lstr += cigar[p];
        }
        else {
            if (lstr.size() > 0){
                length = std::stoi(lstr);
                c      = cigar[p];
                lstr.clear();
                if (first == '\0'){
                    first      = c;
                    last       = c;
                    from_start = length;
                    from_end   = length;
                }
                else {
                    last     = c;
                    from_end = length;
                }
                if (c != 'I' && c != 'S' && c != 'H'){
                    cumulative += length;
                }
            }
        }
        p++;
    }

    //
    // now check if the first or last characters are clippings
    //

    if (first == 'S' || first == 'H'){
        if ((align_pos - from_start) < 0){
            return 1; // left
        }
    }
    if (last == 'S' || last == 'H'){
        if (cumulative + from_end > chrom_length){
            return 2; // right
        }
    }

    return 0; // neither
}

uint
get_query_length(const string &cigar, bool const primary){
    //
    // helper function to get the length of the alignment
    //

    uint p = 0, n = cigar.size(), last_p = 0, size = 0;
    string lstr; // string representation of the length
    char c;
    const char o = primary ? 'S' : 'D';

    uint aligned_len = 0;

    while (p != n){
        c = cigar[p];
        if (isdigit(c)){
            p++;
        }
        else {
            if (c == 'M' || c == 'I' || c == o){
                lstr = cigar.substr(last_p, p - last_p);
                size = stoi(lstr);
                aligned_len += size;
            }
            last_p = p + 1;
            p++;
        }
    }

    return aligned_len;
}

bool
is_not_primary(const string &samflag){
    //
    // check if the alignment is not primary
    //

    int flag = std::stoi(samflag);

    if (flag & 0x100 || flag & 0x800){
        return true; }

    return false;
}

bool
is_unmapped(const string &samflag){
    int flag = std::stoi(samflag);

    return flag & 0x4;
}

int
get_sequence_length(vector<string> &parts, unordered_map<string, uint> &chrom_lengths){
    //
    // add the sequence length to the map
    //

    // start on the 4th index after the : character
    string seq_name         = parts[1].substr(3); 
    uint   seq_length       = std::stoul(parts[2].substr(3));
    chrom_lengths[seq_name] = seq_length;

    return 0;
}

int
find_secondary_alignments(const string &bamfile, unordered_map<string, uint> &chrom_lengths,
                          unordered_map<string, vector<edge_alignment>> &secondaries){
    
    //
    // collect the name of sequences with clips at the ends of the contigs
    //
    string line;
    vector<string> parts;

    //
    // create a file stream
    //
    string cmd = "samtools view -h " + bamfile;
    FILE *fh   = popen(cmd.c_str(), "r");

    int c;
    //
    // We will have to build the line
    //
    while ((c = fgetc(fh)) != EOF){
        if (c == '\n'){
            parse_tabular(line, parts);

            if (parts.empty()){
                line.clear();
                continue;
            }
            if (parts.front() == "@SQ"){
                get_sequence_length(parts, chrom_lengths);
            }
            else if (parts.front()[0] == '@' || is_unmapped(parts[1])){
                line.clear();
                parts.clear();
                continue;
            }
            else {
                string samflag = parts[1];
                string cigar   = parts[5];
                if (is_not_primary(samflag)){
                    uint d = at_edge(parts[3], cigar, chrom_lengths[parts[2]]);
                    if (d != 0){
                        string   read = parts[0];
                        bool      add = true;
                        char      dir = d == 1 ? '5' : '3';
                        uint algn_len = get_query_length(cigar, false);
                        for (int i = 0; i < secondaries[read].size(); i++){
                            if (secondaries[read][i].contig == parts[2]){
                                // these two have already been matched
                                if (secondaries[read][i].length < algn_len){
                                    secondaries[read][i].length = algn_len;
                                    secondaries[read][i].direction = dir;
                                }
                                add = false;
                                break;
                            }
                        }
                        if (add){
                            edge_alignment p = {parts[2], dir, algn_len};
                            secondaries[read].push_back(p);
                        }
                    }
                }
            }
            parts.clear();
            line.clear();
        }
        else {
            line += static_cast<char>(c);
        }
    }

    pclose(fh);

    return 0;
}

int
find_edge_primaries(const string &bamfile, unordered_map<string, uint> &chrom_lengths,
                    unordered_map<string, vector<edge_alignment>> &secondaries){
    //
    // loop back through the file and check if edge sequences for their primary alignments
    //
    string line;
    vector<string> parts;
    std::ofstream ofh("Edge_Sequences.tsv");

    ofh << "#Sequence\tPrimary.Contig\tSecondary.Contigs\n";

    //
    // create a file stream
    //
    string cmd = "samtools view -h " + bamfile;
    FILE *fh   = popen(cmd.c_str(), "r");

    int c;
    //
    // We will have to build the line
    //
    while ((c = fgetc(fh)) != EOF){
        if (c == '\n'){
            parse_tabular(line, parts);

            if (parts.empty()){
                line.clear();
                continue;
            }
            if (parts.front()[0] == '@' || secondaries.count(parts[0]) == 0 || is_not_primary(parts[1]) || is_unmapped(parts[1])){
                line.clear();
                parts.clear();
                continue;
            }
            //
            // primary sequence is at the edge and on a different sequence
            //
            string cigar = parts[5];
            uint d = at_edge(parts[3], cigar, chrom_lengths[parts[2]]);
            if (d != 0){
                bool same    = false;
                string read = parts[0];
                for (int i = 0; i < secondaries[read].size(); i++){
                    if (secondaries[read][i].contig == parts[2]){
                        same = true;
                        break;
                    }
                }
                if (same == false){
                    uint length  = get_query_length(cigar, true);
                    uint aln_len = get_query_length(cigar, false);
                    char dir = d == 1 ? '5' : '3';
                    ofh << read << '\t' << parts[2] << '_' << length << '_' << aln_len << '_' << dir << '\t';
                    int count = 0, total = secondaries[read].size();
                    for (int i = 0; i < total; i++){
                        count++;
                        ofh << secondaries[read][i].contig << '_' << secondaries[read][i].length << '_' << secondaries[read][i].direction;
                        if (count < total){
                            ofh << ',';
                        }
                    }
                    ofh << '\n';
                    secondaries[read].clear();
                    secondaries.erase(read);
                }
            }
            parts.clear();
            line.clear();
        }
        else {
            line += static_cast<char>(c);
        }
    }

    pclose(fh);
    ofh.close();

    return 0;

}

int
main(int argc, char *args[]){

    string bamfile;

    for (int i = 1; i < argc; i++){
        if (string(args[i]) == "-b" && i + 1 < argc){
            bamfile = string(args[i + 1]);
        }
    }

    if (bamfile.size() == 0){
        cerr << "No bam file was provided. Program requires a bam file (-b bam)\n";
        exit(EXIT_FAILURE);
    }
    else if (check_for_index(bamfile) == false){
        cerr << "Could not locate a .bai or .csi index file for " << bamfile << '\n';
        exit(EXIT_FAILURE);
    }

    //
    // step 1: find sequences that are secondary & at an
    // edge of a sequence
    //
    unordered_map<string, uint> chrom_lengths;
    unordered_map<string, vector<edge_alignment>> secondaries;
    find_secondary_alignments(bamfile, chrom_lengths, secondaries);

    //
    // step 2: find the primary sequences that are
    // also at edges, but for different sequences
    //
    find_edge_primaries(bamfile, chrom_lengths, secondaries);

    return 0;
}