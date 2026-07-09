//
// parse a bam file to find scaffolded contigs separated across a different scaffolded assembly
// originating from the same underlying data
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

bool
is_unmapped(vector<string> &parts){
    int flag = std::stoi(parts[1]);

    return flag & 0x4;
}


int
find_separated_contigs(const string &bamfile, unordered_map<string, vector<string>> &mappings){
    
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
            else if (parts.front()[0] == '@' || is_unmapped(parts)){
                line.clear();
                parts.clear();
                continue;
            }
            else {
                if (is_not_primary(parts)){
                    uint d = at_edge(parts[3], parts[5], chrom_lengths[parts[2]]);
                    if (d != 0){
                        string read = parts[0];
                        bool   add  = true;
                        for (int i = 0; i < secondaries[read].size(); i++){
                            if (secondaries[read][i].contig == parts[2]){
                                // these two have already been matched
                                add = false;
                                break;
                            }
                        }
                        if (add){
                            char dir = d == 1 ? '5' : '3';
                            edge_alignment p = {parts[2], dir};
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