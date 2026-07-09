//
// parse a bam file to find scaffolded contigs from ntLink separated across a different scaffolded assembly
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

struct ntLink_Align {
    string contig;
    uint   pos;
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

bool
parse_cigar(string &cigar){

    // determine if this cigar indicates that this alignment could work

    char first = '\0', last = '\0';
    const int size = 1000;

    int p = 0, length;
    char c;
    string lstr; // length string
    bool clipped = false;
    bool large_align = false;

    while (p < cigar.size()){
        if (std::isdigit(cigar[p])){
            lstr += cigar[p];
        }
        else {
            if (lstr.size() > 0){
                length = std::stoi(lstr);
                c      = cigar[p];
                lstr.clear();
                if (length >= size){
                    if (c == 'H' || c == 'S'){
                        clipped = true;
                    }
                    else if (c == 'M'){
                        large_align = true;
                    }
                }
            }
        }
        p++;
    }

    return large_align && clipped;

}

bool
is_unmapped(const string &flag){
    int flag_int = std::stoi(flag);

    return flag_int & 0x4;
}

bool
reverse_strand(const string &flag){
    int flag_int = std::stoi(flag);

    return (flag_int & 0x10) != 0;
}

int
find_separated_contigs(const string &bamfile, unordered_map<string, vector<ntLink_Align>> &mappings){
    
    //
    // find the ntLink scaffolds that map across different loci
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
            else if (parts.front()[0] != 'n' || is_unmapped(parts[1])){
                line.clear();
                parts.clear();
                continue;
            }
            else {
                string cigar     = parts[5];
                bool informative = parse_cigar(cigar);
                if (informative){
                    char dir = reverse_strand(parts[1]) ? 'r' : 'f';
                    uint pos  = static_cast<uint>(std::stoi(parts[3]));
                    ntLink_Align nt;
                    nt.contig    = parts[2];
                    nt.direction = dir;
                    nt.pos       = pos;
                    mappings[parts[0]].push_back(nt);
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
write_mappings(unordered_map<string, vector<ntLink_Align>> &mappings){

    //
    // create the output stream to construct the file
    //
    std::ofstream ofh("Separated.Contigs.txt");

    ofh << "#Contig.ID\tAlignments\n";

    for (auto itr = mappings.begin(); itr != mappings.end(); itr++){
        string scaffold           = itr->first;
        vector<ntLink_Align> &nts = itr->second;

        //
        // sort the vector of nt alignments
        //

        std::sort(nts.begin(), nts.end(), []
            (const ntLink_Align &a, const ntLink_Align &b){
                if (a.contig != b.contig){
                    return a.contig < b.contig; }
                if (a.pos != b.pos){
                    return a.pos < b.pos; }
                else {
                    return true;
                }
            });

        ofh << scaffold << '\t';

        int size = nts.size();
        for (int i = 0; i < size; i++){
            ntLink_Align &nt = nts[i];
            char c = i == size - 1 ? '\n' : ',';
            ofh << nt.contig << '_' << nt.pos << '_' << nt.direction << c;
        }

    } // end of itr

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
    // step 1: find contigs that were mapped across multiple
    // sequences
    //
    unordered_map<string, vector<ntLink_Align>> mappings;
    find_separated_contigs(bamfile, mappings);

    //
    // step 2: write these out to a file
    //
    write_mappings(mappings);


    return 0;
}