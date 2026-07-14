//
// parse a bam file to find scaffolded contigs from ntLink separated across a different scaffolded assembly
// and see if we can reconstruct the contig based on the alignments
//

#include <iostream>
#include <fstream>
#include <algorithm>
#include <sys/stat.h>
#include <string>
#include <unordered_map>
#include <vector>
#include <cstdio>
#include <cctype>
#include <cstdint>

using std::string;
using std::unordered_map;
using std::vector;
using std::fstream;
using std::iostream;
using std::cout;
using std::cerr;
using std::isdigit;
using std::uint8_t;

typedef unsigned int uint;


struct Align {
    string ref_seq;
    string cigar;
    uint   pos;
    char   direction;
    uint8_t id;
};

struct Part {
    uint8_t id;
    uint  start; // relative to the read
    uint  length;
};


class Contig {

private:
    vector<Align>          _aligns; // hold onto each alignment
    vector<vector<Align>> _full_aligns; // contains only alignments that can represent this full contig
    uint8_t               _count = 0;

public:

    string name;
    uint   size;

    Contig(string &n){
        this->name = n;
    }

    ~Contig(void){
        this->_aligns.clear();
        this->_full_aligns.clear();
    }

    void add_align(Align &a){
        this->_aligns.emplace_back(a);
    }

    uint8_t get_count(void){
        this->_count++;
        return this->_count;
    }

    void reconstruct(void){
        //
        // step 1: we need to construct parts
        //
        vector<Part> _parts;
        const uint threshold = 500;

        for (int i = 0; i < this->_aligns.size(); i++){
            Align &aln          = this->_aligns[i];
            const string &cigar = aln.cigar;
            int  start          = -1;
            uint p = 0, n = 0, last_p = 0, size = 0;
            string lstr; // string representation of the length
            char c;

            while (p != n){
                c = cigar[p];
                if (isdigit(c)){
                    p++;
                }
                else {
                    if (start == -1){
                        // first encounter of a feature
                        if (c == 'H' || c == 'S'){
                            lstr = cigar.substr(p, last_p - p);
                            size = std::stoi(lstr);
                            if (aln.direction == 'R'){
                                start = this->size - size;
                            }
                            else {
                                start = size;
                            }
                        }
                        // if first feature is a match
                        else if (c == 'M'){
                            start = 1;
                        }
                    }
                    else if (c != 'H' || c != 'S' || c != 'I'){
                        lstr  = cigar.substr(p, last_p - p);
                        size += std::stoi(lstr);
                    }
                    last_p = p + 1; // skip this char
                }
            }
            
            // now create the part
            _parts.emplace_back(aln.id, start, size);
        } // end of i loop

        std::sort(_parts.begin(), _parts.end(), [](const Part &p1, const Part &p2){
            return p1.start < p2.start;});

        if (_parts.size() == 1){
            if (_parts[0].length >= this->size - threshold){
                vector<Part> v = {_parts[0]};
                this->_full_aligns.emplace_back(v);
            }
        }
        // now see if we can reconstruct a full alignment
        Part prev = _parts[0], next;
        const uint dist = 100;
        vector<vector<Part>> components;
        vector<Part> component = {prev};
        for (int i = 1; i < _parts.size(); i++){
            next = _parts[i];
            uint next_pos = component.back().start + component.back().length - 1;
            if (next.start >= next_pos - dist){
                component.push_back(next);
            }
            else {
                components.push_back(component);
                component = {next};
            }
        }
        

    }

    bool has_reconstruction(void){
        return !this->_full_aligns.empty();
    }


};

int
parse_tabular(string &line, vector<string> &parts){

    int start = 0, end = 0, n = line.size();

    while (end < n){
        if (line[end] == '\t'){
            parts.emplace_back(line.substr(start, end - start));
            start = end + 1;
        }
        end++;
    }

    if (start < n){
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
get_sequence_length(const string &cigar){

    //
    // parse the cigar to determine the sequence length
    //

    uint size = 0, p = 0, n = cigar.size(), last_p = 0;
    string cur_str;

    while (p < n){
        char c = cigar[p];
        if (isdigit(c)){
            p++;
        }
        else{
            if (c != 'I'){
                cur_str = cigar.substr(last_p, p - last_p);
                size += std::stoi(cur_str);
            }
            last_p = p + 1; // skip this char
        }
    }
}

bool
parse_cigar(string &cigar){

    // determine if this cigar indicates that this alignment could work

    char first = '\0', last = '\0';
    const int size = 1000;

    int p = 0, n = cigar.size(), length;
    char c;
    string lstr; // length string
    bool clipped = false;
    bool large_align = false;

    while (p < n){
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
    //
    // check if this is an unmapped record
    //

    return std::stoi(flag) & 0x4;
}

bool
is_primary(const string &flag){
    //
    // check if the alignment is primary
    //

    int enc = std::stoi(flag); // encoding

    if (enc & 0x100 || enc & 0x800){
        return false; }

    return true;
}


bool
reverse_strand(const string &flag){
    //
    // determine if this alignment is on the reverse strand
    //

    return (std::stoi(flag) & 0x10) != 0;
}

int
collect_alignments(const string &bamfile, vector<Contig> &mappings){
    
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
                Align aln;
                string name   = parts[0];
                aln.cigar     = parts[5];
                aln.ref_seq   = parts[2];
                aln.direction = reverse_strand(parts[1]) ? 'R' : 'F';
                aln.pos       = static_cast<uint>(std::stoi(parts[3]));
                for (int i = 0; i < mappings.size(); i++){
                    if (mappings[i].name == name){
                        aln.id = mappings[i].get_count();
                        mappings[i].add_align(aln);
                        name.clear();

                        // check if we need to get the sequence length
                        if (is_primary(parts[1])){
                            mappings[i].size = get_sequence_length(aln.cigar);
                        }
                        break;
                    }
                }
                // did we not add this?
                if (!name.empty()){
                    Contig c(name);
                    aln.id = c.get_count();
                    c.add_align(aln);
                    if (is_primary(parts[1])){
                        c.size = get_sequence_length(aln.cigar);
                    }
                    mappings.push_back(c);
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
    vector<Contig> mappings;
    collect_alignments(bamfile, mappings);

    //
    // step 2: write these out to a file
    //
    write_mappings(mappings);


    return 0;
}