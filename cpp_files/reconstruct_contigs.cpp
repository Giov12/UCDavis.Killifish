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
using std::stoi;
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


class Scaffold {

private:
    vector<Align>          _aligns; // hold onto each alignment
    vector<vector<Align>> _full_aligns; // contains only alignments that can represent this full contig
    uint8_t               _count = 0;

public:

    string name;
    uint   size;

    Scaffold(string &n){
        this->name = n;
    }

    ~Scaffold(void){
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
            uint p = 0, n = cigar.size(), last_p = 0, size = 0, length = 0;
            string lstr; // string representation of the length
            char c;

            while (p != n){
                c = cigar[p];
                if (isdigit(c)){
                    p++;
                }
                else {
                    lstr  = cigar.substr(last_p, p - last_p);
                    size  = stoi(lstr);
                    if (start == -1){
                        // first encounter of a feature
                        if (c == 'H' || c == 'S'){
                            if (aln.direction == 'R'){
                                start = this->size - size;
                            }
                            else {
                                start = size;
                            }
                            size = 0; // this does not contribute to the size of the part
                        }
                        // if first feature is a match
                        else if (c == 'M'){
                            start  = 1;
                            length = size;
                        }
                    }
                    else if (c == 'M' || c == 'I'){
                        length += size;
                    }
                    last_p = p + 1; // skip this char
                    p++;
                }
            }
            
            // now create the part
            _parts.emplace_back(aln.id, start, length);
        } // end of i loop

        //
        // step 2: try to reconstruct using adjacent parts to each other
        // emphasizing the largest parts over smaller parts (parts == aligned portions)
        //
        const uint nparts = _parts.size();
        
        if (nparts == 1){
            if (_parts[0].length >= this->size - threshold){
                vector<Part> v = {_parts[0]};
                this->_full_aligns.emplace_back(v);
            }
        }

        //
        // sort by starts so we can build the scaffold left to right
        //
        std::sort(_parts.begin(), _parts.end(), [](const Part &p1, const Part &p2){
            return p1.start < p2.start;});

        
        //
        // now see if we can reconstruct a full alignment
        //
        const uint   dist = 100;
        vector<uint> best_len(nparts);
        vector<int>  best_from(nparts, -1);

        //
        // this is a chaining problem, since we have a direction
        // lowest start to highest start, but we need to be mindful
        // of possible overlapping alignments
        //

        for (uint i = 0; i < nparts; i++){
            best_len[i] = _parts[i].length; // initialize with its own length
            for (uint j = 0; j < nparts; j++){
                uint j_end = _parts[j].start + _parts[j].length - 1;
                //
                // ensure that the ith part comes after j ends in a reasonable dist
                //
                if (_parts[i].start >= j_end && _parts[i].start <= j_end + dist){
                    uint candidate_len = best_len[j] + _parts[i].length;
                    if (candidate_len > best_len[i]){
                        best_len[i]   = candidate_len; // update the length
                        best_from[i] = j;
                    }
                }
            } // end of j
        } // end of i

        //
        // step 3: find where the best chain ends
        //
        uint best_end = 0;
        for (uint i = 1; i < nparts; i++){
            if (best_len[i] == best_len[best_end]){
                best_end = i;
            }
        }

        // 
        // step 4: go backwards to construct the best
        // chain & then reverse it
        //
        vector<Part> best_chain;
        for (int k = best_end; k != -1; k = best_from[k]){
            best_chain.emplace_back(_parts[k]);
        }

        std::reverse(best_chain.begin(), best_chain.end());

        //
        // step 5: check if this is a full chain
        //
        if (best_len[best_end] >= this->size - threshold){
            vector<Align> reconstruction;
            for (uint i = 0; i < best_chain.size(); i++){
                for (uint j = 0; j < this->_aligns.size(); j++){
                    if (this->_aligns[j].id == best_chain[i].id){
                        reconstruction.emplace_back(this->_aligns[i]);
                        break;
                    }
                }
            }
            this->_full_aligns.emplace_back(reconstruction);
        }

    }

    bool has_reconstruction(void){
        return !this->_full_aligns.empty();
    }

    vector<Align> & get_reconstruction(void){
        return this->_full_aligns.front();
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
            if (c != 'D'){
                cur_str = cigar.substr(last_p, p - last_p);
                size   += stoi(cur_str);
            }
            last_p = p + 1; // skip this char
            p++;
        }
    }

    return size;
}

bool
is_unmapped(const string &flag){
    //
    // check if this is an unmapped record
    //

    return stoi(flag) & 0x4;
}

bool
is_primary(const string &flag){
    //
    // check if the alignment is primary
    //

    int enc = stoi(flag); // encoding

    if (enc & 0x100 || enc & 0x800){
        return false; }

    return true;
}


bool
reverse_strand(const string &flag){
    //
    // determine if this alignment is on the reverse strand
    //

    return (stoi(flag) & 0x10) != 0;
}

int
collect_alignments(const string &bamfile, vector<Scaffold> &mappings){
    
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
                aln.pos       = static_cast<uint>(stoi(parts[3]));
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
                    Scaffold scaf(name);
                    aln.id = scaf.get_count();
                    scaf.add_align(aln);
                    if (is_primary(parts[1])){
                        scaf.size = get_sequence_length(aln.cigar);
                    }
                    mappings.push_back(scaf);
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
write_mappings(vector<Scaffold> &mappings){

    //
    // create the output stream to construct the file
    // containing all the possible full reconstructions
    //
    std::ofstream ofh("Separated.Contigs.txt");
    uint failed = 0, algn_count;

    ofh << "#Contig.ID\tAlignments\n";

    for (int i = 0; i < mappings.size(); i++){
        //
        // try to reconstruct the scaffold based on its
        // alignments
        //
        mappings[i].reconstruct();

        if (!mappings[i].has_reconstruction()){
            failed++;
            continue; // was too fragmented
        }

        vector<Align> &reconstruction = mappings[i].get_reconstruction();
        algn_count = reconstruction.size();

        ofh << mappings[i].name << '\t';

        for (int j = 0; j < reconstruction.size(); j++){
            Align &aln = reconstruction[j];
            char c = j == algn_count - 1 ? '\n' : ',';
            ofh << aln.ref_seq << '_' << aln.pos << '_' << aln.direction << c;
        }

    }

    ofh.close();

    cerr << "Number of scaffolds that failed to reconstruct " << failed << '\n';

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
    vector<Scaffold> mappings;
    collect_alignments(bamfile, mappings);

    //
    // step 2: write these out to a file
    //
    write_mappings(mappings);

    return 0;
}