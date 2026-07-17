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

    //
    // we need to define a constructor
    //
    Part(uint8_t id_, uint start_, uint length_) :
        id(id_), start(start_), length(length_) {}
};


class Scaffold {

private:
    vector<Align>          _aligns; // hold onto each alignment
    vector<vector<Align>> _full_aligns; // contains only alignments that can represent this full contig
    uint8_t               _count = 0;

    void _extend_path(const vector<Part> &parts, uint dist, vector<Part> &current,
                      uint current_len, vector<Part> &best_chain, uint &best_len){
        
        //
        // recursive loop to reconstruct all possible
        // paths
        //

        // is this the new best reconstruction?
        if (current_len > best_len){
            best_len   = current_len;
            best_chain = current;
        }  
        
        uint cur_end = current.back().start + current.back().length - 1;

        uint lower_bound = cur_end > dist ? cur_end - dist : 0;

        for (uint i = 0; i < parts.size(); i++){
            bool used = false;
            for (uint j = 0; j < current.size(); j++){
                // do not include itself
                if (current[j].id == parts[i].id){
                    used = true;
                    break;
                }
            }
            if (used){
                continue;
            }
            if (parts[i].start >= lower_bound && parts[i].start <= cur_end + dist){
                current.push_back(parts[i]);
                _extend_path(parts, dist, current, current_len + parts[i].length, best_chain, best_len);
                current.pop_back(); // backtrack
            }

            }
        }

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
        const uint threshold = 5000;
        // cerr << "Size of " << this->name << " " << this->size << '\n';

        for (int i = 0; i < this->_aligns.size(); i++){
            Align &aln          = this->_aligns[i];
            const string &cigar = aln.cigar;
            uint p = 0, n = cigar.size(), last_p = 0, size = 0, length = 0;
            string lstr; // string representation of the length
            char c;

            uint leading_clip = 0, aligned_len = 0;
            bool first_feature = true;

            while (p != n){
                c = cigar[p];
                if (isdigit(c)){
                    p++;
                }
                else {
                    lstr = cigar.substr(last_p, p - last_p);
                    size = stoi(lstr);
                
                    if (first_feature && (c == 'H' || c == 'S')){
                        leading_clip = size;
                    }
                    else if (c == 'M' || c == 'I'){
                        aligned_len += size;
                    }
                    first_feature = false;
                    last_p = p + 1;
                    p++;
                }
            }

            int start;
            if (aln.direction == 'R'){
                start = this->size - leading_clip - aligned_len;
            }
            else {
                start = leading_clip;
            }
            _parts.emplace_back(aln.id, start, aligned_len);
        }

        //
        // step 2: try to reconstruct using adjacent parts to each other
        // emphasizing the largest parts over smaller parts (parts == aligned portions)
        //
        const uint nparts = _parts.size();
        
        if (nparts == 1){
            if (_parts[0].length >= this->size - threshold){
                vector<Align> v;
                for (uint i = 0; i < this->_aligns.size(); i++){
                    if (this->_aligns[i].id == _parts[0].id){
                        v.push_back(this->_aligns[i]);
                        break;
                    }
                }
                this->_full_aligns.push_back(v);
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
        const uint   dist = 1000;
        uint     best_len = 0;
        vector<Part> best_chain;

        //
        // this will be a recursive function
        //
        for (uint i = 0; i < nparts; i++){
            vector<Part> current = {_parts[i]};
            _extend_path(_parts, dist, current, _parts[i].length, best_chain, best_len);
        }
   

        //
        // step 5: check if this is a full chain
        //
        if (best_len >= this->size - threshold){
            vector<Align> reconstruction;
            for (uint i = 0; i < best_chain.size(); i++){
                for (uint j = 0; j < this->_aligns.size(); j++){
                    if (this->_aligns[j].id == best_chain[i].id){
                        reconstruction.push_back(this->_aligns[j]);
                        break;
                    }
                }
            }
            this->_full_aligns.push_back(reconstruction);
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

    cerr << failed << " out of " << mappings.size() << " scaffolds failed to reconstruct\n";

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