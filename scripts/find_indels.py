#!/bin/env python3

import argparse
import os
import sys
import sys
import gzip

class Indel_Record:
    def __init__(self, type_: str, pos: int, size: int, seq: str):
        self.type = type_
        self.pos  = pos if (type_ == 'I') else pos - size
        self.size = size
        self.end  = pos
        self.seq  = seq

def get_arguments():
    """function to return arguments to find_indels.py"""

    parser = argparse.ArgumentParser(description="Find indels from aligned reads")
    parser.add_argument("-b", "--bam", required=True, type=str, help="path to bam file")
    parser.add_argument("-r", "--ref", required=True, type=str, help="path to reference genome")
    parser.add_argument("-d", "--deletion_length", type = int, default=1, help = "Minimum length to report a deletion [default = 1]")
    parser.add_argument("-i", "--insertion_length", type = int, default=1, help = "Minimum length to report an insertion [default = 1]")
    parser.add_argument("-w", "--write", action="store_true", help="Write reported sequences to a fasta file")

    args = parser.parse_args()
    # check bam
    bam_file = args.bam
    assert os.path.isfile(bam_file), f"Could not locate {bam_file}"

    # check index file
    bai_file1 = f"{bam_file}.bai"
    bai_file2 = bam_file.replace(".bam", ".bai")
    if (os.path.isfile(bai_file1) == False) and (os.path.isfile(bai_file2) == False):
        msg = f"Could not find corresponding *.bai file for {os.path.basename(bam_file)}"
        sys.exit(msg)

    # check that samtools is installed
    which_sam = os.popen("which samtools").read()
    if (which_sam.startswith("which: no samtools")):
        msg = "Samtools was not found in path"
        sys.exit(msg)

    # check integar params
    del_len = args.deletion_length
    ins_len = args.insertion_length
    assert del_len >= 0, "--deletion_length cannot be less than 0"
    assert ins_len >= 0, "--insertion_length cannot be less than 0"

    # check for reference genome
    ref_fa = args.ref
    assert os.path.isfile(ref_fa), f"Could not located {ref_fa}"

    return bam_file, ref_fa, del_len, ins_len, args.write

def make_output_name(bam_file: str):
    """creates the output name"""

    bname   = os.path.basename(bam_file)
    outname = bname.replace(".bam", "_indel_sites.tsv")

    return outname

def read_ref(ref_fa: str):
    """read a fasta file into memory"""

    # variables
    ref_seqs = {} # will return a dict
    seq      = ''
    seq_id   = ''
    fh       = gzip.open(ref_fa, "rt") if ref_fa.endswith(".gz") else open(ref_fa, 'r')

    for line in fh:
        # skip empty or commented lines
        if (len(line) == 0) or (line[0] == '#'):
            continue
        line = line.strip()
        if (line[0] == '>'):
            # for the 1st record
            if (seq == ''):
                seq_id = line[1:]
                continue
            else:
                ref_seqs[seq_id] = seq 
                seq = ''
                seq_id = line[1:]
        else:
            seq += line.upper()

    # for the last record
    if (len(seq) > 0):
        ref_seqs[seq_id] = seq

    return ref_seqs

def rev_seq(seq: str):
    """reverse complement a sequence"""

    seq_rc = ''

    rc_dict = {'A':'T',
               'T':'A',
               'C':'G', 
               'G':'C',
               'N':'N',
               '-':'-',
               'K':'M', # G||T : C||A
               'M':'K', # C||A : G||T
               'W':'W', # A||T : T||A
               'S':'S', # G||C : C||G
               'R':'Y', # A||G : T||C
               'Y':'R', # T||C : A||G
               'B':'V', # C||G||T : A||C||G
               'V':'B', # A||C||G : C||G||T
               'D':'H', # A||G||T : T||C||A
               'H':'D'}

    for nuc in seq[::-1]:
        seq_rc += rc_dict[nuc.upper()]

    return seq_rc

def process_cigar(cigar: str, del_len: int, ins_len: int, pos: int, ref_seq: str, seq: str, rev: bool):

    indel_recs = list() # store records of each indel

    # pointers
    i,j = 0, 0

    # set of operators that consume query and reference
    cq = {'M', 'I', 'S', '=', 'X', 'H'}
    cr = {'M', 'D', 'N', '=', 'X'}

    # keep track of where we are in the sequence and reference
    ref_pos = pos - 1 # first base is inclusive
    seq_pos = 0

    while j < len(cigar):
        if (cigar[j].isdigit()):
            j += 1
        else:
            op   = cigar[j]
            size = int(cigar[i:j])
            if (op in cr):
                ref_pos += size
            if (op in cq):
                seq_pos += size
            if ((op == 'I') and (size >= ins_len)) or \
               ((op == 'D') and (size >= del_len)):
                iseq = seq[seq_pos - size:seq_pos] if (op == 'I') else ref_seq[ref_pos - size:ref_pos]
                rec  = Indel_Record(op, ref_pos, size, iseq)
                indel_recs.append(rec)
            j += 1
            i = j
    
    return indel_recs

def process_bam(bam_file: str, del_len: int, ins_len: int, write_fasta: bool, ref_seqs: dict):
    """process the bam file"""

    # create io streams
    cmd     = f"samtools view {bam_file}"
    fstream = os.popen(cmd, 'r')
    outname = make_output_name(bam_file)
    ostream = None
    fa      = outname.replace(".tsv", ".fa") if write_fasta else None
    sstream = None # for fasta sequences

    for rec in fstream:
        if (rec[0] == '@'):
            continue
        fields = rec.split('\t')
        flag   = int(fields[1])
        # only considering primary alignments
        if ((flag & 256) == 256) or ((flag & 2048) == 2048):
            continue
        rev      = True if ((flag & 16) == 16) else False
        seq_name = fields[0]
        ref      = fields[2]
        cigar    = fields[5]
        pos      = int(fields[3])
        seq      = fields[9]
        ind_recs = process_cigar(cigar, del_len, ins_len, pos, ref_seqs[ref], seq, rev)
        if (len(ind_recs) > 0):
            for rec in ind_recs:
                outline = f"{seq_name}\t{ref}\t{rec.pos}\t{rec.end}\t" + \
                          f"{rec.size}\t{rec.type}\t{rec.seq}\n"
                try:
                    ostream.write(outline)
                except:
                    ostream = open(outname, 'w')
                    ostream.write(outline)
            if (write_fasta):
                # check if need to rev complement
                if (rev):
                    seq = rev_seq(seq)
                outline = f">{seq_name}\n{seq}\n"
                try:
                    sstream.write(outline)
                except:
                    sstream = open(fa, 'w')
                    sstream.write(outline)
        
    # close io streams
    fstream.close()
    if (ostream != None):
        ostream.close()
    if (sstream != None):
        sstream.close()
                    
def main():
    """entry point to the pipeline"""

    # get the inputs
    bam_file, ref_fa, del_len, ins_len, write_fasta = get_arguments()

    # read all the sequences into memory # TODO read only seq in bam header
    ref_seqs = read_ref(ref_fa)

    # start processing
    process_bam(bam_file, del_len, ins_len, write_fasta, ref_seqs)

if __name__ == "__main__":
    main()

