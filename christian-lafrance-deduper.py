#!/usr/bin/env python


import argparse
import re


def get_args():
    parser = argparse.ArgumentParser(description='''
    This script reads in a SAM file containing single-end reads and outputs another SAM file with 
    PCR duplicates removed. It assumes that two reads are duplciates if they are on the same
    strand and chromosome, have the same left-most start position, and have the same UMI. 
    The input SAM file does not need to be sorted. This script stores strings of all UMIs, 
    chromosomes, strands, and positions seen and stores them in a trie data structure to use less memory. 
    Assumes a stranded library prep kit was used. 
    ''')
    parser.add_argument("-f", help="Designates absolute file path to input sam file", required=True, type=str)
    parser.add_argument("-o", help="Designates name and file path of sam file with PCR duplicates removed.", required=True, type=str)
    parser.add_argument("-u", help="Designates file containing the list of expected UMIs", required=True, type=str)
    return parser.parse_args()
args = get_args()


def read_UMIs(umi_file: str) -> set:
    '''
    Takes filename string as an argument, opens the file, and returns 
    a list of the UMI's condained in the file. 
    '''
    with open(umi_file, "r") as f:
        umis = f.read().split("\n")

    return umis


def get_UMI(qname:str) -> str:
    '''
    Takes in a QNAME from SAM header and returns the UMI from the end of the QNAME. 
    '''
    return qname.split(":")[-1]


def get_strand(flag: str) -> bool:
    '''
    Takes the the bitwise flag of the current line of the SAM file as an argument 
    and returns a boolean value representing if the current read is the (+) strand. 
    True means it is the (+) strand and False means it is the (-) strand. 
    Check the bitwise flag (field 2) and if bit 16 is true then the current 
    read is from the (-) strand. 
    '''
    if ((flag & 16) == 16): 
         return False #(it is the reverse complement of the genomic seq). #It is the (-) strand. 
    else:
        return True


def get_pos(pos: int, is_plus_strand: bool, cigar: str) -> int:
    '''
    Takes position, strand, and CIGAR string as arguments, parses the CIGAR 
    string and returns the adjusted 1-based left-most start position of the 
    mapped read. Returns the adjusted start positin as an int. 
    For (+) strand: pos = mapped_pos + S (beginning of CIGAR)
    For (-) strand: pos = mapped_pos + M's + D's + N's + S(end of CIGAR) - 1
    '''
    if is_plus_strand:
        num = ""
        for c in cigar:
            if c == "S":
                sc = int(num)
                return pos - sc
            elif c in {"M", "I", "D", "N", "H", "P", "=", "X"}:
                return pos
            else:
                num += c
        
    else:
        if cigar[-1] == "S":
            s = re.findall('[0-9]+(?=S)', cigar)
            pos += int(s[-1])
        m = sum([int(i) for i in re.findall('[0-9]+(?=M)', cigar)])
        d = sum([int(i) for i in re.findall('[0-9]+(?=D)', cigar)])
        n = sum([int(i) for i in re.findall('[0-9]+(?=N)', cigar)])
        return pos + m + d + n + - 1           


def check_trie(trie: dict, tup: tuple) -> bool:
    '''
    Generates a trie data structure as a nested dictionary to prevent
    redundant prefixes from using more memory, such as reads that have
    the same UMI, chromosome, and strand but have a different position. 
    Returns a boolean representing if the umi-chrom-strand-position combination 
    is unique. 
    '''

    if tup[0] not in trie:
        trie[tup[0]] = {}
        current = trie[tup[0]]
        for i in tup[1:]:
            current[i] = {}
            current = current[i]
        return True
    else:
        current = trie[tup[0]]
        for i in tup[1:]:
            if i not in current:
                current[i] = {}
                if type(i) == int:
                    return True # the current position has not been seen, write to file
            else:
                if type(i) == int: # this might cause problems???
                    return False # this read has been seen before
            current = current[i]


def check_dict(data: dict, value: int|str) -> dict:
    '''
    Checks if a value is in a dictionary and if so, increments
    the counter by 1. If not, it will add it to the dictionary and
    initialize counter to 1. Returns the updated dictionary.
    '''
    if value in data:
        data[value] += 1
    else:
        data[value] = 1

    return data


##### FUNCTION CALLS #####
valid_UMIs = read_UMIs(args.u)
trie = {}
counts = {}

'''
Open input and output same files. 
'''
with open(args.o, "w") as sam_out:
    with open(args.f, "r") as sam_in:
        '''
        Need to transfer header from sam_in to sam_out before iteration. 
        '''
        for line in sam_in:
            if line[0] == "@":
                sam_out.write(line)
            else:
                line_contents = line.split("\t")

                umi = get_UMI(line_contents[0])
                if umi in valid_UMIs:
                    strand = get_strand(int(line_contents[1]))
                    chrom = line_contents[2]
                    cigar = line_contents[5]
                    pos = get_pos(int(line_contents[3]), strand, cigar)


                    is_unique_read = check_trie(trie, (umi, chrom, strand, pos))

                    if is_unique_read == True:
                        sam_out.write(line)
                        counts = check_dict(counts, chrom)

                    else:
                        counts = check_dict(counts, "dups_removed")

                else:
                    counts = check_dict(counts, "invalid_umis")


# This will write a tsv file with all stats per chromosome. 
# with open("deduper_stats.tsv", "w") as stats:
#     for k, v in sorted(counts.items()):
#         stats.write(f"{k}\t{v}\n")