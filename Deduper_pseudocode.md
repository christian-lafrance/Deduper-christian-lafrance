Christian La France  
Deduper PCR duplication removal pseudocode  
Bi624     


### 1. Define the Problem
----
PCR duplicates can occur when PCR preferentially amplifies some sequences more than others. This can occur for many reasons including lower GC content or shorter sequences. This biases read counts in favor of sequences that amplify more efficiencly. This becomes an issue for any type of count analysis, such as differential expression analysis. It can be difficult to tell which sequences are more abundant because those genes are upregulated and which were artifacts of increased PCR efficiency. For this reason it is necessary to remove PCR duplications prior to downstream analysis.   

PCR duplicates will be removed after reference alignment. It will take in a SAM file and use the chromosome number, bitwise flag, start position, and QNAME to determine if a read is a PCR duplicate. A read will be considered a PCR duplicate if it comes from the same strand, has the same UMI, chromosome, and start position as another read. 

### 2. Pseudocode  
----  
Unix:  
The first step will be grep out the lines that start with @. Then sort the SAM file by chromosome (RNAME, field 3), then, maintaining chromosome sort order, sub-sort by pos (1-based leftmost mapping position, field 4). Position will be adjusted based on the CIGAR string but this will hopefully get them close enough to start with. 

Python script:  
- Read in the UMI's and store them in a set. 
- Create an empty dictionary which will be used to store unique SAM entries that occur at the same position/chromosome. 
- Create/open an output SAM file to write to. 
- Open the input SAM file and iterate through line by line.  
    - Check the end of the QNAME for a valid UMI. If it is not valid, skip it and continue to the next read.  
        - If it is valid, using the functions defined below, generate a string with chromosome number, adjusted-start position, strand (+ or -), and UMI sequence. 
            For example, this would look like:
            chr9.123456.+.ATGCATGCA
            This string will be used as the key and the value will be a count.  

        - Check the dictionary for this string. If it is not in the dictionary, add the string to the dictionary with an initial count of 1, and write the current line to the output SAM file. 
        - If the string is already in the dictionary, increment the counter but do not write the line to the output file. Go to the next line in the SAM file.   

Note: need to develop a strategy to prevent dictionary from getting too big. Maybe clear it after each chromosome?

### 3. Functions  
----  

```
def get_chrom(sam_line: str) -> str:
    '''
    Takes the current line of the SAM file as an argument and returns the
    chromosome number (second field, RNAME).
    '''
    return chrom_num
```
Test input: NS500451:154:HWKTMBGXX:1:11101:24260:1121:CTGTTCAC	31	2	100	63	1S3M2D4N1S
Expected output: 2

```
def get_pos(sam_line: str) -> int:
    '''
    Takes the current line of the SAM file as an argument, parses the CIGAR 
    string and returns the adjusted 1-based left-most start position of the 
    mapped read. Will also need to call the get_strand function. 
    For (+) strand: pos = mapped_pos + S (beginning of CIGAR)
    For (-) strand: pos = mapped_pos + M's + D's + N's + S(end of CIGAR) - 1
    '''
    return adjusted_start_pos
```
Test input: NS500451:154:HWKTMBGXX:1:11101:24260:1121:CTGTTCAC	31	2	100	63	1S3M2D4N1S
Expected output: 109

```
def get_UMI(sam_line: str) -> str:
    '''
    Takes the current line of the SAM file as an argument and returns the 
    UMI from the end of the QNAME. 
    '''
    return umi
```
Test input: NS500451:154:HWKTMBGXX:1:11101:24260:1121:CTGTTCAC	31	2	100	63	1S3M2D4N1S
Expected output: CTGTTCAC

```
def get_strand(sam_line: str) -> bool:
    '''
    Takes the current line of the SAM file as an argument and returns a
    boolean value representing if the current read is the (+) strand. 
    True means it is the (+) strand and False means it is the (-) strand. 
    Check the bitwise flag (field 2) and if bit 16 is true then the current 
    read is from the (-) strand. 
    '''
    return is_plus_strand
```
Test input: NS500451:154:HWKTMBGXX:1:11101:24260:1121:CTGTTCAC	31	2	100	63	1S3M2D4N1S
Expected output: False