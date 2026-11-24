
# Script to count events at NUN motif in NUN standard pools. Input is a .bam file that has been aligned to NUN oligo library

import pandas as pd
from collections import Counter
import numpy as np
import pysam
from Bio import SeqIO
from Bio.Seq import Seq
from argparse import ArgumentParser
from collections import defaultdict

def main():
    parser = ArgumentParser(description='Demultiplex reads using custom barcodes')
    parser.add_argument('-b','--bam', help = 'Input bam file')
    args = parser.parse_args()
    bam_filtered = NUN_cleanup(args.bam)
    counts = count_errors('bamfiltered.bam')
    counts.to_csv(args.bam+".count.csv")  
    
    
def NUN_cleanup(infile):
    #Function to remove any reads with barcodes that aren't a perfect match
    bamfile = pysam.AlignmentFile(infile, 'rb')
    bamfiltered = pysam.AlignmentFile("bamfiltered.bam", 'wb', template = bamfile)
    match = 0
    unmatch = 0
    #Iterate through all read
    for read in bamfile.fetch():
      #Get vector of reference base positions in read
      ref_pos = read.get_reference_positions(full_length = False)
      #Check if ends of barcode are present
      if 87 in ref_pos and 94 in ref_pos:
        #Check if barcode in read matches barcode in reference and add to output file
        bar_start = ref_pos.index(87)
        bar_end = ref_pos.index(94)
        if read.get_reference_sequence()[87:94] == read.query_sequence[bar_start:bar_end]:
          bamfiltered.write(read)
          match += 1
        else:
          unmatch += 1
      else: 
          unmatch += 1
    match = match*1.00
    unmatch = unmatch*1.00
    print(match)
    print(unmatch)
    print("Percent reads with correct barcode:"+str((match/(match+unmatch))*100)+"%")
    bamfiltered.close()
    pysam.index('bamfiltered.bam')
    

def count_errors(infile):
    # Command to categorize events at NUN position for each possible NT context
    infile = pysam.AlignmentFile(infile, 'rb')
    contexts = infile.references
    event = ['T','C','G','A','insert','del','total']
     # Initialize counts matrix
    counts = np.zeros((16, 7))
    pileup_bam = pysam.AlignmentFile("pileup.bam", 'wb', template = infile)
    #Iterator for counting reads
    i = 0
    for read in infile.fetch():
      i += 1
      pileup_bam.write(read)
      # infile.pileup() has a max read depth of 10,000 - to solve this make temp .bam file of every 5000 reads        and perform pileup on those - then sum results at the end
      if i % 5000 == 0:
        print('Reads processed: '+ str(i))
        pileup_bam.close()
        pysam.sort("-o", "pileup.sort.bam", "pileup.bam")
        pysam.index("pileup.sort.bam")
        pileup_bam_ind = pysam.AlignmentFile("pileup.sort.bam", 'rb')
        #Iterate through references and count error rates
        for nun in contexts:
          for pile in pileup_bam_ind.pileup(contig=nun, start = 80, stop = 81, truncate = True):
            nun_seq = pile.get_query_sequences(add_indels = True)
            pilecount = Counter(nun_seq)
            #Count indels
            indels = 0
            # To count indels, sum value of every event with 'N' (as indels are split by length, e.g. NN, NNN)
            for key in pilecount:
              if 'N' in key:
                indels += pilecount[key]
            counts[contexts.index(nun), 0] += pilecount['T']
            counts[contexts.index(nun), 1] += pilecount['C']
            counts[contexts.index(nun), 2] += pilecount['G']
            counts[contexts.index(nun), 3] += pilecount['A']
            counts[contexts.index(nun), 4] += indels
            counts[contexts.index(nun), 5] += pilecount['*']
            counts[contexts.index(nun), 6] += len(nun_seq)
            pileup_bam = pysam.AlignmentFile("pileup.bam", 'wb', template = infile)
    counts = pd.DataFrame(counts, columns =event, index = contexts)
    print counts
        
    return counts
        
      
      
        
    
    
      
    
    

# =========================================================
        
if __name__=="__main__":
    main()