# Updated 06-22-2025
# For counting errors in aligned modmapping experiments

# To use: python count_dels_pileup_1.2.py \
    # -n [bam file] \
    # -f [fasta reference] \
    # -o [output prefix (optional)] \
    # -t [threads (optional)] \
    # -w [window size (optional)] \
    # -d [minimum depth (optional)]

import pysam
from argparse import ArgumentParser
from multiprocessing import Pool, cpu_count
import functools
import csv
import time

#pileup function
def process_region(region, bam_path, fasta_path, min_depth):
    """
    Worker function to process a specific genomic window.
    region: tuple (ref_name, start, end)
    """
    ref, r_start, r_end = region
    
    bam = pysam.AlignmentFile(bam_path, "rb")
    ref_fasta = pysam.FastaFile(fasta_path)
    ref_seq = ref_fasta.fetch(ref, r_start, r_end)
    
    results = []
    
    # flag_filter=1540 = exclude unmapped (4) + QC-fail (512) + duplicate (1024)
    # Secondary alignments (256) are intentionally included so that multimapping
    # reads contribute to depth and deletion rates, consistent with IGV display.
    for pileupcolumn in bam.pileup(ref, r_start, r_end, truncate=True, flag_filter=1540):
        # Include indels first so depth counts reads with deletions at this position
        seqs_indels = pileupcolumn.get_query_sequences(mark_matches=True, add_indels=True)
        depth = len(seqs_indels)

        if depth < min_depth:
            continue

        current_pos = pileupcolumn.reference_pos
        #Load pileups to count bases (excludes reads with deletion at this position)
        seqs_bases = pileupcolumn.get_query_sequences(mark_matches=True)
        seq_str = "".join(seqs_bases).upper()

        Amut = seq_str.count('A')
        Cmut = seq_str.count('C')
        Tmut = seq_str.count('T')
        Gmut = seq_str.count('G')
        insert = seqs_indels.count('-')
        dels = seqs_indels.count('*')
        
        local_pos = current_pos - r_start
        base = ref_seq[local_pos].upper()
        
        results.append([
            ref, current_pos, current_pos + 1, base, depth, 
            dels, float(dels)/depth, insert, float(insert)/depth, 
            Amut, float(Amut)/depth, Cmut, float(Cmut)/depth, 
            Tmut, float(Tmut)/depth, Gmut, float(Gmut)/depth
        ])
    
    bam.close()
    ref_fasta.close()
    return results

#split into windows
def get_windows(bam_path, window_size=5_000_000):
    windows = []
    with pysam.AlignmentFile(bam_path, "rb") as bam:
        for ref, length in zip(bam.references, bam.lengths):
            for start in range(0, length, window_size):
                end = min(start + window_size, length)
                windows.append((ref, start, end))
    return windows


def main():
    parser = ArgumentParser(description="For counting errors in aligned modmapping experiments")
    parser.add_argument('-n','--nt', default="nt.bam", help="sorted bam file")
    parser.add_argument('-f', '--fasta', default="reference.fa", help="reference fasta file")
    parser.add_argument('-o','--out', default="error_rate", help="output prefix")
    parser.add_argument('-t','--threads', type=int, default=cpu_count(), help="number of threads to use")
    parser.add_argument('-w','--window', type=int, default=5_000_000, help="window size in bp")
    parser.add_argument('-d','--min_depth', type=int, default=1, help="minimum read depth")
    
    args = parser.parse_args()

    windows = get_windows(args.nt, window_size=args.window)
    total_windows = len(windows)
    print(f"Divided genome into {total_windows} windows.")
    print(f"Using {args.threads} threads. Filtering for depth >= {args.min_depth}...")

    worker = functools.partial(
        process_region, 
        bam_path=args.nt, 
        fasta_path=args.fasta, 
        min_depth=args.min_depth
    )

    headers = [
        'ref','start','end','base','depth','dels','del_rate','insert','ins_rate', 
        'Amut', 'A_rate', 'Cmut', 'C_rate', 'Tmut', 'T_rate', 'Gmut', 'G_rate'
    ]

    csv_file = f"{args.out}.csv"
    bed_file = f"{args.out}.bed"
    
    # Calculate exactly how many windows make up 5%
    five_percent_step = max(1, total_windows // 20)
    start_time = time.time()

    with open(csv_file, 'w', newline='') as f_csv, open(bed_file, 'w', newline='') as f_bed:
        csv_writer = csv.writer(f_csv)
        bed_writer = csv.writer(f_bed, delimiter='\t')
        
        csv_writer.writerow(headers)
        
        with Pool(processes=args.threads) as pool:
            for i, result in enumerate(pool.imap(worker, windows), 1):
                if result:  
                    csv_writer.writerows(result)
                    bed_writer.writerows(result)
                
                if i % five_percent_step == 0 or i == total_windows:
                    percent = (i / total_windows) * 100
                    elapsed_sec = time.time() - start_time
                    speed = i / elapsed_sec
                    
                    # flush=True forces the print to write to the log file immediately
                    print(f"Progress: {percent:02.0f}% | "
                          f"Windows: {i}/{total_windows} | "
                          f"Elapsed: {elapsed_sec/60:.1f} min | "
                          f"Speed: {speed:.1f} windows/sec", flush=True)

    print(f"\nDone! Results saved to {csv_file} and {bed_file}")

if __name__ == "__main__":
    main()