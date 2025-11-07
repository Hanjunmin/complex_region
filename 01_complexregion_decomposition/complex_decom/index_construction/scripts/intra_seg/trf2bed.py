import pandas as pd
import sys
import argparse
import glob

def trf_to_bed(input_files, output_file):
    trf = []
    header = '#chr start end PeriodSize CopyNumber ConsensusSize PercentMatches PercentIndels Score A C G T Entropy Motif Sequence'.split()
    chrom = None
    with open(input_files, 'r') as dat:
        for line in dat:
            splitline = line.split()
            
            if line.startswith("Sequence:"):
                chrom = int(line.split()[1].strip())
            elif line.startswith("@"):
                chrom = splitline[0][1:].strip()
            else:
                try:
                    try:
                        int(splitline[0])
                    except ValueError:
                        continue
                    trf.append([chrom] + splitline[0:(len(header)-1)])
                except IndexError:
                    pass
    
    trf_df = pd.DataFrame(trf, columns=header)
    trf_df["start"] = trf_df["start"].astype(int) - 1
    trf_df.sort_values(by=["#chr", "start"], inplace=True)
    trf_df.to_csv(output_file, sep="\t", index=False)
    return trf_df

def main():
    parser = argparse.ArgumentParser(description='TRF TO BED FILE')
    parser.add_argument('-i', '--input', required=True,
                       help='input')
    parser.add_argument('-o', '--output', required=True,
                       help='output')
    
    args = parser.parse_args()

    trf_df = trf_to_bed(args.input, args.output)
   
if __name__ == "__main__":
    main()