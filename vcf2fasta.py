import argparse
import pandas as pd
import pyfaidx
from sys import stdin as std_in
from sys import stdout as std_out

def is_valid_sequence(seq):
    """"
    Check if sequence is valid
    """
    return all(char in "AGTCN" for char in seq)

def get_variant_sequences(vcf_file, fasta_file, out_file, flank_size=500):
    """
    Given a VCF file and a refrence FASTA file, returns sequences containing the genomic sequence
    for each variant, with the variant allele inserted into the sequence at the appropriate
    position.
    
    Parameters:
    vcf_file (str): Path to the VCF file.
    fasta_file (str): Path to the FASTA file.
    out_file (str): Output FASTA file name.
    flank_size (int): Number of bases to include before and after the variant site.
    
    Returns:
    None
    """
    # Load the reference genome
    genome = pyfaidx.Fasta(fasta_file)
    
    # Initialize FASTA file to store the variant sequences
    fasta_file = open(f"{out_file}.fa", 'w')
    
    # Loop over each variant in the VCF file
    for index, row in pd.read_csv(vcf_file, sep='\t', comment='#', header=None, names=['chrom', 'pos', 'id', 'ref', 'alt', 'info']).iterrows():
        # Extract the chromosome, position, reference allele, and alternative allele
        chrom = row['chrom']
        pos = row['pos'] - 1
        ref = row['ref']
        alt = row['alt']
        variant_id = row['id']
        
        # Determine the size of the allele change
        ref_len = len(ref)
        alt_len = len(alt)
        change_len = alt_len - ref_len
        
        # Extract the left flanking sequence
        left_seq = genome[chrom][pos - flank_size : pos]
        
        # Create the variant sequence by replacing the reference allele with the alternative allele
        # Can only handle SNPs and simple indels for now
        if change_len == 0 and ref_len == 1 and alt_len == 1:
            # SNP
            ref_true = genome[chrom][pos : pos + 1]            
            right_seq = genome[chrom][pos + 1 : pos + 1 + flank_size - 1]            
            ref_seq = left_seq.seq.upper() + ref_true.seq.upper() + right_seq.seq.upper()
            alt_seq = left_seq.seq.upper() + alt.upper() + right_seq.seq.upper()
            
        elif change_len > 0 and ref_len == 1:
            # Insertion
            ref_true = genome[chrom][pos : pos + 1]
            right_seq = genome[chrom][pos + 1 : pos + 1 + flank_size - 1]
            ref_seq = left_seq.seq.upper() + ref_true.seq.upper() + right_seq.seq.upper()
            right_seq = genome[chrom][pos + alt_len : pos + alt_len + flank_size - alt_len]
            alt_seq = left_seq.seq.upper() + alt.upper() + right_seq.seq.upper()
            
        elif change_len < 0 and alt_len == 1:
            # Deletion
            ref_true = genome[chrom][pos : pos + ref_len]
            right_seq = genome[chrom][pos + ref_len : pos + ref_len + flank_size - ref_len]
            ref_seq = left_seq.seq.upper() + ref_true.seq.upper() + right_seq.seq.upper()
            right_seq = genome[chrom][pos + 1 : pos + 1 + flank_size - 1]
            alt_seq = left_seq.seq.upper() + alt.upper() + right_seq.seq.upper()
            
        else:
            raise ValueError("Invalid change length")
        
        # Save FASTA files
        if len(ref_seq) == flank_size*2 and len(alt_seq) == flank_size*2 and ref_true == ref and ref_seq != alt_seq and is_valid_sequence(ref_seq) and is_valid_sequence(alt_seq):
            
            # Make sequence header
            sequence_header_ref = f"{chrom}:{pos+1}_{ref_true.seq.upper()}_{variant_id}"
            fasta_file.write(">" + sequence_header_ref + "\n")
            fasta_file.write(ref_seq + "\n")
            sequence_header_alt = f"{chrom}:{pos+1}_{str(alt).upper()}_{variant_id}"
            fasta_file.write(">" + sequence_header_alt + "\n")
            fasta_file.write(alt_seq + "\n")
                
    # Close the reference genome file
    genome.close()
    
    # Close the output FASTA file
    fasta_file.close()


def get_options():

    parser = argparse.ArgumentParser(description='vcf2fasta')
    
    parser.add_argument('-I', metavar='input', nargs='?', default=std_in, required=True,
                        help='path to the input VCF file, defaults to standard in')
    parser.add_argument('-R', metavar='reference', required=True,
                        help='path to the reference genome fasta file')
    parser.add_argument('-O', metavar='output', nargs='?', default=std_out, required=True,
                        help='path to the output FASTA file, defaults to standard out')
    parser.add_argument('-F', metavar='flank', nargs='?', default=500,
                        type=int)
    args = parser.parse_args()

    return args

def main():

    # Get options
    args = get_options()
    
    # Get FASTAs
    get_variant_sequences(vcf_file=args.I, fasta_file=args.R, out_file=args.O, flank_size=args.F)
    
if __name__ == '__main__':
    main()    

