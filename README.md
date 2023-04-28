# VCF2FASTA
Given a VCF file and a reference FASTA file, returns two sequences for each variant, with the reference and variant (alternative) allele inserted into the sequence at the appropriate position.

Usage:
```
python vcf2fasta.py -I input.vcf -R reference.fasta -O output.fa -F 500
```
