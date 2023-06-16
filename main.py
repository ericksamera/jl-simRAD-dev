from Bio import SeqIO

reader = SeqIO.parse('GCA_940337035.1_PGI_AGRIOTES_LIN_V1_genomic.fna', 'fasta')
for record in reader:
    pass