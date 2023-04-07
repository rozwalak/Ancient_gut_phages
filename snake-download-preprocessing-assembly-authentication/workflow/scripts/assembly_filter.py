import sys
from Bio import SeqIO

input_file = sys.argv[1]
output_file = sys.argv[2]
min_length = 4000
min_coverage = 20

with open(output_file, "w") as output_handle:
    for record in SeqIO.parse(input_file, "fasta"):
        id_parts = record.id.split("_")
        print(id_parts)
        length = int(id_parts[3])
        cov = float(id_parts[5])
        print(length, cov)
        if length >= min_length and cov >= min_coverage:
            SeqIO.write(record, output_handle, "fasta")