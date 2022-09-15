
import numpy as np
from Bio import SeqIO

mRNA=SeqIO.read("NTRK1.gb","genbank")
aa_seq=mRNA.translate()
print(aa_seq)
SeqIO.write(aa_seq,"aa_seq.fasta","fasta")
