import numpy as np
from Bio import SeqIO

mRNA=SeqIO.read("NTRK1.gb","genbank")
print(mRNA)

'''
aa_seq=mRNA.translate(to_stop=True)

print(aa_seq)
SeqIO.write(aa_seq,"aa_seq.fasta","fasta")
'''
