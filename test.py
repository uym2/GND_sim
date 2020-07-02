from sequence_lib import read_fasta, p_distance
from sys import argv

d=argv[1]

names_mutated,seqs_mutated = read_fasta(d+"/mutated.fasta")
names_org,seqs_org = read_fasta("AG-359-G18_contigs.fasta")

genome_org = {}
for name,seq in zip(names_org,seqs_org):
    genome_org[name] = seq

genome_mu = {}
for name,seq in zip(names_mutated,seqs_mutated):
    genome_mu[name] = seq

mutated = 0
for x in genome_org:
    m = p_distance(genome_org[x],genome_mu[x])
    mutated += m*len(genome_org[x])

total = sum(len(s) for s in seqs_org)
print(mutated/total)
