#! /usr/bin/env python

from sequence_lib import read_fasta, write_fasta
from random import *
from numpy.random import gamma, binomial, choice

def read_codon_table():
    codon = {}
    with open("codon_table.txt",'r') as fin:
        for line in fin:
            code = line.strip().split()
            codon[code[0]] = code[1:]
    return codon

def read_blosum62():
    M = {}
    with open("blosum62sym.csv") as fin:
        aa = [x[1] for x in fin.readline().strip().split(",")]
        for line in fin:
            row = line.strip().split(",")
            a1 = row[0][1]
            M[a1] = {}
            for a2,p in zip(aa[1:],row[1:]):
                M[a1][a2] = float(p) if a2 != a1 else 0
    return M

def normalize(v):
    s = sum(v)
    return [x/s for x in v]

def randomize_rates(n_gene,alpha):    
    rates = gamma(alpha,1/alpha,n_gene)
    return normalize(rates)

def evolve(gseqs,g_weights,n_mus,M):
    # n_mus: number of mutations

    print("Computing population ...")
    population = [(i,j) for i in range(len(gseqs)) for j in range(len(gseqs[i]))]

    print("Computing weights ...")
    weights = normalize([g_weights[i][j] for i in range(len(g_weights)) for j in range(len(g_weights[i]))])

    print("Choosing mutation sites ...")
    mutated = choice(range(len(weights)),size=n_mus,replace=False,p=weights)

    for k in mutated:
        i,j = population[k]
        seq = gseqs[i]
        c = choice(list(M[seq[j]].keys()),p=normalize(M[seq[j]].values()))
        gseqs[i] = seq[:j] + c + seq[j+1:]

def extract_locations(gnames):
# parse the gene names to extract the start and end points
    g_locations = []    
    #g_locations_adjusted = []
    #i_locations = []

    for i,g in enumerate(gnames):
        node,start,end,direction = g.strip().split("#")[:4]
        node = "_".join(node.strip().split("_")[:-1])
        start = int(start)-1
        end = int(end)-1
        rev = float(direction) < 0
        g_locations.append((node,start,end,rev))
    
    return g_locations

def assign_weights(gseqs,g_locations,alpha):
# randomly assign a (relative) mutation rate for each gene
# the rate multipliers are drawn from a gamma distribution
# all the loci in each gene are assigned the rate multipliers 
# of that gene as their weights; instead for the start and 
# end codons and the overlapping regions between the genes 
# are assigned weights 0 (i.e. not allowed to mutate)
# alpha control the variance of the rate (i.e. shape of the
# gamma distribution) 
    n_genes = len(gseqs)
    g_rates = randomize_rates(n_genes,alpha)
    g_weights = []
    prv_node, prv_start, prv_end = (None,None,None)

    for seq,rate,(node,start,end,_) in zip(gseqs,g_rates,g_locations):
        w = [rate]*len(seq)
        # weight 0 for start and end codons
        w[0] = 0 
        w[-1] = 0
        # check for overlapping with previous genes
        # here we assume that the genes are ordered
        # and a gene can only overlap with at most
        # one other gene
        if prv_node == node:
            if prv_end > start:
                prv_w = g_weights[-1]
                # set zero-weights for prv_w
                e = prv_end
                i = -1
                while e > start:
                    prv_w[i] = 0
                    e -= 3
                    i -= 1
                # set zero-weights for w    
                s = start
                i = 1
                while s < prv_end:
                    w[i] = 0
                    s += 3
                    i += 1                
        # add w to g_weights
        g_weights.append(w)    

    return g_weights

def reverse_complement(seq):
# compute the DNA reverse complement
    complement = {'A':'T','T':'A','G':'C','C':'G'}
    rseq = ''
    for x in seq[::-1]:
        rseq += complement[x]
    return rseq

def back_translate(peptide,ref,rev):
# back translation is not unique
# here we do the following:
# if a codon exactly matches the original, 
# then use it. Otherwise, randomly choose a codon 
# with the weight is the probability 
# that the reference dna mutated into that codon 
# eg: prob(ACT --> ACG) = 1/4*(3/4)**2; prob(ACT --> AAG) = (1/4)**2*(3/4)
    if rev:
        ref = reverse_complement(ref)
    codon = read_codon_table()
    dna = ''
    for i,aa in enumerate(peptide):
        ref_i = ref[3*i:3*i+3]
        w_total = 0
        translated = None
        for tr in codon[aa]:
            d = sum(a != b for a,b in zip(tr,ref_i))
            if d == 0:
                translated = tr
                break
            #w = (0.25**d)*(0.75**(3-d))
            w = 4-d
            w_total += w
            if random() <= w/w_total: # choose this codon with probability w/w_total
                translated = tr
        dna += translated
    D = sum(x!=y for x,y in zip(dna,ref))/len(dna)
    if rev:
        dna = reverse_complement(dna)
    return dna,D 

def stitch_back(names,seqs,g_locations,gseqs):
    # hash genome names
    genome = {}
    for node,seq in zip(names,seqs):
        genome[node] = seq

    # stitch the genes to the genome
    for gseq,(node,start,end,rev) in zip(gseqs,g_locations):
        ref = genome[node][start:end+1]
        dna,D = back_translate(gseq,genome[node][start:end+1],rev)
        new_seq = genome[node][:start] + dna + genome[node][end+1:]
        genome[node] = new_seq

    return list(genome.keys()),list(genome.values())

def compute_n_mus(n_sites,p,p_g=0.95):
# n_sites is the number of DNA sites, while
# n_mus is the number of aa mutations
# p is the proportion of DNA mutations in the genome
# p_g is the proportion of genes to the full genome
    nt2aa = 1/2
    return round(n_sites*p*nt2aa)

from sys import argv
from os import mkdir,getcwd,rmdir,listdir

p = float(argv[1]) # p is the proportion of aa mutations (i.e. 1-AAI)
alpha = float(argv[2])

M = read_blosum62()
names,seqs = read_fasta("AG-359-G18_contigs.fasta") # full genome
gnames,gseqs = read_fasta("AG-359-G18_contigs_genes.faa") # genes
g_locations = extract_locations(gnames)
g_weights = assign_weights(gseqs,g_locations,alpha)

#n_sites = sum(len(s) for s in seqs) 
#n_mus = compute_n_mus(n_sites,p)
n_sites = sum(len(s) for s in gseqs) 
n_mus = round(n_sites*p)

evolve(gseqs,g_weights,n_mus,M)
mutated_names, mutated_seqs = stitch_back(names,seqs,g_locations,gseqs)

outdir = "mutated_a" + str(int(alpha)) + "_aad" + str(round(p*100)).rjust(3,'0')
mkdir(outdir)
write_fasta(outdir+"/mutated_genes.faa",gnames,gseqs)
write_fasta(outdir+"/mutated.fasta",mutated_names,mutated_seqs)
