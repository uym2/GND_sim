# GND_sim
  * Used a gamma distribution (with an alpha parameter) to draw the relative rate for each gene. 
  * Given a desired GND, compute the corresponding number of DNA mutations needed. The number of mutations in aa, nmus, is roughly 2/3 that of DNA. Empirically, 5/8 works better than 2/3, so this simulation uses 5/8 instead.
  * Randomly select nmus amino acids (i.e. sampling without replacement) to mutate using BLOSUM62 model; each amino acid is selected with a probability determined by the rate of the gene it belongs to. 
  * Avoid adding/removing start/end codons: don't allow the start and end condons to mutate
  * Avoid interrupting the reading frame: when encounter a pair of genes that overlap each other (possibly with different reading frames), don't mutate the overlapping region.
  * Back translation: when there are multiple codons for an amino acid, choose the one that is closest to that of the original genome (i.e. minize the number of mutations).
  * Intergenic regions (~5% of the genome): We don't mutate them.
