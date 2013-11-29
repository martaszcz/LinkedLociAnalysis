LinkedLociAnalysis
==================

General overview:
-----------------
A module to BioPython that perform analysis of linked loci and clusters them in blocks of associated SNP markers based on linkage disequilibrium (LD) (Hartl, Clark, 2009).
The analysis is performed using the Minority Allel Frequency (MAF) criterion and R^2 factor.

Module is proposed to be included to BioPython as Bio.PopGen.LinkedLociAnalysis.

Example usage:
--------------
chrLengths=[100]
filename='chr1_100'
minMAF = 0.05
minR2 = 1

chromosomes=LinkedLociAnalysis.splitAllels(filename, chrLengths)
selectedLociFromAllChromosomes=LinkedLociAnalysis.selectGenes(chromosomes,minMAF,minR2)
clusters=LinkedLociAnalysis.clusterChromosomes(selectedLociFromAllChromosomes)
print clusters

OR using convenience function:

clusters=LinkedLociAnalysis.clusterLinkedLoci(filename, chrLengths, minMAF=0.05, minR2=1)
print clusters

Example usage with more than one chromosome:
--------------------------------------------

chrLengths=[100,200,300]

filename='chr1_2_3_100_200_300'
minMAF = 0.05
minR2 = 1

chromosomes=LinkedLociAnalysis.splitAllels(filename, chrLengths)
selectedLociFromAllChromosomes=LinkedLociAnalysis.selectGenes(chromosomes,minMAF,minR2)
clusters=LinkedLociAnalysis.clusterChromosomes(selectedLociFromAllChromosomes)
print clusters

OR using convenience function:

clusters=LinkedLociAnalysis.clusterLinkedLoci(filename, chrLengths, minMAF=0.05, minR2=1)
print clusters

Example files and results:
--------------------------
Files including 1, 2 and 3 example chromosomes and results computed for this data are located in 'example_data' directory.
Above example data is generated with LDSO programme (Ytournel, 2008).
