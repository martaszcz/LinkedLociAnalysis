#!/usr/bin/python
# -*- coding: utf-8  -*-

def splitAllels(filename, chrLengths):
	"""
	Return a list of all loci observations divided into list containing separate chromosomes loci.
	Example return value for 3 chromosomes with 5 loci, 3 loci, 10 loci respectively, 3 entities (ie. 3 people) each:
	[
	[[1,2,2,2,1,2,1,2,2,1], [1,2,2,2,1,2,1,2,2,1], [1,2,2,2,1,2,1,2,2,1]],
	[[1,1,2,2,2,1],[1,1,2,2,2,1],[1,1,2,2,2,1]],
	[[1,2,2,2,1,2,1,2,2,1,1,2,2,2,1,2,1,2,2,1],[1,2,2,2,1,2,1,2,2,1,1,2,2,2,1,2,1,2,2,1],[1,2,2,2,1,2,1,2,2,1,1,2,2,2,1,2,1,2,2,1]]
	]
	"""
	file=open(filename, "r")
	lines=file.readlines()
	chromosomes=[[]]
	chrStartPositions=[1]
	for i in range(1,len(chrLengths)):
		chromosomes.append([])
		chrStartPositions.append(chrStartPositions[-1]+chrLengths[i-1]*2)
	chrStartPositions.append(chrStartPositions[-1]+chrLengths[-1]*2)
	for line in lines:
		line=line.strip().split() 
		intLine = [ int(i) for i in line]
		line=intLine
		for chr in range(len(chromosomes)):
			chromosomes[chr].append(line[chrStartPositions[chr]:(chrStartPositions[chr+1])])
	return chromosomes 



def selectGenes(chromosomes, minMAF, minR2):
	"""
	Returns set of edges between loci which have MAF greater than minMAF and R^2 coefficient on these edges greater than minR2. Results from each chromosome are in a separate list. All these lists are returned in one, bigger list.
	Example returned value for 3 chromosomes:
	[
	[{}] # edges in first chromosome
	# edges in second chromosome
	# edges in third chromosome
	]
	"""
	import collections
	linkedLociFromAllChromosomes=[]
	for chr in chromosomes:
		edges=collections.defaultdict(set)
		locusIndices=[]
		for c in range(0, len(chr[0]), 2):
			total=0
			for r in range(len(chr)): 
				total+=chr[r][c]+chr[r][c+1]
			maf=2-(float(total)/(len(chr)*2))
			if maf>=minMAF and maf<=1-minMAF: 
				locusIndices.append((c,maf))
		linkedLoci=[]
		for a in range(len(locusIndices)):
			locAind=locusIndices[a][0]
			locAmaf=locusIndices[a][1]
			aMAF2=(1-locAmaf) * locAmaf
			for b in range(a+1, len(locusIndices)):
				locBind=locusIndices[b][0]
				locBmaf=locusIndices[b][1]
				bMAF2=(1-locBmaf) * locBmaf
				a1b1=0.0
					a2b1=0.0
				a1b2=0.0
				a2b2=0.0
				for r in range(len(chr)):
					for x in (0,1):
						for y in (0,1):
							if chr[r][locAind+x] == 1:
								if chr[r][locBind+y] == 1:
									a1b1+=1
								else:
									a1b2+=1
							else:
								if chr[r][locBind+y] == 1:
									a2b1+=1
								else:
									a2b2+=1
				d=(a1b1 * a2b2 - a1b2 * a2b1)/(len(chr)*len(chr))
				r2=(d*d)/(aMAF2 * bMAF2)
				if r2 > minR2:
					edges[locAind].add(locBind)
					edges[locBind].add(locAind)
		linkedLociFromAllChromosomes.append(edges)
	return linkedLociFromAllChromosomes


def cluster(edges):
	"""
	Returns list of lists, each containing loci (of one chromosome) located in the same cluster.
	Unclustered loci are not included in the returned value.
	Example returned value:
	[ #list of 3 blocks found in clustered chromosome
	[1,15,16,17] # block of 4 loci
	[2,3,4,5,6,7,8] #block of 7 loci
	[9,11,12] # block of 3 loci
	]
	"""
	lociStillToCluster=set(edges.keys())
	clusters=[]
	while lociStillToCluster:
		neighbourPairsScore={}
		bestSeed=None
		bestSeedScore=-1
		for locus in lociStillToCluster:
			for neighbour in edges[locus].intersection(lociStillToCluster):
				tmpScore=len(edges[locus].intersection(edges[neighbour]))
				if tmpScore > bestSeedScore:
					bestSeedScore = tmpScore
					bestSeed = set((locus,neighbour))
		if bestSeedScore==-1:
			break
		newCluster=set()
		seed=bestSeed
		seed2 = set([loc/2 for loc in seed])
		newCluster=newCluster.union(seed2)
		lociStillToCluster=lociStillToCluster.difference(seed)
		neighbourSet=(edges[tuple(seed)[0]].intersection(edges[tuple(seed)[1]]).intersection(lociStillToCluster))
		while neighbourSet:
			neighbourScores={}
			for neighbour in neighbourSet:
				neighbourScores[neighbour]=len(neighbourSet.intersection(edges[neighbour]))
			newMax=max(neighbourScores,key=lambda x: neighbourScores[x])
			newCluster.add(newMax)
			lociStillToCluster.remove(newMax)
			neighbourSet=neighbourSet.intersection(edges[newMax])
		clusters.append(newCluster)
	return clusters

def clusterChromosomes(chromosomes):
	"""
	Convenience function that allows to cluster many chromosomes at once. Returned value is a list of lists, each containing lists of clusters.
	"""
	listOfClusters=[]
	for chr in chromosomes:
		listOfClusters.append(cluster(chr))
	return listOfClusters

def clusterLinkedLoci(filename, chrLengths, minMAF=0.05, minR2=1):
	"""
	Convenience function that allows the automation of whole process of clustering.
	Input: name of file with observated loci versions (filename), list or tuple of lengths of each chromosome (chrLengths), minimal MAF and R^2 used to select loci to cluster (minMAF, default 0.05; minR2, default 1).
	Lower minMAF make the process longer, but give more sensitivity.
	Lower minR2 values makes the process longer, but gives bigger clusters with lower linkage between loci.
	Output: the same as of function clusters;
	List of lists, each containing loci (of one chromosome) located in the same cluster.
	Unclustered loci are not included in the returned value.
	Example returned value:
	[ #list of 3 blocks found in clustered chromosome
	[1,15,16,17] # block of 4 loci
	[2,3,4,5,6,7,8] #block of 7 loci
	[9,11,12] # block of 3 loci
	]
	"""
	chromosomes=splitAllels(filename, chrLengths)
	selectedLociFromAllChromosomes=selectGenes(chromosomes,minMAF,minR2)
	clusters=clusterChromosomes(selectedLociFromAllChromosomes)
	return clusters
