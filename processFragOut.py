#! /usr/bin/python
import sys
from scripts.genomeUtils import *

if len(sys.argv) < 4:
	sys.exit("Usage: processFragOut.py <fragout> <fasta> <outPrefix>")

fastadict = readFASTAinDict(sys.argv[2])

f = open(sys.argv[1])
snFile = open(sys.argv[3]+".sn","w")
asFile = open(sys.argv[3]+".as","w")

data = readFragOutFmt(sys.argv[1],sys.argv[2])

for key in data.keys():
	genes = data[key][1]
	if len(genes) > 0:
		if genes[0].strand == 0:
			snFile.write(">"+key + "\n")
			for gene in genes:
				snFile.write(gene.info[:-6]+"\n")
		else:
			asFile.write(">"+key + "\n")
			asFile.write(fastadict[key]+"\n")

f.close()
snFile.close()
asFile.close()
