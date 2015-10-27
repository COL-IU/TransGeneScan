#! /usr/bin/python
import sys
from genomeUtils import *
import re

if len(sys.argv) < 4:
 sys.exit("Usage: getNegTranscripts.py <consensusFile> <outPrefix> <genome>")
 
consensusFilePath = sys.argv[1]
fastaFile = open(sys.argv[2]+".fasta","w")
genomeFilePath = sys.argv[3]

fi = open(consensusFilePath)
lines = fi.readlines()
firstWords = lines[0].split(None)

prevPos = int(firstWords[1])
seq = firstWords[3]
start = int(firstWords[1])
prevNuc = firstWords[3]

transcripts = []

for i in range(1,len(lines)):
 line = lines[i]
 words = line.split(None)
 curPos = int(words[1])
 refBase = words[2]
 curNuc = words[3]
 if curNuc != '-':
  if ((curPos-1) == prevPos) and prevNuc != '-':
   seq = seq + curNuc
  else:
   if start != 0:
    if prevPos-start > 120:
     transcripts.append([start,prevPos,seq])
   start = curPos
   seq = curNuc
  prevPos = curPos
 else:
  if seq != "":
   if prevPos-start > 120:
    transcripts.append([start,prevPos,seq])
  start = 0
  seq = ""
 prevNuc = curNuc
   
for transcript in transcripts:
 seq = reverse(complement(transcript[2]))
 fastaFile.write(">rtranscript:"+str(transcript[0])+":"+str(transcript[1])+"\n"+seq+"\n")
 
fi.close()
fastaFile.close()
