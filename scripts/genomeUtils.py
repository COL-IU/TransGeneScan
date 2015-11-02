from RBTree import *
import re
import math
import sys

class Gene:
 def __init__(self,start,end,strand,index,seq):
  self.start = start
  self.end = end
  self.strand = strand # 0 for forward and 1 for reverse
  self.m1 = None
  self.m2 = None
  self.m3 = None
  self.info = None
  self.index = index
  self.seq = seq
  
 #def __cmp__(self, other):
 # if other == None:
 #  return False
 # return self.start == other.start and self.end == other.end and self.strand == other.strand
  
 def __hash__(self):
  return hash(str(self.start)+" "+str(self.end)+" "+str(self.strand))
  
 def write(self):
  print str(self.start) + "\t" + str(self.end),
  if self.strand == 1:
   print "\t-\t" + str(self.getFrame()-2),
  else:
   print "\t+\t" + str(self.getFrame()+1),
  print "" 
  
 def toString(self):
  retString = str(self.start) + "\t" + str(self.end)
  if self.strand == 0:
   retString = retString + "\t+\t" + str(self.getFrame()+1)
  else:
   retString = retString + "\t-\t" + str(self.getFrame()+1)
  return retString
	
 def getFrame(self):
  if self.seq == "":
   return -2
  seq = self.seq
  if seq[-3:] == "TAG" or seq[-3:] == "TGA" or seq[-3:] == "TAA":
   if self.strand == 0:
    return self.end%3
   else:
    return ((self.start-1)%3)+3
  elif seq[:3] == "ATG" or seq[:3] == "GTG" or seq[:3] == "TTG":
   if self.strand == 0:
    return (self.start-1)%3
   else:
    return (self.end%3)+3
  else:
   return -1
   
def compareFrames(frame1,frame2):
 if frame1 == -1 or frame2 == -1:
  return False
 return frame1 == frame2

def readFragOutFmt(filename,seqFile):
 if seqFile != None:
  seqDict = readFASTAinDict(seqFile)
 else:
  seqDict = {}
 fastalines = {}
 f = open(filename)
 lines = f.readlines()
 f.close()
 if lines[0][0] != '>':
  sys.exit("Invalid FASTA file: "+filename)
 geneNodes = []
 geneObjects = []
 header = ""
 for line in lines:
  if len(line) > 0:
   if line[0]=='>':
    line = line.rstrip()
    if header != "":
     fastalines[header] = [RBTree(geneNodes),geneObjects]
     geneNodes = []
     geneObjects = []
    header = line.strip()[1:]
   else:
    words = line.split(None)
    start = int(words[0])
    end = int(words[1])
    if header in seqDict:
     seq = seqDict[header][start-1:end]
    else:
     seq = ""
    if words[2] == '+':
     strand = 0
    else:
     strand = 1
     seq = reverse(complement(seq))
    gene = Gene(start,end,strand,len(geneNodes),seq)
    gene.info = line.strip()
    geneNodes.append([start,end,gene])
    geneObjects.append(gene)
 fastalines[header] = [RBTree(geneNodes),geneObjects]
 return fastalines

def complement(sequence):
#This function defines the nucleotide complement dictionary and utilizes it to return a complementary DNA sequence.
 dict = {"A":"T", "T":"A", "C":"G", "G":"C", "X":"X", "O":"O", "?":"?", "N":"N"}
 complement = ""
 for i in range(0,len(sequence)):
  complement = complement + dict[sequence[i]]
 return complement

def reverse(sequence):
#This function simply returns the reverse of a sequence. This function coupled with the above complement function
#will give a reverse complementary DNA sequence.
 return sequence[::-1]

def getMatch(gene,tree,overlapThreshold,ignoreFrame=False):
 queryStart = gene.start
 queryEnd = gene.end
 queryFrame = gene.getFrame()

 matchNode = tree.intervalSearch([queryStart,queryEnd])
 matchedGene = matchNode.obj
 if matchedGene == None:
  return -1
 matchedStart = matchedGene.start
 matchedEnd = matchedGene.end
 matchedFrame = matchedGene.getFrame()

 overlapPercent = 0.0
 if matchedStart != 0 and matchedEnd != 0:
  overlapNumber = 0.0 + min(matchedEnd,queryEnd)-max(matchedStart,queryStart)
  overlapPercent = (overlapNumber/(queryEnd-queryStart))*100
  if overlapPercent > overlapThreshold and (compareFrames(queryFrame,matchedFrame) or ignoreFrame):
   return [matchedGene,overlapPercent]
 return -1
 
def getExactMatch(gene,tree):
 queryStart = gene.start
 queryEnd = gene.end
 queryFrame = gene.getFrame()

 matchNode = tree.intervalSearch([queryStart,queryEnd])
 matchedGene = matchNode.obj
 if matchedGene == None:
  return -1
 matchedStart = matchedGene.start
 matchedEnd = matchedGene.end
 matchedFrame = matchedGene.getFrame()
 if matchedStart != 0 and matchedEnd != 0:
  if matchedStart == queryStart and matchedEnd == queryEnd and compareFrames(queryFrame,matchedFrame):
   return matchedGene
 return -1

def readFASTA(filename):
 fastalines = []
 f = open(filename)
 lines = f.readlines()
 f.close()
 if lines[0][0] != '>':
  sys.exit("Invalid FASTA file: "+filename)
 seq = ""
 header = ""
 for line in lines:
  if len(line) > 0:
   if line[0]=='>':
    if seq != "":
     fastalines.append([header,seq.upper()])
     seq = ""
    header = line.strip()[1:]
   else:
    line = re.sub(r'\s','',line)
    seq = seq + line
 if seq != "":
  fastalines.append([header,seq.upper()])
 return fastalines

def readFASTAinDict(filename):
 fastalines = {}
 f = open(filename)
 lines = f.readlines()
 f.close()
 if lines[0][0] != '>':
  sys.exit("Invalid FASTA file: "+filename)
 seq = ""
 header = ""
 for line in lines:
  if len(line) > 0:
   if line[0]=='>':
    if seq != "":
     words = header.split(None)
     fastalines[words[0][1:]] = seq.upper()
     seq = ""
    header = line
   else:
    line = re.sub(r'\s','',line)
    seq = seq + line
 if seq != "":
  words = header.split(None)
  fastalines[words[0][1:]] = seq.upper()
 return fastalines
