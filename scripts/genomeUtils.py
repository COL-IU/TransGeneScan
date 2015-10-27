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

def testPrint():
 print "hi"

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

def readRegulonFile(regFile,seqFile):
 singleSeqDs = readFASTA(seqFile)
 genome = singleSeqDs[0][1]
 geneNodes = []
 geneObjects = []
 f = open(regFile)
 lines = f.readlines()
 f.close()
 for i in range(1,len(lines)):
  line = lines[i]
  words = line.split(",")
  if len(words) > 2:
   start = int(words[2])
   end = int(words[3])
   if words[6] == 'forward':
    strand = 0
    seq = genome[start-1:end]
   else:
    strand = 1
    seq = reverse(complement(genome[start-1:end]))
   gene = Gene(start,end,strand,len(geneNodes),seq)
   gene.info = line.strip()
   geneNodes.append([start,end,gene])
   geneObjects.append(gene)
 return [RBTree(geneNodes),geneObjects] 
 
def readPttFile(pttFile,seqFile):
 singleSeqDs = readFASTA(seqFile)
 genome = singleSeqDs[0][1]
 geneNodes = []
 geneObjects = []
 f = open(pttFile)
 lines = f.readlines()
 f.close()
 for i in range(3,len(lines)):
  line = lines[i]
  words = line.split("\t")
  if len(words) > 2:
   startEnd = words[0].split("..")
   start = int(startEnd[0])
   end = int(startEnd[1])
   if words[1] == '+':
    strand = 0
    seq = genome[start-1:end]
   else:
    strand = 1
    seq = reverse(complement(genome[start-1:end]))
   gene = Gene(start,end,strand,len(geneNodes),seq)
   gene.info = line.strip()
   geneNodes.append([start,end,gene])
   geneObjects.append(gene)
 return [RBTree(geneNodes),geneObjects]

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
'''
def getLargerMatch(gene,tree,overlapThreshold):
 overlap = True
 start = gene.start
 end = gene.end
 frame = gene.frame

 match = tree.intervalSearch([start,end])

 matchedGene = match.key
 low = matchedGene["low"]
 high = matchedGene["high"]
 matchFrame = match.frame

 if frame == matchFrame:
  overlapPercent = 0.0
  if low != 0 and high != 0:
   overlapNumber = 0.0 + min(high,end)-max(low,start)
   overlapPercent = (overlapNumber/min((high-low),(end-start)))*100
   if overlapPercent > overlapThreshold:
    return match.info[0]
 return -1
'''
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

def readClustalW(filename):
 data = {}
 f = open(filename)
 lines = f.readlines()
 f.close()
 header = ""
 for i in range(1,len(lines)):
  line = lines[i]
  words = line.split(None)
  if len(words) != 0:
   if words[0] not in data.keys():
    if header == "":
     header = words[0]
    data[words[0]] = words[1].strip()
   else:
    data[words[0]] = data[words[0]] + words[1].strip()
 return [header,data] 

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
     fastalines.append([header,seq])
     seq = ""
    header = line.strip()[1:]
   else:
    line = re.sub(r'\s','',line)
    seq = seq + line
 if seq != "":
  fastalines.append([header,seq])
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
     fastalines[words[0][1:]] = seq
     seq = ""
    header = line
   else:
    line = re.sub(r'\s','',line)
    seq = seq + line
 if seq != "":
  words = header.split(None)
  fastalines[words[0][1:]] = seq
 return fastalines
'''
def getCodonMutCount(gene,xo):
 start = gene.start
 end = gene.end
 frame = gene.frame
 n1 = 0.0
 n2 = 0.0
 n3 = 0.0
 if start == 0 and end == 0:
  return [n1,n2,n3,0]
 if frame < 3:
  extract = xo[start-1:end]
 else:
  extract = reverse(complement(xo[start-1:end]))
 mutscore = 0.0
 for i in range(0,len(extract),3):
  if extract[i] == 'O':
   n1 = n1 + 1
   mutscore = mutscore + math.log (0.2) -math.log(0.33)
  if i+1 < len(extract):
   if extract[i+1] == 'O':
    n2 = n2 + 1
    mutscore = mutscore + math.log (0.2) -math.log(0.33)
  if i+2 < len(extract):
   if extract[i+2] == 'O':
    n3 = n3 + 1
    mutscore = mutscore + math.log (0.6) -math.log(0.33)
 #total = n1+n2+n3
 #print str(n1*100/total) + "\t" + str(n2*100/total) + "\t" + str(n3*100/total)
 #print str(start) + "\t" + str(end) + "\t" + str(n1) + "\t" + str(n2) + "\t" + str(n3)
 n4 = extract.count('?')
 return [n1,n2,n3,n4,mutscore]
'''
