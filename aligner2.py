#!/usr/bin/env python3

from Bio import SeqIO
import sys

#Support Class to create alignments and check relationships between alignments
class Alignment:
  def __init__(self, readId, chr, pos, strand, hammingDistance):
    self.readId = readId
    self.chr = chr
    self.pos = pos
    self.strand = strand
    self.hammingDistance = hammingDistance
  def __repr__(self):
			return "\t".join([self.readId, self.chr, str(self.pos), self.strand, str(self.hammingDistance)])

#Support class to execute the handiwork of read alignment to a reference sequence
class Aligner2:
  def __init__(self, ref, klen=13):
    self.refname = ref.id
    self.refseq = ref.seq
    self.kmerlength = klen
    self.dictionary = self.buildKmerDict() 
  def buildKmerDict(self):
    
    kmers = {}
    
    for i in range(len(self.refseq)-self.kmerlength+1):
      kmer = self.refseq[i:i+self.kmerlength]
      if (kmer not in kmers):
        kmers[kmer] = [i]
      else:
        kmers[kmer].append(i)
      
    #print(kmers)

    return kmers

  def align(self, read):
    #print(read.seq)
    #print(self.dictionary)
    matches = []
    for i in range(len(read)-self.kmerlength+1):
      matches.append(self.dictionary.get(read.seq[i:i+self.kmerlength],[]))

    #print(matches)

    runs = {}
    for m in range(len(matches)):
      for a in matches[m]:
        i = 1
        while m+i < len(matches) and a+i in matches[m+i]:
          matches[m+i].remove(a+i)
          i += 1

        runs[a]  = i

    for i in range(len(read)-self.kmerlength+1):
      matches.append(self.dictionary.get(read.seq.reverse_complement()[i:i+self.kmerlength],[]))

    reverseruns = {}
    for m in range(len(matches)):
      for a in matches[m]:
        i = 1
        while m+i < len(matches) and a+i in matches[m+i]:
          matches[m+i].remove(a+i)
          i += 1

        reverseruns[a]  = i

    print(runs)
    print(reverseruns)  

    #the length of the second best run
    if len(runs) > 0:
      s =  sorted(list(runs.values()) + list(reverseruns.values()))[max(-1*len(runs),-2)]
    
    best = 3
    bestpos = 0
    beststrand = '*'
    for pos, run in runs.items():
      if run >=s:        
        for i in range(pos-1,pos+1):
          hamdist = 0
          for j in range(len(read)):
            if read.seq[j] != self.refseq[i+j]:
              hamdist += 1

          if hamdist < best:
            best = hamdist
            bestpos = i
            beststrand = '+'
          if hamdist == best:
            bestpos = min(i, bestpos)

    for pos, run in reverseruns.items():
      print(str(pos) + ' ' + str(run))
      if run >=s:
                
        for i in range(pos-1,pos+1):
          hamdist = 0
          for j in range(len(read)):
            if read.reverse_complement().seq[j] != self.refseq[i+j]:
              hamdist += 1
              
          if hamdist < best:
            best = hamdist
            bestpos = i
            beststrand = '-'
          if hamdist == best and beststrand == '-':
            bestpos = min(i, bestpos)
          
    if best < 3:  
      return Alignment(read.id, self.refname, bestpos+1, beststrand, best)
          
    
    
    
    return Alignment(read.id,  "*", 0, "*", 0) 


#Check the command line arguments
if len(sys.argv) < 3:
	print("Usage: <reference file (fasta)> <read file (fasta)> ")
	sys.exit(0)


#Read the reference sequence and initiate the aligner
try:
  for s in SeqIO.parse(sys.argv[1], "fasta"):
    #print(s)
    if(len(sys.argv)>3):
      aligner = Aligner1(s, int(sys.argv[3]))
    else:
      aligner = Aligner1(s)
    break #Stop after the fist sequence in the reference
except IOError as e:
  print("Could not read reference sequence file (see below)!")
  print(e)
  sys.exit(1)


#Open the read data and get to work
try: 
  for read in SeqIO.parse(sys.argv[2], "fastq"):
    #print(read)
    alignment = aligner.align(read) 
    print(str(alignment))
except IOError as e:
  print("Could not read fastq input file (see below)!")
  print(e)
  sys.exit(1)

