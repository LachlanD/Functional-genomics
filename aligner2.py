#!/usr/bin/env python

from Bio import SeqIO
import sys
import time

before = time.time()

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

# Builds a dictionary of kmers in the reference  
  def buildKmerDict(self):
    
    kmers = {}
    
    # Slide along the reference string extracting kmers
    for i in xrange(len(self.refseq)-self.kmerlength+1):
      kmer = self.refseq[i:i+self.kmerlength]
      if (kmer not in kmers):
	# New kmer
        kmers[kmer] = [i]
      else:
	# add to existing kmer
        kmers[kmer].append(i)
      
    return kmers

# Method for 
  def align(self, read):
    matches = []
    readlength = len(read)
    uread = read.seq.upper()
    ureverse = read.seq.reverse_complement().upper()

    # extract kmers from the forward read and match them to the reference dictionary
    for i in xrange(readlength-self.kmerlength+1):
      matches.append(self.dictionary.get(uread[i:i+self.kmerlength],[]))

    # find the runs in the forward strand by finding consecutive matches
    runs = {}
    for m in xrange(len(matches)):
      for a in matches[m]:
        i = 1
        while m+i < len(matches) and a+i in matches[m+i]:
          matches[m+i].remove(a+i)
          i += 1
        runs[a]  = i

    # extract kmers from the reverse read and match them to the reference dictionary 
    for i in xrange(len(read)-self.kmerlength+1):
      matches.append(self.dictionary.get(ureverse[i:i+self.kmerlength],[]))

    # find runs in the reverse strand by finding consecutive matches
    reverseruns = {}
    for m in xrange(len(matches)):
      for a in matches[m]:
        i = 1
        while m+i < len(matches) and a+i in matches[m+i]:
          matches[m+i].remove(a+i)
          i += 1
        reverseruns[a]  = i

    #seed the forward stand with any run equal to the best forward run length
    if len(runs) > 0:
      forwardseed =  max(list(runs.values()))
    
    #seed the backward strand with any run equal to the best reverse run length
    if len(reverseruns) > 0:
      backwardseed = max(list(reverseruns.values()))
    
    best = 3
    bestpos = 0
    beststrand = '*'
    # find all the runs with the maximum forward run length and check the hamming distance
    for pos, run in runs.iteritems():
      if run >=forwardseed:
	# Overhang is the readleangth - (the run length + the kmer length)        
        for i in xrange(max(0,pos-(readlength-self.kmerlength)),min(pos+1, len(self.refseq)-readlength)):
          hamdist = 0
          for j in xrange(len(read)):
            if uread[j] != self.refseq[i+j]:
              hamdist += 1
          # update the best hamming distance if a new best is found
          if hamdist < best:
            best = hamdist
            bestpos = i
            beststrand = '+'
          # update if an equally good alignment is found at an earlier position
          if hamdist == best:
            bestpos = min(i, bestpos)
    # find all the runs with the maximum backwards run length and check the hamming distance
    for pos, run in reverseruns.iteritems():
      if run >=backwardseed:
	# Overhang is the readleangth - (the run length + the kmer length)                       
        for i in xrange(max(0,pos-(readlength-self.kmerlength)),min(pos+1, len(self.refseq) - readlength)):
          hamdist = 0
          for j in range(len(read)):
            if ureverse[j] != self.refseq[i+j]:
              hamdist += 1
              
          # update the best hamming distance if a new best is found
          if hamdist < best:
            best = hamdist
            bestpos = i
            beststrand = '-'
	  # update if an equally good position is found at an earlier position
          if hamdist == best and beststrand == '-':
            bestpos = min(i, bestpos)
    
    # If the best is less than 3 hamming distance return that      
    if best < 3:  
      return Alignment(read.id, self.refname, bestpos+1, beststrand, best)
          
    return Alignment(read.id,  "*", 0, "*", 0) 


#Check the command line arguments
if len(sys.argv) < 3:
	print "Usage: <reference file (fasta)> <read file (fasta)> "
	sys.exit(0)


#Read the reference sequence and initiate the aligner
try:
  for s in SeqIO.parse(sys.argv[1], "fasta"):
    # Lets me add a kmer length to the command line for something other than 13 
    if(len(sys.argv)>3):
      aligner = Aligner2(s, int(sys.argv[3]))
    else:
      aligner = Aligner2(s)
    break #Stop after the fist sequence in the reference
except IOError as e:
  print "Could not read reference sequence file (see below)!"
  print e
  sys.exit(1)


#Open the read data and get to work
try:
  before = time.time()
   
  for read in SeqIO.parse(sys.argv[2], "fastq"):
    alignment = aligner.align(read) 
    print(str(alignment))
except IOError as e:
  print "Could not read fastq input file (see below)!"
  print e
  sys.exit(1)

after = time.time()

print "\nruntime " + str(after-before)

