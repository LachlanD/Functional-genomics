#!/usr/bin/env python

from Bio import SeqIO
import sys
import time


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
    self.seqlen = len(self.refseq)
    self.dictionary = self.buildKmerDict()
    
# Builds a dictionary of kmers in the reference  
  def buildKmerDict(self):
    
    kmers = {}
    
    # Slide along the reference string extracting kmers
    for i in xrange(self.seqlen-self.kmerlength+1):
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
    fmatches = []
    readlength = len(read)
    uread = read.seq.upper()
    ureverse = read.seq.reverse_complement().upper()

    # extract kmers from the forward read and match them to the reference dictionary
    for i in xrange(readlength-self.kmerlength+1):
      fmatches.append(self.dictionary.get(uread[i:i+self.kmerlength],[]))

    # find the best forward run by finding the most consecutive matches
    bestfrun = 0
    bestfpos = []
    for m in xrange(len(fmatches)):
      for a in fmatches[m]:
        i = 1
        while m+i < len(fmatches) and a+i in fmatches[m+i]:
          fmatches[m+i].remove(a+i)
          i += 1
        if i > bestfrun:
          bestfrun = i
          bestfpos = [a]
        elif i == bestfrun:
          bestfpos.append(a)
    
    rmatches = []
    # extract kmers from the reverse read and match them to the reference dictionary 
    for i in xrange(len(read)-self.kmerlength+1):
      rmatches.append(self.dictionary.get(ureverse[i:i+self.kmerlength],[]))

    # find the best reverse runs by finding the most consecutive matches
    bestrrun = 0  
    bestrpos = []
    for m in xrange(len(rmatches)):
      for a in rmatches[m]:
        i = 1
        while m+i < len(rmatches) and a+i in rmatches[m+i]:
          rmatches[m+i].remove(a+i)
          i += 1
        if i > bestrrun:
          bestrrun = i 
          bestrpos = [a]
	elif i == bestrrun:
          bestrpos.append(a)		

    best = 3
    bestpos = 0
    beststrand = '*'

    for pos in bestfpos:
      if best == 0:
        break
      # Explore the run from 2 spots before its start until the end of the read overhangs atleast 2 positions or if the read length is longer than the run 2 spots past the run start
      start = max(0, pos -readlength + bestfrun)
      end = min(max(pos+2, pos + bestfrun - readlength +2), self.seqlen-readlength)
      for i in xrange(start, end):
        hamdist = 0
        for j in xrange(len(read)):
          if uread[j] != self.refseq[i+j]:
            hamdist += 1
        if hamdist < best:
          best = hamdist
          bestpos = i
          beststrand = '+'
          

    for pos in bestrpos:
      # Explore the run from 2 spots before its start until the end of the read overhangs atleast 2 positions or if the read length is longer than the run 2 spots past the run start
      start = max(0, pos -readlength + bestrrun)
      end = min(max(pos+2, pos + bestrrun - readlength + 2), self.seqlen-readlength)
      for i in xrange(start, end):
        hamdist = 0
        for j in xrange(len(read)):
          if ureverse[j] != self.refseq[i+j]:
            hamdist += 1
        if hamdist < best:
          best = hamdist
          bestpos = i
          beststrand = '-'
        # Use this if statement if earlier positions in the reverse strand are reported before the forward strand.  Spec and examples are ambiguous
        elif hamdist == best and i < bestpos:
	  bestpos = i
          beststrand = '-'
        
 
    
    # If the best is less than 3 hamming distance return that      
    if best < 3:  
      return Alignment(read.id, self.refname, bestpos+1, beststrand, best)
          
    return Alignment(read.id,  "*", 0, "*", 0) 

before = time.time()

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
  for read in SeqIO.parse(sys.argv[2], "fastq"):
    alignment = aligner.align(read) 
    print(str(alignment))
except IOError as e:
  print "Could not read fastq input file (see below)!"
  print e
  sys.exit(1)

after = time.time()

print "\nruntime " + str(after-before)
