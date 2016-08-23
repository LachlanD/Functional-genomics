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
class Aligner1:
  def __init__(self, ref):
    self.refname = ref.id
    self.refseq = ref.seq
  def align(self, read):
    readLength = len(read)
    refLength = len(self.refseq)
    best = 3
    bestPos = 0
    bestStrand = '*'
    uread = read.upper()
    ureverse = read.reverse_complement().upper()
    # Note. This implementation will report alignments closer to zero on the reverse strand before later alignments on forward strand spec/examples ambiguous
    for i in xrange(refLength-readLength+1):
      # calculate hamming distance on forward strand
      ham = 0
      for j in xrange(readLength):
        if uread[j] != self.refseq[i+j]:
          ham += 1
          if ham > best:
            break
      if ham < best:
        best = ham
        bestPos = i+1 
        bestStrand = '+'

      # calculate hamming distance on reverse strand
      ham = 0
      for j in xrange(readLength):
        if ureverse[j] != self.refseq[i+j]:
          ham += 1
          if ham > best:
            break
      if ham < best:
        best = ham
        bestPos = i+1 
        bestStrand = '-'
      # if we have already found a perfect alignment no need to look further
      if best == 0:
        break          
      
     

    if best > 2:	
      alignment = Alignment(read.id, "*", 0, "*", 0)
    else:
      alignment = Alignment(read.id, self.refname, bestPos, bestStrand, best)
    return alignment


#Check the command line arguments
if len(sys.argv) < 3:
	print "Usage: <reference file (fasta)> <read file (fasta)> " 
	sys.exit(0)


#Read the reference sequence and initiate the aligner
try:
  for s in SeqIO.parse(sys.argv[1], "fasta"):
    #print(s)
    aligner = Aligner1(s)
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
