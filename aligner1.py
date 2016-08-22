#!/usr/bin/env python

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
    strand = '+'
    for r in [read, read.reverse_complement()]:
      for i in range(refLength-readLength):
        ham = 0
        for j in range(readLength):
          if r[j] != self.refseq[i+j]:
            ham += 1
            if ham >= best:
              break
          
        if ham < best:
          best = ham
          bestPos = i+1 
          bestStrand = strand
        if ham == 0:
          break
      
      strand = '-'  

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
    #print(read)
    alignment = aligner.align(read) 
    print(str(alignment))
except IOError as e:
  print "Could not read fastq input file (see below)!" 
  print e 
  sys.exit(1)

