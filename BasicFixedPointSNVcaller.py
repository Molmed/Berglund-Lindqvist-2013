# --------------------------------------------------------------------------
# Python script to get mutant allele fractions on prespecified positions
# --------------------------------------------------------------------------
#
# This script has been tested with python version 2.6.6
#
# Author: Carl Maarten Lindqvist  (but build upon script cowritten with Johan Dahlberg)
#
# Run example: 
# python BasicFixedPointSNVcaller.py bamfile_to_call.bam reference_samtools_indexed.fa ListWithPostionsToCheck.vcf phredScoreCutoff mapqCutoff 
# python BasicFixedPointSNVcaller.py bamfile_to_call.bam reference_samtools_indexed.fa ListWithPostionsToCheck.vcf 13 0 




from __future__ import division
import pysam
import sys
import re
from collections import defaultdict

#using arguments
bamFile = sys.argv[1]
reference = sys.argv[2]
vcf = sys.argv[3]
PhredCutoff = sys.argv[4]
#It is important to convert to int
PhredCutoff = int(PhredCutoff)
MapqCutoff = sys.argv[5]
MapqCutoff = int(MapqCutoff)

NameOfSample = "Sample"

MaxDepthCutoff = 1000000

Bases = ('A','C','G','T')

#using a Global dictionary, just to easily print the column-names to output...
baseCountGlobal = {'A': 0,'C': 0,'G': 0,'T': 0,'.': 0, 'N': 0 , 'RawDepth' : 0}
# Get base pileup from positon
def getBase(chrom,pos,bamFile):
     baseCount = {'A': 0,'C': 0,'G': 0,'T': 0,'CoverDel': 0, 'N': 0, 'RawDepth' : 0, 'FlankIns' : 0 ,'FlankDel' : 0}
     startOfRead = defaultdict(int)

     for pileupColumn in bamFile.pileup(chrom, pos - 2, pos + 2,max_depth=MaxDepthCutoff):
          if(pileupColumn.pos == pos - 1):
               for pileupRead in pileupColumn.pileups:
                    startOfRead[(pileupRead.qpos)] += 1
                    phreadQual = ord(pileupRead.alignment.qual[pileupRead.qpos]) - 33
                    #baseCount[pileupRead.alignment.seq[pileupRead.qpos]] += 1
                    #Rawdepth is reads with MapQ over cutoff (0)
                    if pileupRead.alignment.mapq > MapqCutoff:
                         baseCount['RawDepth'] += 1
                    if pileupRead.indel > 0:
                         baseCount['FlankIns'] += 1
                    if pileupRead.indel < 0:
                         baseCount['FlankDel'] += 1
                    #print str(pileupRead.is_del) + " " + str(pileupRead.alignment.seq[pileupRead.qpos])
                    if pileupRead.is_del:
                         baseCount['CoverDel'] += 1
                    #add check that base is ACGT ...
                    elif phreadQual > PhredCutoff and pileupRead.alignment.mapq > MapqCutoff and pileupRead.alignment.seq[pileupRead.qpos] in Bases:
                         baseCount[pileupRead.alignment.seq[pileupRead.qpos]] += 1
                    
                    
                         
                    
     return baseCount

# Input files
samFile = pysam.Samfile( bamFile, "rb" )
fastaFile = pysam.Fastafile(reference)
vcfFile = open(vcf,'r')



print "##fileformat=VCFv4.1"
print "##this is called with BasicFixedPointSNVcaller.py"
print "##headers are not commented, otherwise fulfilling the vcf-standard"
print "CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + NameOfSample +  "\trawdepth\tref_reads\talt_reads\tallele_freq"

for line in vcfFile:                                                                                                                                                                                  
     if not (line.startswith('#') or line.startswith('CHROM')):
          chrome = line.split('\t')[0]
          pos  = line.split('\t')[1]
          pos=int(pos)
          #using ref from fasta and calculates alt from baseCounts
          if "INDEL" in line:
               continue
          refFromFasta = fastaFile.fetch(chrome,pos-1,pos)
          baseCounts = getBase(chrome,pos,samFile)
          if baseCounts == None:
               print "Warning!!!"
               continue
          #Calculate which (if any) of the non-reference-bases has the highest frequency
          altBase = '.'
          altCount = 0
          for i in Bases:
               if i != refFromFasta:
               #print str(baseCounts[i]) + "\t" + i
                    if baseCounts[i] > altCount:
                         altCount = baseCounts[i]
                         altBase = i


          allele_freq = 0
          if baseCounts[refFromFasta] + altCount > 0:
               allele_freq = altCount/(baseCounts[refFromFasta] + altCount)
          print str(chrome) + "\t" + str(pos) + "\t" + "."  + "\t" +  refFromFasta + "\t" + altBase + "\t" + "." + "\t" + "PASS" + "\t" + "." + "\t" + "DP4;Flankinsert,FlankDel,CoverDel"  + "\t" + str(baseCounts['A']) + "," + str(baseCounts['C']) + ","+ str(baseCounts['G']) + ","+ str(baseCounts['T'])  + ";" + str(baseCounts['FlankIns']) + "," + str(baseCounts['FlankDel']) + "," + str(baseCounts['CoverDel'])  + "\t" + str(baseCounts['RawDepth'])  + "\t" +  str(baseCounts[refFromFasta]) + "\t" +  str(altCount)  + "\t" + str(allele_freq)

samFile.close()
fastaFile.close()
vcfFile.close()

