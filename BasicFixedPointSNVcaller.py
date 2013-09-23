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


# Setup
#bamFile="<input bam file>"
#reference="<referece>"
#vcf = "<vcf reference>"
#PhredCutoff = only bases above this will be considered for allele-freq-calculation (AFC)
#MapqCutoff = only reads with a mapq over this will be considered when calculating raw depth and AFC

#using arguments (in a lazy way) for this
bamFile = sys.argv[1]
reference = sys.argv[2]
vcf = sys.argv[3]
PhredCutoff = sys.argv[4]
#It is important to convert to int
PhredCutoff = int(PhredCutoff)
MapqCutoff = sys.argv[5]
MapqCutoff = int(MapqCutoff)
#special for extracting name from filename, in pilot
#NameOfSample = regex.findall(bamFile)[0]
#NameOfSample = "RNAseqSample"
NameOfSample = "Sample"

#print str(PhredCutoff) + " " + str(MapqCutoff)

MaxDepthCutoff = 1000000

Bases = ('A','C','G','T')

#using a Global dictionary, just to easily print the column-names to output...
baseCountGlobal = {'A': 0,'C': 0,'G': 0,'T': 0,'.': 0, 'N': 0 , 'RawDepth' : 0}
# Get base pileup from positon
def getBase(chrom,pos,bamFile):
     baseCount = {'A': 0,'C': 0,'G': 0,'T': 0,'CoverDel': 0, 'N': 0, 'RawDepth' : 0, 'FlankIns' : 0 ,'FlankDel' : 0}
     startOfRead = defaultdict(int)
     #for pileupColumn in bamFile.pileup(chrom, pos - 2, pos + 2):
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
                         #print str(pileupRead.alignment.seq[pileupRead.qpos])
                    #print str(phreadQual) + " " + str(pileupRead.alignment.mapq)
                    #print pileupRead.alignment.seq[pileupRead.qpos]
                    #if pileupRead.is_del:
                         #print "deletion with phred" + "\t" + str(phreadQual)
                         
                    
     return baseCount

# Input files
samFile = pysam.Samfile( bamFile, "rb" )
fastaFile = pysam.Fastafile(reference)
vcfFile = open(vcf,'r')

# Actual parsing and writing
#print "\tpos" + "\t" + "ref" + "\t" "alt" + "\t" + "hapmap" + "\t" + "A" + "\t" + "C" + "\t" + "G" + "\t" + "T" + "\t" + "." + "\t" + "N"
#changing so columns follows column-names
#print "\tpos" + "\t" + "ref" + "\t" "alt" + "\t" + "hapmap" + "\t" + "\t".join(str(x) for x in baseCountGlobal.keys())

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
          #ref =  line.split('\t')[3]
          #alt =  line.split('\t')[4]
          #alt = re.sub("\,.*","",alt)
          if "INDEL" in line:
               continue

          #genotypeInRef =  line.split('\t')[10].strip("\n")
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


          #Print-out for test               
          #print str(chrome) + "\t" + str(pos)  + "\t" +  str(baseCounts[refFromFasta]) + "\t" + refFromFasta + "\t" + str(altCount)  + "\t" +  altBase + "\t" + str(baseCounts['CoverDel'])  + "\t" + "Cover Deletion"  + "\t" + str(baseCounts['RawDepth']) + "\t" +  "RawDepth" + "\t" +  str(baseCounts['FlankIns']) + "\t" + "Flanking insertion" + "\t" +  str(baseCounts['FlankDel']) + "\t" + "Flanking deletion"
          
          allele_freq = 0
          if baseCounts[refFromFasta] + altCount > 0:
               allele_freq = altCount/(baseCounts[refFromFasta] + altCount)
          print str(chrome) + "\t" + str(pos) + "\t" + "."  + "\t" +  refFromFasta + "\t" + altBase + "\t" + "." + "\t" + "PASS" + "\t" + "." + "\t" + "DP4;Flankinsert,FlankDel,CoverDel"  + "\t" + str(baseCounts['A']) + "," + str(baseCounts['C']) + ","+ str(baseCounts['G']) + ","+ str(baseCounts['T'])  + ";" + str(baseCounts['FlankIns']) + "," + str(baseCounts['FlankDel']) + "," + str(baseCounts['CoverDel'])  + "\t" + str(baseCounts['RawDepth'])  + "\t" +  str(baseCounts[refFromFasta]) + "\t" +  str(altCount)  + "\t" + str(allele_freq)

#allele freq does not work if both are zero
#str(float(altCount)/float(altCount + baseCounts[refFromFasta]))

#"\trawdepth\tref_reads\talt_reads\tallele_freq\n"

samFile.close()
fastaFile.close()
vcfFile.close()

