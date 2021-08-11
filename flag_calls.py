###############################################################################
### COPYRIGHT #################################################################
#
# Variant Bio, Inc.
#
# SOFTWARE COPYRIGHT NOTICE AGREEMENT
# This software and its documentation are copyright (2021) by Variant Bio.
# All rights are reserved. This software is supplied without any warranty
# or guaranteed support whatsoever. Variant Bio cannot be responsible for
# its use, misuse, or functionality. 
#
# Author: Anne-Katrin Emde
#
# version 1.0 flag_calls.py
#
###############################################################################
############################################################################### 
# 
# Script to merge post-imputation VCF with unfiltered pre-imputation VCF
# Required are two input files in the specified order: 
#  1) unfiltered.vcf.gz that contains variant calls pre-GQ-filtering
#  2) imputed.vcf.gz is the Beagle output VCF that contains imputed GT only
# The two vcf files are expected to contain the same variant sites in the same 
# order and will error out otherwise. Sample columns are expected to be in the
# same order.
#
###############################################################################
import sys
import os
import gzip
from itertools import izip


# Check the number of command line arguments
if not len(sys.argv)==3:
    sys.stderr.write("\nError:\tincorrect number of command-line arguments\n")
    sys.stderr.write("Syntax:\tflag_calls.py unfiltered.vcf.gz imputed.vcf.gz\n")
    sys.stderr.write("\tOutput vcf written to stdout\n")
    sys.exit()


# open both files
with gzip.open(sys.argv[1], "r") as file1, gzip.open(sys.argv[2], "r") as file2: 

    # take care of header first
    headerline_added = False
    for x in file1:
        if x.startswith("##"):
            if (not headerline_added) and x.startswith("##FORMAT"):
                sys.stdout.write("##FORMAT=<ID=IM,Number=1,Type=Integer,Description=\"Indicates whether GT was original call (0), imputed (1) or in disagreement with filtered heterozygous (2) or homozygous call (3).\">\n")
                headerline_added = True
        sys.stdout.write(x)
        # last VCF header, just output file1 line
        if x.startswith("#") and not x.startswith("##"):
            break

    # skip through header in file2 
    for y in file2:
        if y.startswith("#") and not y.startswith("##"):
            break

    # go through each line
    for x, y in izip(file1, file2):
        x = x.strip()
        y = y.strip()
       
        xcols = x.split("\t")
        ycols = y.split("\t")

        # if x contains sites that were filtered from y -> skip ahead in x
        while (xcols[1] < ycols[1]):
            x=next(file1)
            x = x.strip()
            xcols = x.split("\t")


        if not xcols[1] == ycols[1]:

            sys.stderr.write("Variants in VCF input files differ:\n")
            sys.stderr.write("File 1:" + "\n")
            sys.stderr.write(x + "\n")
            sys.stderr.write("File 2:" + "\n")
            sys.stderr.write(y + "\n\n")
            break

        # output line up to FORMAT column, add :IM to FORMAT
        format = xcols[8] + ":IM"
        varinfoX = xcols[:8]
        GQ_index = xcols[8].split(":").index("GQ")
        GT_index = xcols[8].split(":").index("GT")

        for v in varinfoX:
            sys.stdout.write(v + "\t")
        sys.stdout.write(format)

        # extract sample info - check number of samples are the same 
        samplesX = xcols[9:]
        samplesY = ycols[9:]
        if not len(samplesX) == len(samplesY):
            sys.stderr.write("File1 and file2 have different number of samples at position: {}:{}\n\n".format(xcols[0], xcols[1]))
            break
        
        # compare genotypes for unfiltered and imputed GT for each sample
        for s,b in izip(samplesX, samplesY):
            imputed = 0
            # s is sequence based info
            GTs = s.split(":")[GT_index]
            GQs = s.split(":")[GQ_index]
            if GTs == ".":
                GTs = "./."
            if GQs == ".":
                GQs = 0
            else:
                GQs = int(GQs)
            # i is imputed GT
            GTi = b.split(":")[GT_index]
            
            # strip phase from the GATK VCF
            GTs = GTs.replace("|","/")
            if GTs == "1/0": GTs = "0/1"
            # Note: assuming up to a maximum of 9 alleles
            GTistr = ''.join(sorted(GTi[0] + GTi[2])) # turn 0|1 and 1|0 into 01  
            GTsstr = ''.join(sorted(GTs[0] + GTs[2]))
 
            # use non-imputed calls to correct/flag imputed calls
            # wgs called a het and imputed and wgs call are not the same - potentially seq/mapping error/contamination
            if (GTsstr[0] != GTsstr[1]) and (GTistr[0]!=GTsstr[0] or GTistr[1]!=GTsstr[1]):
                GTo = GTi    # this used to be set to het
                imputed = 2  # set IM flag to 2 - original lowGQ het call in disagreement with imputed call
            # wgs call is hom and neither imputed allele agrees
            elif (GTsstr[0] != ".") and (GTsstr[0] == GTsstr[1]) and not (GTistr[0]==GTsstr[0] or GTistr[1]==GTsstr[0]):
                # could use a further distinction based on GQs whether to keep imputed, or set to het 
                GTo = GTi    # if phase is not needed, setting to het seems to be most accurate 
                imputed = 3  # set IM flag to 3 - original lowGQ hom call in disagreement with imputed call (usually opposite hom)
            else:   
                # imputed and wgs call are the same
                if (GTsstr[0] == GTistr[0]) and (GTsstr[1] == GTistr[1]):
                    imputed = 0 
                    # should check if GQ value is below or above threshold and split the 0 cases into 
                # all other cases (either wgs was uncalled or no conflict e.g. 0/0 -> 0/1)
                else:
                    imputed = 1
                # we do this anyway for all cases    
                GTo = GTi
 
           
            # write out this sample's info with IM tag added
            indexStart = s.find(':')
            sys.stdout.write("\t" + GTo + s[int(indexStart):] + ":" + str(imputed))
    
            # if uncalled GT then always use IM=1
            #if s.startswith("."):
            #    sys.stdout.write("\t" + GTo + s[int(indexStart):] + ":1")
            #else:
            #    sys.stdout.write("\t" + GTo + s[int(indexStart):] + ":" + str(imputed))
        sys.stdout.write("\n")



