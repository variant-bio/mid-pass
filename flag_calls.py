####################################################################################
# flag_calls.py
# Simple script to merge post-imputation VCF with unfiltered pre-imputation VCF
# Required are two input files in the specified order: 
#  1) unfiltered.vcf.gz that contains variant calls pre-GQ-filtering
#  2) imputed.vcf.gz is the Beagle output VCF that contains imputed GT only
# The two vcf files are expected to contain the same variant sites in the same order and will error out otherwise
# Sample columns are expected to be in the same order, and only biallelic sites will be handled correctly.
#
#
#    Copyright (C) 2020 Anne-Katrin Emde
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
#
####################################################################################
import sys
import os
import gzip
from itertools import izip


# Check the number of command line arguments
if not len(sys.argv)==3:
    sys.stderr.write("\nError:\tincorrect number of command-line arguments")
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

            # use non-imputed calls to correct/flag imputed calls
            if (GTi == "0|0" or GTi == "1|1") and (GTs == "0/1"):
                GTo = GTi    # this used to be set to 0/1
                imputed = 2  # set IM flag to 2 - original lowGQ het call in disagreement with imputed call
            elif (GTi == "0|0" and GTs == "1/1") or (GTi == "1|1" and GTs == "0/0"):
                # could use a further distinction based on GQs whether to keep imputed, or set to het 
                GTo = GTi    # if phase is not needed, setting to 0/1 seems to be most accurate 
                imputed = 3  # set IM flag to 3 - original lowGQ hom call in disagreement with imputed opposite hom call
            else:   
                # Note: this script is not meant for multi-allelic sites but won't complain if they exist
                # if multi-allelic and allele > 1 is involved it'll always end up in this case and mark the call as imputed (IM=1) unless it agrees with the seq-based call (IM=0)
                # assuming up to a maximum of 9 alleles
                GTistr = ''.join(sorted(GTi[0] + GTi[2])) # turn 0|1 and 1|0 into 01  
                GTsstr = ''.join(sorted(GTs[0] + GTs[2]))
                if (GTsstr[0] == GTistr[0]) and (GTsstr[1] == GTistr[1]):
                    imputed = 0 
                else:
                    imputed = 1
                GTo = GTi
            
            # write out this sample's info with IM tag added
            indexStart = s.find(':')
            sys.stdout.write("\t" + GTo + s[int(indexStart):] + ":" + str(imputed))
    
            #if s.startswith("."):
            #    sys.stdout.write("\t" + GTo + s[int(indexStart):] + ":1")
            #else:
            #    sys.stdout.write("\t" + GTo + s[int(indexStart):] + ":" + str(imputed))
        sys.stdout.write("\n")



