#!/usr/bin/python
#  Copyright (c) 2014, Boston University. All rights reserved.
#  
#  Redistribution and use in source and binary forms, with or without
#  modification, are permitted provided that the following conditions are met: 
#  
#  1. Redistributions of source code must retain the above copyright notice, this
#     list of conditions and the following disclaimer. 
#  2. Redistributions in binary form must reproduce the above copyright notice,
#     this list of conditions and the following disclaimer in the documentation
#     and/or other materials provided with the distribution. 
#  
#  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
#  ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
#  WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
#  DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR
#  ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
#  (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
#  LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
#  ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
#  (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
#  SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
#  
#  The views and conclusions contained in the software and documentation are those
#  of the authors and should not be interpreted as representing official policies, 
#  either expressed or implied, of Boston University.
#  
#  Authors:
#    Daniel Gusenleitner [1,2], Vinay Kartha [1,2], Francesca Mulas [2], 
#    Yuxiang Tan [1,2], Liye Zhang [2], Stefano Monti [1,2]
#
#  [1] Bioinformatics Program, Boston University
#  [2] Center for Computational Biomedicine, Boston University  
#  

import sys
import re
import gzip

from Bio.SeqIO.QualityIO import FastqGeneralIterator #Biopython 1.51 or later

##########################################################
#
# Change the following settings to suit your needs
## Read paired-end fastq filenames from commandline arguments provided
#NOTE : THESE FILENAMES INCLUDE THE FULL PATH TO THE FILES
input_forward_filename=sys.argv[1]    
input_reverse_filename=sys.argv[2]
output_paired_forward_filename=sys.argv[3]
output_paired_reverse_filename=sys.argv[4]

#output_forward_filename_prefix=sys.argv[3]
#output_reverse_filename_prefix=sys.argv[4]
#output_orphan_filename_prefix=sys.argv[5]
#output_directory=sys.argv[6]



####################################       ORIGINAL SCRIPT BEGINS HERE ################################

#output_paired_forward_filename = "%s_matched.fastq.gz" %(output_directory+output_forward_filename_prefix)
#output_paired_reverse_filename = "%s_matched.fastq.gz" %(output_directory+output_reverse_filename_prefix)
#output_orphan_filename = "%s_unpaired_orphans.fastq" %(output_directory+output_orphan_filename_prefix)

# Illumina's fastq format specification in the header to differentiate between the forward and the reverse strands 


f_suffix = ""
r_suffix = ""

##########################################################

if f_suffix:
    f_suffix_crop = -len(f_suffix)
    def f_name(title):
        """Remove the suffix from a forward read name."""
        name = title.split()[0]
        assert name.endswith(f_suffix), name
        return name[:f_suffix_crop]
else:
    def f_name(title):
        return title.split()[0]

if r_suffix:
    r_suffix_crop = -len(r_suffix)
    def r_name(title):
        """Remove the suffix from a reverse read name."""
        name = title.split()[0]
        assert name.endswith(r_suffix), name
        return name[:r_suffix_crop]
else:
    def r_name(title):
        return title.split()[0]

print "Scanning reverse file to build list of names..."    
reverse_ids = set()
paired_ids = set()
for title, seq, qual in FastqGeneralIterator(gzip.open(input_reverse_filename,'rb')):
    reverse_ids.add(r_name(title))

print "Processing forward file..."
forward_handle = gzip.open(output_paired_forward_filename, "wb")
#orphan_handle = gzip.open(output_orphan_filename, "wb")
for title, seq, qual in FastqGeneralIterator(gzip.open(input_forward_filename,'rb')):
    name = f_name(title)
    if name in reverse_ids:
        #Paired
        paired_ids.add(name)
        reverse_ids.remove(name) #frees a little memory
        forward_handle.write("@%s\n%s\n+\n%s\n" % (title, seq, qual))
    else:
        #Orphan
#        orphan_handle.write("@%s\n%s\n+\n%s\n" % (title, seq, qual))
#Modified code
        continue
forward_handle.close()
del reverse_ids #frees memory, although we won't need more now

print "Processing reverse file..."
reverse_handle = gzip.open(output_paired_reverse_filename, "wb")
for title, seq, qual in FastqGeneralIterator(gzip.open(input_reverse_filename,'rb')):
    name = r_name(title)
    if name in paired_ids:
        #Paired
        reverse_handle.write("@%s\n%s\n+\n%s\n" % (title, seq, qual))
    else:
        #Orphan
#Modified code
        continue

#        orphan_handle.write("@%s\n%s\n+\n%s\n" % (title, seq, qual))
#orphan_handle.close()
reverse_handle.close()
print "Done"
