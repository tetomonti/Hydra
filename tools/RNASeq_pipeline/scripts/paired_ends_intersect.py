#!/usr/bin/python
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
