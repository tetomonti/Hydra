#!/usr/bin/python
import sys
import gzip

from Bio.SeqIO.QualityIO import FastqGeneralIterator #Biopython 1.51 or later

##########################################################
#
# Change the following settings to suit your needs
## Read paired-end fastq filenames from commandline arguments provided
#NOTE : THESE FILENAMES INCLUDE THE FULL PATH TO THE FILES
INPUT_FORWARD_FILENAME = sys.argv[1]
INPUT_REVERSE_FILENAME = sys.argv[2]
OUTPUT_PAIRED_FORWARD_FILENAME = sys.argv[3]
OUTPUT_PAIRED_REVERSE_FILENAME = sys.argv[4]

def f_name(title):
    return title.split()[0]

def r_name(title):
    return title.split()[0]

print "Scanning reverse file to build list of names..."
reverse_ids = set()
paired_ids = set()
for title, seq, qual in FastqGeneralIterator(gzip.open(INPUT_REVERSE_FILENAME, 'rb')):
    reverse_ids.add(r_name(title))

print "Processing forward file..."
forward_handle = gzip.open(OUTPUT_PAIRED_FORWARD_FILENAME, "wb")
#orphan_handle = gzip.open(output_orphan_filename, "wb")
for title, seq, qual in FastqGeneralIterator(gzip.open(INPUT_FORWARD_FILENAME, 'rb')):
    nam = f_name(title)
    if nam in reverse_ids:
        #Paired
        paired_ids.add(nam)
        reverse_ids.remove(nam) #frees a little memory
        forward_handle.write("@%s\n%s\n+\n%s\n" % (title, seq, qual))
    else:
        #Orphan

        continue
forward_handle.close()
del reverse_ids #frees memory, although we won't need more now

print "Processing reverse file..."
reverse_handle = gzip.open(OUTPUT_PAIRED_REVERSE_FILENAME, "wb")
for title, seq, qual in FastqGeneralIterator(gzip.open(INPUT_REVERSE_FILENAME, 'rb')):
    nam = r_name(title)
    if nam in paired_ids:
        #Paired
        reverse_handle.write("@%s\n%s\n+\n%s\n" % (title, seq, qual))
    else:
        #Orphan
        continue

reverse_handle.close()
print "Done"
