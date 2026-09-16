from Bio.SeqIO.QualityIO import FastqGeneralIterator
from gzip import open as gzopen

import argparse

ap = argparse.ArgumentParser()
ap.add_argument("-i", "--input", required=True, help="input fq.gz of read1")
ap.add_argument("-o1", "--output_txt", required=True, help="read name output file")
ap.add_argument("-o2", "--output_fq", required=True, help="fastq output file")

args = vars(ap.parse_args())

input_r1 = args["input"]
output_file_txt = args["output_txt"]
output_file_fq = args["output_fq"]

with gzopen(input_r1, "rt") as in_handle, open(output_file_fq, "w") as out_handle, open(output_file_txt, "w") as fo:
	for title, seq, qual in FastqGeneralIterator(in_handle):
		ligation = seq[32:83]
		if ligation == 'GTTCGCAACATGTCTGGCGTCATAGAATTCCGCAGTCCAGTACGACTCACT':
			read2name = title.split(" ")[0] + '\n'
			out_handle.write("@%s\n%s\n+\n%s\n" % (title, seq, qual))
			fo.write(read2name)	