from Bio.SeqIO.QualityIO import FastqGeneralIterator
from Bio.Seq import Seq
from gzip import open as gzopen

import argparse

ap = argparse.ArgumentParser()
ap.add_argument("-i", "--input", required=True, help="input fq file")
ap.add_argument("-o", "--output", required=True, help="output txt file")

args = vars(ap.parse_args())

input_file = args["input"]
output_file = args["output"]

bc_start=0
bc_end=32

# umi_start=83
# umi_end=93

with gzopen(input_file, "rt") as in_handle, open(output_file, "w") as out_handle:
	for title, seq, qual in FastqGeneralIterator(in_handle):
		barcode = Seq(seq[bc_start:bc_end])
		barcode = barcode.reverse_complement()
		out_handle.write("%s\n" % barcode)
