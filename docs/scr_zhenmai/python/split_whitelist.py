from Bio.SeqIO.QualityIO import FastqGeneralIterator
from gzip import open as gzopen

import argparse

ap = argparse.ArgumentParser()
ap.add_argument("-i", "--input", required=True, help="input fq file")
ap.add_argument("-o1", "--output1", required=True, help="output txt file")
ap.add_argument("-o2", "--output2", required=True, help="output txt file")

args = vars(ap.parse_args())

input_file = args["input"]
output_file1 = args["output1"]
output_file2 = args["output2"]

with gzopen(input_file, "rt") as in_handle, open(output_file1, "w") as out_handle1, open(output_file2, "w") as out_handle2:
	for title, seq, qual in FastqGeneralIterator(in_handle):
		x_cor = int(title.split(" ")[0].split(":")[-2])
		y_cor = int(title.split(" ")[0].split(":")[-1])
		barcode = seq[:-3]
		out_handle1.write("%s\t%d\t%d\n" % (barcode, x_cor, y_cor))
		out_handle2.write("%s\n" % barcode)