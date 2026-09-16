# cd /s1/tangzhidong/projectV1/SpatialOmicsZhenMai/
# SAMPLES = ["G00", "G06"]
SAMPLES = ["D00"]
REFERENCE = "doc/reference/STAR_index/"
GTF = "doc/reference/gencode.v47.basic.annotation.gtf"
WHITELIST = "doc/location/1/2D250811075_B1_barcode.txt"

rule all:
	input:
		expand("data/rawdata_cellranger/{sample}/{sample}_S1_L001_R1_001.fastq.gz", sample=SAMPLES),
		expand("data/rawdata_cellranger/{sample}/{sample}_S1_L001_R2_001.fastq.gz", sample=SAMPLES),
		expand("data/mapping/{sample}/{sample}_Aligned.sortedByCoord.out.bam", sample=SAMPLES)

rule filter_R1:
	input:
		R1_fq = "data/rawdata/1/{sample}_R1.fq.gz"
	output:
		R2_name_list = temp("data/rawdata/1/{sample}_R2name.list"),
		R1_fq_f = temp("data/rawdata/1/{sample}_R1_filter.fastq"),
		R1_fq_f_gz = temp("data/rawdata/1/{sample}_R1_filter.fastq.gz")
	params:
		scripts_path = "scr/python/filter_R1_byliga.py"
	shell:"""
		python {params.scripts_path} -i {input.R1_fq} -o1 {output.R2_name_list} -o2 {output.R1_fq_f}
		pigz -p 100 -k -f {output.R1_fq_f}
	"""

rule filter_R2:
	input:
		R2_fq = "data/rawdata/1/{sample}_R2.fq.gz",
		R2_name = rules.filter_R1.output.R2_name_list
	output:
		R2_fq_f = temp("data/rawdata/1/{sample}_R2_filter.fastq"),
		R2_fq_f_gz = temp("data/rawdata/1/{sample}_R2_filter.fastq.gz")
	shell:"""
		seqtk subseq {input.R2_fq} {input.R2_name} > {output.R2_fq_f}
		pigz -p 100 -k -f {output.R2_fq_f}
	"""

rule split_barcode:
	input:
		fq_R1 = "data/rawdata/1/{sample}_R1_filter.fastq.gz",
		fq_R2 = "data/rawdata/1/{sample}_R2_filter.fastq.gz"
	output:
		fq_R1_un = temp("data/rawdata_cellranger/{sample}/{sample}_S1_L001_R1_001.fastq"),
		fq_R1_gz = "data/rawdata_cellranger/{sample}/{sample}_S1_L001_R1_001.fastq.gz",
		fq_R2_gz = "data/rawdata_cellranger/{sample}/{sample}_S1_L001_R2_001.fastq.gz",
	params:
		split_script = "scr/python/split_bc.py"
	shell:"""
		python {params.split_script} -i {input.fq_R1} -o {output.fq_R1_un}
		pigz -p 100 -k -f {output.fq_R1_un}
		cp {input.fq_R2} {output.fq_R2_gz}
	"""

rule STARsolo:
	input:
		fq_R1 = rules.split_barcode.output.fq_R1_gz,
		fq_R2 = rules.split_barcode.output.fq_R2_gz
	output:
		bam = "data/mapping/{sample}/{sample}_Aligned.sortedByCoord.out.bam"
	params:
		Prefix = "data/mapping/{sample}/{sample}_"
	shell:"""
		STAR --runThreadN 16 \
		--genomeDir {REFERENCE} \
		--sjdbGTFfile {GTF} \
		--readFilesIn {input.fq_R2} {input.fq_R1} \
		--outFileNamePrefix {params.Prefix} \
		--readFilesCommand gunzip -c \
		--outSAMtype BAM SortedByCoordinate \
		--soloType CB_UMI_Complex \
		--soloCBmatchWLtype Exact \
		--soloCBwhitelist {WHITELIST} \
		--soloCBposition 0_0_0_31 \
		--soloUMIposition 0_32_0_41 \
		--soloFeatures Gene GeneFull
	"""