
from os.path import join
from os import listdir
import os
import re
import yaml
import boto3

configfile: "config/config.yaml"

RESOURCESYAML=config['resources']

with open(RESOURCESYAML) as json_file:
    CLUSTER = yaml.safe_load(json_file)
getthreads=lambda rname:int(CLUSTER[rname]["threads"]) if rname in CLUSTER and "threads" in CLUSTER[rname] else int(CLUSTER["__default__"]["threads"])
getmemg=lambda rname:CLUSTER[rname]["mem"] if rname in CLUSTER and "mem" in CLUSTER[rname] else CLUSTER["__default__"]["mem"]
getmemG=lambda rname:getmemg(rname).replace("g","G")
getmem_mb=lambda rname:int(getmemg(rname).replace("g",""))*1000


from snakemake.remote.S3 import RemoteProvider as S3RemoteProvider
S3 = S3RemoteProvider()

s3 = boto3.resource('s3')
my_bucket = s3.Bucket('ccr-genomics-testdata')

for object_summary in my_bucket.objects.filter(Prefix="chr"):
    print(object_summary.key)


SAMPLES = ["Test1_R_T","Test2_R_T"]

#workpath="./"


rule all:
	input:
#		expand(join("{sample}","fastqc","{sample}_R1_fastqc.zip"),sample=SAMPLES),
#		expand(join("{sample}","multiqc","multiqc_report.html"),sample=SAMPLES),
#		expand(join("{sample}","star_out","{sample}.star.bam"),sample=SAMPLES),
#		expand(join("{sample}","star_out","{sample}.star_transcriptome.bam"),sample=SAMPLES),
		expand(join("{sample}","rsem_out","{sample}.genes.results"),sample=SAMPLES)


rule fastqc:
	input:
		R1 = S3.remote("ccr-genomics-testdata/testdata/{sample}_R1.fastq.gz", keep_local=True),
		R2 = S3.remote("ccr-genomics-testdata/testdata/{sample}_R2.fastq.gz", keep_local=True),
	params:
		out = join("{sample}","fastqc")
	output:
		join("{sample}","fastqc","{sample}_R1_fastqc.zip"),
		join("{sample}","fastqc","{sample}_R2_fastqc.zip")

	conda: "env.yaml"

	threads: getthreads("fastqc")

	container: "docker://nciccbr/ccrgb_qctools:latest"

	shell: """

	fastqc {input} -o  {params.out}

	"""

rule multiqc:
	input:
		join("{sample}","fastqc","{sample}_R1_fastqc.zip"),
		join("{sample}","fastqc","{sample}_R2_fastqc.zip")
	params:
		out = join("{sample}","multiqc")
	output:
		join("{sample}","multiqc","multiqc_report.html")

	conda: "env.yaml"

	threads: getthreads("fastqc")
	
	container: "docker://nciccbr/ccrgb_qctools:latest"

	shell: """

	mkdir -p {params.out}
	multiqc {input} -o {params.out}

	"""

rule star:
	input:
		R1 = S3.remote("ccr-genomics-testdata/testdata/{sample}_R1.fastq.gz", keep_local=True),
		R2 = S3.remote("ccr-genomics-testdata/testdata/{sample}_R2.fastq.gz", keep_local=True),
		STARgenome = S3.remote("ccr-genomics-testdata/References/index-STAR_2.7.9a.tar", keep_local=True),
		gtf = S3.remote("ccr-genomics-testdata/References/gencode.v37lift37.annotation_ERCC92.gtf", keep_local=True),
	params:
		out = join("{sample}","star_out"),
#		STARgenome = S3.remote("ccr-genomics-testdata/References/index-STAR_2.7.9a.tar", keep_local=True),
#		gtf = S3.remote("ccr-genomics-testdata/References/gencode.v37lift37.annotation_ERCC92.gtf", keep_local=True),
	output:
		G_bam = join("{sample}","star_out","{sample}.star.bam"),
		T_bam = join("{sample}","star_out","{sample}.star_transcriptome.bam")
	
	container: "docker://nciccbr/ncigb_star_v2.7.10a:latest"

	threads: getthreads("star")
	resources:
		mem_mb=getmem_mb("star")

	conda: "env.yaml"

	shell: """
	
	indexdir="index-STAR_2.7.9a"

	if  [ ! -f ${{indexdir}}/SA ];
	then
	tar -xvf {input.STARgenome} 
	fi

	STAR --genomeDir $indexdir \
		--readFilesIn  {input.R1} {input.R2} \
		--sjdbGTFfile {input.gtf} \
		--readFilesCommand zcat \
		--runThreadN {threads} \
		--twopassMode Basic \
		--outSAMunmapped Within \
		--chimSegmentMin 12 \
		--chimJunctionOverhangMin 12 \
		--alignSJDBoverhangMin 10 \
		--alignMatesGapMax 100000 \
		--chimSegmentReadGapMax 3 \
		--outFilterMismatchNmax 2 \
		--outSAMtype BAM SortedByCoordinate \
		--quantMode TranscriptomeSAM \
		--outBAMsortingThreadN 6 \
		--limitBAMsortRAM 80000000000
	mv *Aligned.sortedByCoord.out.bam {output.G_bam}
	mv *Aligned.toTranscriptome.out.bam {output.T_bam}
	
	"""
	
rule rsem:
	input:
		bam = join("{sample}","star_out","{sample}.star_transcriptome.bam"),
		rsemindex = S3.remote("s3://ccr-genomics-testdata/References/rsem_1.3.2.tar", keep_local=True),
		
	params:
		out = join("{sample}","rsem_out"),
		sample = "{sample}"
	output:
		genes = join("{sample}","rsem_out","{sample}.genes.results"),
		isoform = join("{sample}","rsem_out","{sample}.isoforms.results"),
	
	threads: getthreads("star")
	resources:
		mem_mb=getmem_mb("star")

	container: "docker://nciccbr/ccbr_rsem_1.3.3:v1.0"

	conda: "env.yaml"

	shell: """

	rsemindexdir="rsem_1.3.2"
	
	if  [ ! -f ${{rsemindexdir}}/${{rsemindexdir}}.transcripts.fa ];
	then
	tar -xvf {input.rsemindex} 
	fi
	ls rsem_1.3.2
	#cd {params.out}
	rsem-calculate-expression --no-bam-output --paired-end -p {threads}  --estimate-rspd  --bam {input.bam} $rsemindexdir/$rsemindexdir {params.sample}
	cp {params.sample}.genes.results {output.genes}
	cp {params.sample}.isoforms.results {output.isoform}
	"""


