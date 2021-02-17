import pandas as pd

configfile: "config.yaml"

samples_fp = config['samples']


sample_table = pd.read_csv(samples_fp, sep='\t', header=0)
sample_table.set_index('Sample', inplace=True)

samples=sample_table.index
reads=['R1', 'R2']

def get_read(sample, read):
    return(sample_table.loc[sample, read])


rule all:
    input:
        "output/qc/multiqc.html"

rule fastqc:
    input:
        lambda wildcards: get_read(wildcards.sample,
                                   wildcards.read)
    output:
        html="output/qc/fastqc/{sample}.{read}.html",
        zip="output/qc/fastqc/{sample}.{read}_fastqc.zip" # the suffix _fastqc.zip is necessary for multiqc to find the file. If not using multiqc, you are free to choose an arbitrary filename
    params: ""
    log:
        "logs/fastqc/{sample}.{read}.log"
    threads: 1
    wrapper:
        "0.72.0/bio/fastqc"

rule cutadapt_pe:
    input:
        lambda wildcards: sample_table.loc[wildcards.sample,
                                      'R1'],
        lambda wildcards: sample_table.loc[wildcards.sample,
                                      'R2']
    output:
        fastq1="output/trimmed/{sample}.R1.fastq.gz",
        fastq2="output/trimmed/{sample}.R2.fastq.gz",
        qc="output/trimmed/{sample}.qc.txt"
    params:
        "-a {} {}".format(config["params"]["cutadapt"]['adapter'],
                          config["params"]["cutadapt"]['other'])
    log:
        "output/logs/cutadapt/{sample}.log"
    threads: 4
    wrapper:
        "0.17.4/bio/cutadapt/pe"


rule fastqc_post_trim:
    input:
        "output/trimmed/{sample}.{read}.fastq.gz"
    output:
        html="output/trimmed/fastqc/{sample}.{read}.trimmed.html",
        zip="output/trimmed/fastqc/{sample}.{read}.trimmed_fastqc.zip" # the suffix _fastqc.zip is necessary for multiqc to find the file. If not using multiqc, you are free to choose an arbitrary filename
    params: ""
    log:
        "logs/fastqc_post_trim/{sample}.{read}.log"
    threads: 1
    wrapper:
        "0.72.0/bio/fastqc"


rule multiqc:
    input:
        lambda wildcards: expand(rules.fastqc.output.zip,
                                 sample=samples,
                                 read=reads),
        lambda wildcards: expand(rules.fastqc_post_trim.output.zip,
                                 sample=samples,
                                 read=reads),
        lambda wildcards: expand(rules.cutadapt_pe.output.qc,
                                 sample=samples)

    output:
        "output/qc/multiqc.html"
    params:
        config['params']['multiqc']  # Optional: extra parameters for multiqc.
    log:
        "logs/multiqc.log"
    wrapper:
        "0.72.0/bio/multiqc"
