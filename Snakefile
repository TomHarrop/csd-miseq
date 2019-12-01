#!/usr/bin/env Python3

import multiprocessing
import pandas
from pathlib import Path

def resolve_raw_fastq(wildcards):
    return {'r1': sample_data.loc[wildcards.sample_name, 'r1_path'],
            'r2': sample_data.loc[wildcards.sample_name, 'r2_path']}

sample_table_file = 'data/sample_table.csv'
ref = 'data/GCF_003254395.2_Amel_HAv3.1_genomic.fna'

# for filtering and trimming
bbduk_ref = '/phix174_ill.ref.fa.gz'
bbduk_adaptors = '/adapters.fa'

# set up directories
outdir = Path('output')
logdir = Path(outdir, 'logs')
tmpdir = Path(outdir, 'tmp')
mergedir = Path(outdir, '010_merge')
mapdir = Path(outdir, '020_map')

# singularity
bbmap = 'shub://TomHarrop/singularity-containers:bbmap_38.00'
bwa = 'shub://TomHarrop/singularity-containers:bwa_0.7.17'
samtools = 'shub://TomHarrop/singularity-containers:samtools_1.9'


sample_data = pandas.read_csv(sample_table_file,
                              index_col='sample_name')

all_samples = sorted(set(sample_data.index))
all_samples = [x for x in all_samples if x != 'Undetermined']

rule target:
    input:
        Path(outdir, 'merged.bam'),
        Path(outdir, 'merged.bam.bai')

rule merge_bam:
    input:
        bam = expand(Path(mapdir, '{sample_name}.bam').as_posix(),
                     sample_name=all_samples),
    output:
        Path(outdir, 'merged.bam')
    log:
        Path(logdir, 'merge_bam.log')
    threads:
        multiprocessing.cpu_count()
    singularity:
        samtools
    shell:
        'samtools merge '
        '-l 9 '
        '-O BAM '
        '-@ {threads} '
        '{output} '
        '{input.bam} '
        '2> {log}'


rule sort:
    input:
        Path(tmpdir, '{sample_name}.sam')
    output:
        Path(mapdir, '{sample_name}.bam')
    threads:
        multiprocessing.cpu_count()
    log:
        Path(logdir, '{sample_name}_sort.log')
    singularity:
        samtools
    shell:
        'samtools sort '
        '-@ {threads} '
        '{input} '
        '>> {output} '
        '2> {log}'


rule bwa_map:
    input:
        fq = Path(mergedir, '{sample_name}.fastq'),
        index = expand(
            Path(outdir, '000_ref', 'ref.fasta.{suffix}').as_posix(),
            suffix=['amb', 'ann', 'bwt', 'pac', 'sa'])
    output:
        temp(Path(tmpdir, '{sample_name}.sam'))
    params:
        prefix = Path(outdir, '000_ref', 'ref.fasta'),
        rg = '\'@RG\\tID:{sample_name}\\tSM:{sample_name}\''
    threads:
        multiprocessing.cpu_count()
    log:
        Path(logdir, '{sample_name}_bwa-map.log')
    singularity:
        bwa
    shell:
        'bwa mem '
        '-t {threads} '
        '-p '
        '-R {params.rg} '
        '{params.prefix} '
        '{input.fq} '
        '> {output} '
        '2> {log}'


rule merge:
    input:
        Path(tmpdir, '{sample_name}_trim.fastq')
    output:
        fastq = Path(mergedir, '{sample_name}.fastq'),
        singles = Path(mergedir, '{sample_name}_singles.fastq'),
        ihist = Path(mergedir, '{sample_name}_ihist.txt')
    log:
        Path(logdir, '{sample_name}_merge.log')
    singularity:
        bbmap
    shell:
        'bbmerge.sh '
        'in={input} '
        'int=t '
        'verystrict=t '
        'out={output.fastq} '
        'outu={output.singles} '
        'ihist={output.ihist} '
        '2> {log}'


rule trim:
    input:
        Path(tmpdir, '{sample_name}_filter.fastq')
    output:
        pipe = temp(Path(tmpdir, '{sample_name}_trim.fastq')),
        stats = Path(mergedir, '{sample_name}_trim.txt')
    log:
        Path(logdir, '{sample_name}_trim.log')
    params:
        trim = bbduk_adaptors
    singularity:
        bbmap
    shell:
        'bbduk.sh '
        'in={input} '
        'int=t '
        'out=stdout.fastq '
        'ref={params.trim} '
        'ktrim=r k=23 mink=11 hdist=1 tpe tbo '
        'forcetrimmod=5 '
        'stats={output.stats} '
        '>> {output.pipe} '
        '2> {log}'

rule filter:
    input:
        Path(tmpdir, '{sample_name}_repair.fastq')
    output:
        pipe = temp(Path(tmpdir, '{sample_name}_filter.fastq')),
        stats = Path(mergedir, '{sample_name}_filter.txt')
    log:
        Path(logdir, '{sample_name}_filter.log')
    params:
        filter = bbduk_ref,
    singularity:
        bbmap
    shell:
        'bbduk.sh '
        'in={input} '
        'int=t '
        'out=stdout.fastq '
        'ref={params.filter} '
        'hdist=1 '
        'stats={output.stats} '
        '>> {output.pipe} '
        '2> {log}'

rule check_pairing:
    input:
        unpack(resolve_raw_fastq)
    output:
        pipe = temp(Path(tmpdir, '{sample_name}_repair.fastq')),
    log:
        Path(logdir, '{sample_name}_repair.txt')
    singularity:
        bbmap
    shell:
        'repair.sh '
        'in={input.r1} '
        'in2={input.r2} '
        'out=stdout.fastq '
        '>> {output.pipe} '
        '2> {log}'

# bwa index rule
rule bwa_index:
    input:
        ref
    output:
        expand(Path(outdir, '000_ref', 'ref.fasta.{suffix}').as_posix(),
               suffix=['amb', 'ann', 'bwt', 'pac', 'sa'])
    params:
        prefix = Path(outdir, '000_ref', 'ref.fasta')
    log:
        Path(logdir, 'bwa_index.log')
    singularity:
        bwa
    shell:
        'bwa index '
        '-p {params.prefix} '
        '{input} '
        '2> {log}'


# generic bamfile index
rule index_bamfile:
    input:
        Path('{folder}', '{file}.bam')
    output:
        Path('{folder}', '{file}.bam.bai')
    log:
        Path(logdir, '{folder}', '{file}_index-bamfile.log')
    threads:
        2
    shell:
        'samtools index -@ {threads} {input} 2> {log}'
