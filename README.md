# Decode-seq

## Description

The Decode-seq data analysis pipeline contains 3 python scripts: decode_barcode.py, decode_gene.py, decode_quant.py. We need to run these three scripts in sequential order on the input fastq files and the output is the count based quantification matrix, which is ready for downstream analysis by edgeR or DEseq2.

- ```decode_barcode.py```: take Read1 as the input, generate a barcode table (read name, USI, and UMI).
- ```decode_gene.py```: take Read2 as the input, generate a gene table (read name, transcript).
- ```decode_quant.py```: take two tables as the input, generate the count matrix (transcript X sample).

## Dependencies

- Python 3
- [STAR](https://github.com/alexdobin/STAR) 

## Usage

```
decode_barcode.py -i|--input <input> -u|--usi <usifile> [-q|--qscutoff <int>] [-b|--boundary GGG] [-c|--countonly]

decode_gene.py -f|--fastq <input> -d|--outdir <outdir> [-x|--STARIDX <STARIDX>] [-g|--gtf gtf] [-t|--thread threads] [-b|--bam bam]

decode_quant.py -g|--genetable <gene_table> -b|--barcodetable <barcode_table> -u|--usi <usifile>
```

## Arguments

decode_barcode.py
- -i,--input <input>    Read1 fastq file
- -u,--usi <usifile>    A plain text file caontaining 2 columns:  USI name and 6bp sequence.
- -q,--qscutoff <int>    Minimum sequencing quality score of 1-6bp of Read1 (the USI position). Default 20
- -b,--boundary <NNN>    3-6 base between UMI sequence and cDNA sequence, usually three Guanines. Default GGG
- -c,--countonly    Only output the barcode filter summary in the standard output. Default output barcode sequence, cDNA sequence and barcode filter summary to standard output

decode_gene.py
- Note: take either '-f-d[-x-g-t]' (run star then process bam) or '-b' (process bam directly) as the input
- -f, --fastq <input>    Read2 fastq file
- -d, --outdir <dir>     Run star directory
- -x, --STARIDX <starindex>   Reference genome STAR index
- -g, --gtf <file.gtf>   Reference genome annotation
- -t, --thread <int>     Number of threads to use for STAR mapping. Default 1.
- -b, --bam <file.bam>   STAR output bam file, default file name is ```Aligned.toTranscriptome.out.bam```. This bam output requires the parameter ``` --quantMode TranscriptomeSAM GeneCounts``` when running STAR.

decode_quant.py
- -g, --genetable <gene_table>  output of decode_barcode.py
- -b, --barcodetable <barcode_table> output of decode_gene.py
- -u, --usi <usifile> A plain text file containing 2 columns:  USI name and 6bp sequence.

## Examples

```
decode_barcode.py  -i sample_R1.fq -u usi.txt > sample.barcode.tab

decode_gene.py     -f sample_R2.fq -d star_output_dir -t 30 -x star_index -g annotation.gtf > sample.gene.tab

decode_quant.py    -b sample.barcode.tab -g sample.gene.tab -u usi.txt > sample.quant.tab
```
