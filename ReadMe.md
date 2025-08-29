<img src="./media/banner.jpg"  width="500"><br>
# Fastq (Mini) Data Generator

## Overview:
<i>Did you ever need a mini fastq dataset for testing? This will generate one for you!

## Features:
- Reads will show differential expression in your region of choice for (RNA-seq/ChIP-Seq pipeline testing)
- Barcodes can be included (UMI and/or Cell BC)
- Target sequence, read length, target genome, degree of up/down regulation, can all be adjusted

## Flags:
<b>python accioFastq.py --args
| Arg | Feature |
|---|---|
| --fasta | your target sequence |
|--fasta: |your target sequence|
|--outdir: |your output folder name|
|--layout: |single or paired-end reads (SE, PE)|
|--reads: |Number of reads (e.g. 5000)|
|--read-len |Read length (e.g. 75)|
|--umi-len |UMI length (e.g. 10)|
|--cb-len |Cell Barcode length (e.g. 10)|
|--up-region |Up-regulated region/gene nuc range|
|--down-region |Down-regulated region/gene nuc. range|
|--up-name |Up-regulated gene/region name|
|--down-name |Down-regulated gene/region name|
|--up-mult |Up-regulated gene fold multiplier|
|--down-mult |Down regulated gene fold multiplier|
|--seed |Any random number|
|--verbose|Reports generation progress. Recommend: on|

# Outputs: (PE flag)
- 3 Control samples (e.g. ctrl1_R1.fastq.gz / ctrl1_R1.fastq.gz)
- 3 Test samples (e.g. test1_R1.fastq.gz / test1_R1.fastq.gz)
- Bed file with the coordinates and names of the up/down regulated genes
- Metrics files (1 per sample) show the up/down-regulated gene name and the number of reads in control/test<br><br>
![genReport.jpg](./media/genReport.jpg)<br>
![files](./media/files.jpg)

# Usage example:
python accioFastq.py \\<br>
    --fasta chrM.fa \\<br>
    --outdir fastq \\<br>
    --layout PE \\<br>
    --reads 5000 \\<br>
    --read-len 75 \\<br>
    --umi-len 10 \\<br>
    --cb-len 0 \\<br>
    --up-region 3307:4262 \\<br>
    --down-region 9207:9990 \\<br>
    --up-name MT-ND1 \\<br>
    --down-name MT-CO3 \\<br>
    --up-mult 4 \\<br>
    --down-mult 0.05 \\<br>
    --seed 42 \\<br>
    --verbose