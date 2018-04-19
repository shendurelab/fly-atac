# Single cell ATAC-seq Pipeline
NOTE: This pipeline was written for Python2, and so won't work if you are using Python3.

## Step 1: Convert fastq files to deduplicated bam files.
Use wrapper script to go through steps of generating deduplicated bam file.
e.g.
```bash
python sc_atac_fastq2bam.py -R1 READ1 -R2 READ2 -O OUTDIR -P PREFIX -G GENOME
```

## Step 2: Deconvolute bam file into barcodes from different experiments.
Use sc_atac_library_deconvoluter to split out reads associated with each experiment.  
e.g.
```bash
cat 2to4.indextable.txt 6to8.indextable.txt 10to12.indextable.txt | sed -i 's/_P*//g' > master.indextable.txt
bash sc_atac_library_deconvoluter [Input Bam file] master.indextable.txt [Output prefix] [Output extension, e.g. '.nodpus.bam']
```
The `*.indextable.txt` files referenced above can be downloaded from this repository. These files include two columns - one which is a potential cell barcode and one which is an experimental annotation for that barcode. The barcodes themselves are 36bp indices that have the following structure: [8bp Nextera N7 barcode] + [10bp PCR P7 barcode] + [10bp PCR P5 barcode] + [8bp Nextera N5 barcode]. If you run the `sc_atac_fastq2bam.py` script to generate bam files, the format of the read names in the output bam will be [Cell Barcode]:[Read Number]#[An indicator of how many edits were required to identify an acceptable barcode for the current read] (e.g. `TCTCGCGCCTCGTCGTAGTCGTCCTTCGTAATCTTA:192901334#0110`). The index table is used to match these cell barcodes with the experiment they came from.

## Step 3: Determine read depth cutoff for cells and then generate window matrix.  
Wrapper script to count reads, call cells and generate some quick diagnostic plots and then count reads overlapping defined windows for each cell barcode. For the "READCUTOFF" parameter, you can either set a specific read depth value or use "auto" instead (the default), which will determine an appropriate cutoff automatically.  
e.g.
```bash
python sc_atac_bam2matrix.py -B BAMFILE -I master.indextable.txt -O OUTDIR -P PREFIX -C READCUTOFF -W WINDOWBED
```

## Step 4: You can generate additional matrices for other defined windows in the genome (e.g. called peaks).
Given defined windows in the genome, generate a matrix that scores each cell for insertion in each window. 
e.g.
```bash
python sc_atac_window_counter.py [Input Bam file] [Input index table of cells] [Window BED] [Output file] [Include sites with no reads? (True/False)]
```

In this instance, the index table should only list actual cells. Such a file is generated automatically when using `sc_atac_bam2matrix.py` above. The standard format we use for index tables is a two column, tab separated text file that lists the cell barcode in the first column and some sort of experimental annotation for that cell in the second column. For this step in the pipeline, the second column is not required though. This repository contains examples of these tables (2to4.readdepth.cells.indextable.txt, 6to8.readdepth.cells.indextable.txt, and 10to12.readdepth.cells.indextable.txt).