# Single cell ATAC-seq Pipeline

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
cat 2to4.indextable.txt 6to8.indextable.txt 10to12.indextable.txt | sed -i 's/_P*//g' > master.indextable.txt'
bash sc_atac_library_deconvoluter [Input Bam file] master.indextable.txt [Output prefix] [Output extension, e.g. '.nodpus.bam']
```


## Step 3: Determine read depth cutoff for cells and then generate window matrix.  
Wrapper script to count reads, call cells and generate some quick diagnostic plots and then count reads overlapping defined windows for each cell barcode.  
e.g.
```bash
python sc_atac_bame.matrix.py -B BAMFILE -I master.indextable.txt -O OUTDIR -P PREFIX -C READCUTOFF -W WINDOWBED
```


## Step 4: You can generate additional matrices for other defined windows in the genome (e.g. called peaks).
Given defined windows in the genome, generate a matrix that scores each cell for insertion in each window.  
e.g.
```bash
python sc_atac_window_counter.py [Input Bam file] [Input Index table] [Window BED] [Output file] [Include sites with no reads? (True/False)]
```