import argparse
import os
import subprocess

parser = argparse.ArgumentParser(description='A program to convert NextSeq BCL files to fastq files for scATAC-seq analysis.')
parser.add_argument('-R1','--read1', help='Run directory containing BCL files',dest='read1',required=True)
parser.add_argument('-R2','--read2', help='Run directory containing BCL files',dest='read2',required=True)
parser.add_argument('-O','--outdir', help='Output directory',dest='outdir',required=True)
parser.add_argument('-P','--prefix',help='Output file prefix',dest='prefix',required=True)
parser.add_argument('-G','--genome',help='Bowtie genome',dest='genome',required=True)
args = parser.parse_args()

def submitter(commander):
        """Submits commands directly to the command line and waits for the process to finish."""
        submiting = subprocess.Popen(commander,shell=True)
        submiting.wait()

if args.outdir[-1] != '/':
	args.outdir = args.outdir + '/'

try:
	os.makedirs(args.outdir)
except OSError:
	print 'Outdir already exists...'

print "Fixing barcodes..."
scriptdir = os.path.dirname(os.path.abspath(__file__))
splitter = 'python ' + scriptdir + '/sc_atac_10bpbarcode_split.py -1 ' + args.read1 + ' -2 ' + args.read2 + ' -O1 ' + args.outdir + args.prefix + '.split.1.fq -O2 ' + args.outdir + args.prefix + '.split.2.fq -L ' + args.outdir + args.prefix + '.split.log -Z -X'

submitter(splitter)

print "Trimming adapters..."
trimmer = 'java -Xmx1G -jar ' + scriptdir + '/trimmomatic-0.32.jar PE ' + args.outdir + args.prefix + '.split.1.fq.gz ' + args.outdir + args.prefix + '.split.2.fq.gz ' + args.outdir + args.prefix + '.split.1.trimmed.paired.fastq.gz ' + args.outdir + args.prefix + '.split.1.trimmed.unpaired.fastq.gz ' + args.outdir + args.prefix + '.split.2.trimmed.paired.fastq.gz ' + args.outdir + args.prefix + '.split.2.trimmed.unpaired.fastq.gz ILLUMINACLIP:' + scriptdir + '/NexteraPE-PE.fa:2:30:10:1:true TRAILING:3 SLIDINGWINDOW:4:10 MINLEN:20 2> ' + args.outdir + args.prefix + '.split.trimmomatic.log'
submitter(trimmer)

print "Cleaning up..."
cleaner = 'rm ' + args.outdir + args.prefix + '.split.1.fq.gz; rm ' + args.outdir + args.prefix + '.split.2.fq.gz; rm ' + args.outdir + args.prefix + '.split.1.trimmed.unpaired.fastq.gz; rm ' + args.outdir + args.prefix + '.split.2.trimmed.unpaired.fastq.gz'
submitter(cleaner)

print "Mapping reads..."
mapper = "bowtie2 -p 8 -X 2000 -3 1 -x " + args.genome + " -1 $1.split.1.trimmed.paired.fastq.gz -2 $1.split.2.trimmed.paired.fastq.gz 2> $1.split.bowtie2.log | samtools view -bS - > $1.split.bam; samtools view -h -f3 -F12 -q10 $1.split.bam | grep -v '[0-9]'$'\t'chrM | grep -v '[0-9]'$'\t'chrU | grep -v _CTF_ | grep -v _AMBIG_ | samtools view -Su - | samtools sort -@ 8 - $1.split.q10.sort; samtools index $1.split.q10.sort.bam"
submitter(mapper)

print "Deduplicating reads..."
dedup = "python " + scriptdir + "/sc_atac_true_dedup.py $1.split.q10.sort.bam $1.true.nodups.bam"
submitter(dedup)

