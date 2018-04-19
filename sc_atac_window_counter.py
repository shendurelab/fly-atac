import sys
import pysam

if len(sys.argv) != 6:
	sys.exit('Usage: python sc_atac_window_counter.py [Input Bam file] [Input Index table] [Window BED] [Output file] [Include sites with no reads? (True/False)]')

inbam = sys.argv[1]
indextable = sys.argv[2]
type1bed = sys.argv[3]
outfile = sys.argv[4]
includezeroes = sys.argv[5]

descer = open(indextable,'r')
cells = [x.strip().split()[0] for x in descer.readlines() if '@' not in x]
descer.close()

cellsdic = {}
for x,cell in enumerate(cells):
	cellsdic[cell] = x

def lister(bedfile):
	currfile = open(bedfile,'r')
	currrecout = [line.strip().split()[0:3] for line in currfile]
	currfile.close()
	return currrecout

print "Building window map..."
rec1list = lister(type1bed)
bamfile = pysam.Samfile(inbam,'rb')

def counter(bedtuple,outsfile,first=False):
	templen = len(cells)
	if first:
		print >> outsfile, "chr\tstart\tend\tannot\t" + "\t".join(cells)
	if includezeroes:
		for rec in bedtuple:
			recname = rec[0] + "_" + rec[1] + "_" + rec[2]
			currcounts = [0]*templen
			reads = bamfile.fetch(rec[0], int(rec[1]), int(rec[2]))
			for read in reads:
				readname = read.qname.split(':')[0]
				try:
					currcounts[cellsdic[readname]] += 1
				except KeyError:
					pass
			print >> outsfile, "\t".join(rec[0:3]) + "\t" + recname + "\t" + "\t".join([str(x) for x in currcounts])

	else:
		for rec in bedtuple:
			recname = rec[0] + "_" + rec[1] + "_" + rec[2]
			currcounts = [0]*templen
			reads = bamfile.fetch(rec[0], int(rec[1]), int(rec[2]))
			for read in reads:
				readname = read.qname.split(':')[0]
				try:
					currcounts[cellsdic[readname]] += 1
				except KeyError:
					pass
			if sum(currcounts) > 0:
				print >> outsfile, "\t".join(rec[0:3]) + "\t" + recname + "\t" + "\t".join([str(x) for x in currcounts])


outmat = open(outfile,'w')
print "Counting window reads for each cell..."
counter(rec1list,outmat,first=True)
outmat.close()
