import sys
import pysam

if len(sys.argv) != 5:
        sys.exit('Usage: python sc_atac_library_deconvoluter.py [Input Bam file] [Input Index table] [Output prefix] [Output extension]')

inbam = pysam.Samfile(sys.argv[1],'rb')
inindex = open(sys.argv[2],'r')

libdic = {}
outdic = {}
for line in inindex:
	liner = line.strip().split()
	libdic[liner[0]] = liner[1]
	try:
		outdic[liner[1]]
	except KeyError:
		outdic[liner[1]] = pysam.Samfile(sys.argv[3] + "." + liner[1] + sys.argv[4],'wb',template=inbam)

inindex.close()

readdic = {}
reads = inbam.fetch()
for line in reads:
	try:
		currlib = libdic[line.qname.split(':')[0]]
		outdic[currlib].write(line)	
		try:
			readdic[currlib] += 1
		except KeyError:
			readdic[currlib] = 1
	except KeyError:
		continue

inbam.close()
for read in readdic.keys():
	print read + "\t" + str(readdic[read])
	outdic[read].close()
