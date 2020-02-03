#!/usr/bin/python

import sys
import os
import getopt
import subprocess
import HTSeq
import re

fqfile    = ''

#STARIDX   = "/home/ysli/db/genome/starindex.r84+.mkdm1.mkdmy.gtf.gsdf/"  # v1
STARIDX   = "/home/ysli/db/genome/starindex.v2.2plus.gtf"  #v2.2

gtf = "/home/ysli/db/annotation/medaka.v2.2plus.gtf"  #v2.2
#gtf = "/home/ysli/data/gene.model/v2.2/trans.blast/all.super.LG.long.ensonly.gtf"  #contain UTR

outdir    = ''
thread    = 1
bamfile   = ''

############################################################
# read command line arguments

def usage():
    print "Error!"
    print "Usage: " + sys.argv[0] + " -f fastq -d outdir [-x STARIDX] [-g gtf] [-t threads] -b bam"
    print "Note: take either '-f-d[-x-g-t]' (run star then process bam) or '-b' (process bam directly) as input"
    return()

try:
    opts, args = getopt.getopt(sys.argv[1:], 'f:x:g:d:t:b:',
                               ['fastq=','STARIDX=','gtf=','outdir=','thread=','bam='])

except getopt.GetoptError:
    usage()
    sys.exit(2)

for opt, arg in opts:
    if opt in ('-f', '--fastq'):
        fqfile = os.path.abspath(arg)
    elif opt in ('-x','--STARIDX'):
        STARIDX = arg
    elif opt in ('-g', '--gtf'):
        gtf = os.path.abspath(arg)
    elif opt in ('-d', '--outdir'):
        outdir = os.path.abspath(arg)
    elif opt in ('-t', '--thread'):
        thread = arg
    elif opt in ('-b', '--bam'):
        bamfile = os.path.abspath(arg)
        
if (fqfile=='' and bamfile==''):   # input: either fastq or bam
    usage();
    sys.exit(2)

if (fqfile!='' and outdir==''):    # fqfile and outdir must co-exist
    usage();
    sys.exit(2)


############################################################
# run STAR

if (fqfile != ''):

    if os.path.exists(outdir):
        print "outdir", outdir, "already exists!"
        sys.exit(2)

    subprocess.call('mkdir %s' %outdir, shell=True)
    subprocess.call('time STAR --runThreadN %s      \
                               --genomeDir %s       \
                               --readFilesIn %s     \
                               --outReadsUnmapped Fastx \
                               --quantMode TranscriptomeSAM GeneCounts \
                               --outSAMtype BAM SortedByCoordinate \
                               1>&2'        # redirect stdout to stderr, so stdout is data only
                    %(thread, STARIDX, fqfile), shell=True, cwd=outdir)

    bamfile = outdir + '/Aligned.toTranscriptome.out.bam'


############################################################
# process bam

if not os.path.exists(bamfile):
   print "cannot find bamfile", bamfile
   sys.exit(2)

statfile = os.path.dirname(os.path.abspath(bamfile)) + '/ReadsPerGene.out.tab'

if not os.path.exists(statfile):
   print "Warning: Can't analyze the mapping stat because ", statfile, " not exist"

bam_reader = HTSeq.BAM_Reader(bamfile)

total = 0

print 'readname\ttranscript'
for align in bam_reader:    
    total += 1
    myread = align.read.name
    mytrpt = align.iv.chrom
    print '{}\t{}'.format(myread, mytrpt)

############################################################
# mapping stat

if not os.path.exists(statfile):
   print "cannot find mapping stat file", statfile
   sys.exit(2)

elif os.path.exists(statfile):   

   list0,list1,list2= [],[],[]

   fi = open(statfile,'r')
   for line in fi.readlines():
       line = line.strip()
       matcher = re.match(r'ERCC',line)
       matchno = re.match(r'N_',line)
       if matchno:
          list0.append(line.split('\t')[0])
          list1.append(int(line.split('\t')[1]))
       elif matcher:
            continue
       else:   
          list2.append(int(line.split('\t')[1]))
     
   total = sum(list2) + sum(list1)

   print '#############################################################'
   print '# mapping stat'
   for i in range(0,4):
       print '#\t{}\t{}\t{:.1f}%'.format(list0[i], list1[i], 100.0*list1[i]/total)
   print '#\t{}\t{}\t{:.1f}%'.format('N_Feature', sum(list2), 100.0*sum(list2)/total)
   print '#\t{}\t{}'.format('Total', total)
   print '#############################################################'

