#!/usr/bin/python

import sys
import getopt
import HTSeq
import re

from collections import defaultdict
def hashofhashes(): return defaultdict(hashofhashes)

inputfile  = ''
usifile    = ''
qscutoff   = 20
countonly  = False
boundary   = 'GGG'

############################################################
# read command line arguments
def usage():
    print "Error!"
    print "Usage: " + sys.argv[0] + " -i|--input <input> -u|--usi <usifile> -q|--qscutoff 20 -b|--boundary GGG -c|--countonly"
    return()

try:
    opts, args = getopt.getopt(sys.argv[1:], 'ci:u:q:b:',
                    ['input=', 'usi=', 'qscutoff=', 'boundary=', 'countonly'])
except getopt.GetoptError:
    usage()
    sys.exit(2)

for opt, arg in opts:
    if opt in ('-i', '--input'):
        inputfile = arg
    elif opt in ('-u', '--usi'):
        usifile = arg
    elif opt in ('-q', '--qscutoff'):
        qscutoff = int(arg)
    elif opt in ('-b', '--boundary'):
        boundary = arg
    elif opt in ('-c', '--countonly'):
        countonly = True
        
if ((inputfile=='') or (usifile=='')):
    usage();
    sys.exit(2)


############################################################
# read USI table, one code per line;
# write into a hash of hash.

usitable = hashofhashes()
usiio    = open(usifile, 'r')
for line in usiio:
    line = line.strip()           # remove leading/tailing whitespaces
    if len(line)==0:
        continue
    if re.match('^#', line):
        continue
    usi, sample = line.split()
    usitable[usi]['before'] = 0
    usitable[usi]['after']  = 0


############################################################
# process fastq

fileio  = HTSeq.FastqReader(inputfile)

bdrylen = len(boundary)
total   = 0
badusi  = 0
badqs   = 0
badbdry = 0
badusibdry = 0
passnum = 0
usidict = hashofhashes()

if not countonly:
    print 'readname\tusi\tumi'

for read in fileio:    
    total += 1

    myusi = read.seq[0:6]
#    myumi = read.seq[6:16]
    myumi = read.seq[6:23]
#    mybdry = read.seq[16:16+bdrylen]
    mybdry = read.seq[23:23+bdrylen]
    myname = read.name
    myname = myname[0:myname.find(" ")]

    if usidict.has_key(myusi):
        usidict[myusi] += 1
    else:
        usidict[myusi] = 1 ##

    if min(read.qual[0:6])<=qscutoff:
        badqs   += 1

    if usitable.has_key(myusi):
        usitable[myusi]['before'] += 1
    else:
        badusi += 1

    if mybdry != boundary:
        badbdry += 1

    if ( (not usitable.has_key(myusi)) and (mybdry != boundary) ):
        badusibdry += 1

    if ( usitable.has_key(myusi)) and (min(read.qual[0:6])>qscutoff and (mybdry==boundary) ):
        passnum += 1
        usitable[myusi]['after'] += 1
        if not countonly:
            print '{}\t{}\t{}'.format(myname, myusi, myumi)


print '#############################################################'
print '# USI   \tTotal\t\tPass\t\tFailed'
print '#------------------------------------------------------------'
for myusi in sorted(usitable):
    mybefore = usitable[myusi]['before']
    myafter  = usitable[myusi]['after']
    print '# {}:\t{}\t{:.1f}%\t{}\t{:.1f}%\t{}\t{:.1f}%'.format(
        myusi,
        mybefore, 100.0*mybefore/total,
        myafter, 100.0*myafter/total,
        mybefore-myafter, 100.0*(mybefore-myafter)/total
    )
print '#------------------------------------------------------------'
print '# Sum:  \t{}\t\t{}\t{:.1f}%\t{}\t{:.1f}%'.format(
    total, passnum, 100.0*passnum/total,
    total-passnum, 100.0*(total-passnum)/total)
print '# Bad QS(<={}):\t\t\t\t\t{}\t{:.1f}%'.format(
    qscutoff, badqs, 100.0*badqs/total)
print '# Bad USI:\t\t\t\t\t{}\t{:.1f}%'.format(
    badusi, 100.0*badusi/total)
print '# Bad Bdry:\t\t\t\t\t{}\t{:.1f}%'.format(
    badbdry, 100.0*badbdry/total)
print '# Bad USI/Bdry:\t\t\t\t\t{}\t{:.1f}%'.format(
    badusibdry, 100.0*badusibdry/total)
print '#------------------------------------------------------------'
print '# Top USIs (>1%)'
for myusi in sorted(usidict):
    myusict = usidict[myusi]
    if usitable.has_key(myusi):
        myusi = myusi + "(*)"
    else:
        myusi = myusi + "   "
    if myusict >= total * 0.01:
        print '#\t{}\t{}\t{:.1f}%'.format(myusi, myusict, 100.0*myusict/total)

print '#############################################################'
