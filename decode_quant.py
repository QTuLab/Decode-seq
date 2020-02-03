#!/usr/bin/python

import re
import sys
import getopt

from collections import defaultdict
def hashofhashes(): return defaultdict(hashofhashes)

bamfile = ''
bcfile  = ''
usifile = ''

def usage():
    print "Error!"
    print "Usage: " + sys.argv[0] + " -g <gene_table> -b <barcode_table> -u <usifile>"
    return()

try:
    opts, args = getopt.getopt(sys.argv[1:], 'g:b:u:',
                               ['genetable=', 'barcodetable=','usi='])
except getopt.GetoptError:
    usage()
    sys.exit(2)

for opt, arg in opts:
    if opt in ('-g', '--genetable'):
        bamfile = arg
    elif opt in ('-b', '--barcodetable'):
        bcfile = arg
    elif opt in ('-u', '--usi'):
        usifile = arg

if ((bamfile=='') or (bcfile=='') or (usifile=='')):
    usage();
    sys.exit(2)

bamio = open(bamfile, 'r')
bcio  = open(bcfile, 'r')
usiio = open(usifile, 'r')


bamtab = {}
for line in bamio:
    if re.match('^readname', line):
        continue
    if re.match('^#', line):
        continue
    readname, gene = line.split()
    bamtab[readname] = gene
bamio.close()


mytab = hashofhashes()
for line in bcio:
    if re.match('^readname', line):
        continue
    if re.match('^#', line):
        continue
    readname, usi, umi = line.split()
    if bamtab.has_key(readname):
        gene = bamtab[readname]
        mytab[gene][usi][umi] = 1
bcio.close()


usitable = {}
for line in usiio:
    line = line.strip()
    if len(line)==0:
        continue
    if re.match('^#', line):
        continue
    usi, sample = line.split()
    usitable[usi] = sample

print "transcript",
for myusi in sorted(usitable):
    print "\t{}".format(usitable[myusi]),
print ""

for trpt in sorted(mytab.keys()):
    print '{}'.format(trpt),
    for myusi in sorted(usitable):
        if mytab[trpt].has_key(myusi):
            umict = len(mytab[trpt][myusi].keys())
        else:
            umict = 0
        print '\t{}'.format(umict),
    print
