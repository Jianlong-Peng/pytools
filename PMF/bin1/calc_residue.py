'''
#=============================================================================
#     FileName: calc_residue.py
#         Desc: 
#       Author: jlpeng
#        Email: jlpeng1201@gmail.com
#     HomePage: 
#      Version: 0.0.1
#   LastChange: 2014-08-22 12:44:03
#      History:
#=============================================================================
'''

import sys

def main(argv=sys.argv):
	if len(argv) != 4:
		print "\n  Usage: %s v2013_affinity.log pmf.score output"%argv[0]
		print "\n  v2013_affinity.log: name and -LOG(IC50) are in 2nd and last column"
		print "                      1st line is header"
		print "  pmf.score: name and score are in 1st and last column"
		print "             score will be transformed using `0.00863564 * pmf_score - 5.74134`"
		print "  output: where to save `LOG(IC50,M) - transform(pmf.score)`"
		print ""
		sys.exit(1)
	

	#read v2013_affinity.log
	logIC50 = {}
	inf = open(argv[1],'r')
	line = inf.readline()
	for line in inf:
		line = line.split()
		logIC50[line[1]] = -1.*float(line[-1])
	inf.close()

	#read PMF score
	inf = open(argv[2],'r')
	outf = open(argv[3],'w')
	for line in inf:
		line = line.split()
		try:
			val = logIC50[line[0]] - (0.00863564*float(line[-1])-5.74134)
		except KeyError:
			print line[0]
		else:
			print >>outf, line[0],val
	inf.close()
	outf.close()

main()

