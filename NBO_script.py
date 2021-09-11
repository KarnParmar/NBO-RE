#
#Karnjit Parmar
#2021
#
#!/usr/bin/env python

import glob, os

comfilename = "mynewfile.txt"

def writefile(comfilename):
	charge = 0
	multiplicity = 1



	new_file_path = os.path.join(os.getcwd(), comfilename + ".com")
	comfile = open(new_file_path,'w')
	xyz_path = os.path.join(os.getcwd(), comfilename)
	xyzfile = open(xyz_path, 'r')

	comfile.write("%chk=" + comfilename + ".chk")
	comfile.write('\n%mem=00\n%nproc=32\n')
	comfile.write("# opt M062X/6-31+g(d,p) pop=(full,nbo6read) 	gfoldprint")
	comfile.write('\n\n' + comfilename + '\n\n')
	#write charge + multiplicity
	comfile.write(str(charge) + ' ' + str(multiplicity) + '\n')



	xyzlen = xyzfile.readlines()
	type(xyzlen)
	for i in range(0, len(xyzlen)):
		if (i>1):
			comfile.write(xyzlen[i])

	comfile.write('\n$nbo ARCHIVE PLOT BNDIDX file=' + comfilename + ' CMO NLMO STERIC NRT NRTMEM=1 NRTREF=150 NRTTHR=10 MEMORY=1561449115  $end\n\n')



os.chdir(os.getcwd())
for file in glob.glob("*.xyz"):
	print("writing file:")
	print(file)
	writefile(file)
	
