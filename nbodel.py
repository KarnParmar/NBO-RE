#
#Karnjit Parmar
#2021
#
#!/usr/bin/env python
import glob, os
import sys

comfilename = "mynewfile.txt"

def writefile(comfilename, filetype):
    charge = 0
    multiplicity = 1



    xyz_path = os.path.join(os.getcwd(), comfilename)
    xyzfile = open(xyz_path, 'r')
    if(namefile == True):
        comfilename = sys.argv[3]
    new_file_path = os.path.join(os.getcwd(), comfilename + ".com")
    comfile = open(new_file_path,'w')     
    
    comfile.write("%chk=" + comfilename + ".chk")
    comfile.write('\n%mem=00\n%nproc=32\n')
    comfile.write("# M062X/6-31+G(d,p) pop=nbo6del nosymm 	gfoldprint")
    comfile.write('\n\n' + comfilename + '\n\n')
    #write charge + multiplicity
    comfile.write(str(charge) + ' ' + str(multiplicity) + '\n')
  
  
  
    xyzlen = xyzfile.readlines()
    type(xyzlen)
    for i in range(0, len(xyzlen)):
        if (i>1):
            comfile.write(xyzlen[i])
    
    if(int(filetype) == 1):
        comfile.write('\n$nbo BNDIDX $end\n $del lewis $end\n\n')
    if(int(filetype) == 2 and namefile == False):
        length = (len(sys.argv)-2)/2
        comfile.write('\n$nbo BNDIDX $end\n $del \n delete ' + str(int(length)) + ' elements')
  	
        for i in range(2, len(sys.argv)-1):
            if(i%2 == 0):
                comfile.write('\n   ' + sys.argv[i] + ' ' + sys.argv[i+1])      

        comfile.write('\n$end\n\n')      
    if(int(filetype) == 2 and namefile == True):
        length = (len(sys.argv)-3)/2
        comfile.write('\n$nbo BNDIDX $end\n $del \n delete ' + str(int(length)) + ' elements')
  	
        for i in range(3, len(sys.argv)-1):
            if(i%2 == 0):
                comfile.write('\n   ' + sys.argv[i] + ' ' + sys.argv[i+1])      

        comfile.write('\n$end\n\n') 
        
def usage():
    print("Usage:")
    print("--------------------------------------\n")
    print("-l\t: Delete all non-Lewis type interactions (pi->pi* NBO, etc). Removes all \n\t    stereoelectronic interactions.")
    print("\t\t -Requires no inputfile and will process all .xyz files in the current folder")
    print("\n\t\t example: > nbo.py -l\n\n")
    print("-o\t: Delete all specified NBO orbital interactions.")
    print("\t\t -Requires input orbitals to be specified (Format: donor, acceptor)")
    print("\t\t -Requires no input file and will process all .xyz files in the current folder")
    print("\n\t\t\n example: nbo.py -l 1 3 6 7 12 15\n")
    print("\n\t\t will delete the following donor->acceptor pairs of NBOs: 1->3, 6->7, 12->15\n\n")

noinput = True
namefile = False
if(len(sys.argv) == 1):
    usage()
if(len(sys.argv) > 1):
    if(sys.argv[2] == '-n'):
        namefile = True 
    if(sys.argv[1] == '-o'):
        noinput = False
        os.chdir(os.getcwd())
        if((len(sys.argv)-2)%2 == 0 and len(sys.argv) > 2):
            print("Deleting all non-Lewis type interactions (pi->pi* NBO, etc)")
            for file in glob.glob("*.xyz"):
            	print("writing file:")
            	print(file)
            	writefile(file,2)  
        else:
            print("\n<<Odd number of orbitals/no orbitals specified!>>")
    if(sys.argv[1] == '-l'):     
        noinput = False     
        print("Writing generic file to delete all non-Lewis type interactions (pi->pi* NBO, etc)")
        os.chdir(os.getcwd())
        for file in glob.glob("*.xyz"):
            print("writing file:")
            print(file)
            writefile(file,1)
    if(noinput):
        print("\n\n<<Unrecognized input>> \n\n")
        usage()
