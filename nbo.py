#
#Karnjit Parmar
#2021
#
#!/usr/bin/env python

import csv, os, sys
import scipy as np
import pprint as pp

i = 1
j = 1

pi_NBO = []
pi_NBO_temp = []
line_num = 1
bndidx_arr = []
perturbnrg_arr = []
CMO_arr = []
CMO_list = []
nicspizz = []
bondindexes  = [[],[]]
spaces = []
sections = 1
atoms = 0
rows = 9
x=0
m=0
wiberg_idx = np.array([])
perturbnrg_list = []
perturb_min = 1.00
filetext = []
NBO6 = 0
nicspizz_orbs = []

# Read the file contents into a large array
def read_filetext():
    with open(filename, 'r') as f:
        for (i,line) in enumerate(f):
            filetext.append(line)
            
# Find all line numbers with phrase "Wiberg bond index ... "
def call_wiberg():
    phrase = "Wiberg bond index matrix in the NAO basis"
    with open(filename, 'r') as f:
        for (i,line) in enumerate(f):
            if phrase in line:
                bndidx_arr.append(i)

def call_nicspizz():
    phrase = "Principal components of the tensor (ppm) for atom gh("
    phrase_loc = 0
    phrase2 = "This tensor is non-symmetric.  The antisymmetric part will be printed."
    phrase2_loc = 0
    
    with open(filename, 'r') as f:
        for (i, line) in enumerate(f):
            if phrase in line:
                phrase_loc = i
            if phrase2 in line:
                phrase2_loc = i
            if(int(phrase2_loc) == int(phrase_loc) + 1):
                nicspizz.append(phrase_loc)
                phrase_loc = 0
                phrase2_loc = 0
                

def call_CMO():
    phrase = "CMO: NBO Analysis of Canonical Molecular Orbitals"
    with open(filename, 'r') as f:
        for (i, line) in enumerate(f):
            if phrase in line:
                CMO_arr.append(i)
                    

def call_perturbnrg():
    foundphrase = 0
    phrase = "Second Order Perturbation Theory Analysis of Fock Matrix in NBO Basis"
    with open(filename, 'r') as f:
        for (i,line) in enumerate(f):
            if phrase in line:
                foundphrase = 1
                perturbnrg_arr.append(i)

    if (foundphrase == 0):
        phrase = "SECOND ORDER PERTURBATION THEORY ANALYSIS OF FOCK MATRIX IN NBO BASIS"
        global NBO6 
        NBO6 = 1
        print(NBO6)
        with open(filename, 'r') as f:
            print("looking...")
            for (i,line) in enumerate(f):
                if phrase in line:
                    foundphrase = 1
                    perturbnrg_arr.append(i)
    if (foundphrase == 1):
        for (i, line) in enumerate(filetext): # append all lines of the perturbation list into perturbnrg_list
            if (i > int(perturbnrg_arr[-1]) + 7):
                perturbnrg_list.append(line) # finding the line number where the very last perturbation energy list is
                if(line[0] == '\n'):
                    break        

def atom_count():
    atom_count = "NAtoms=   "
    atom_marker = 11
    atomstr = []
    counter = ''
    linerange = 0
    tempstr = ''
    
    with open(filename, 'r') as f:
        for (i,line) in enumerate(f):        
            if atom_count in line:
                tempstr = line
                atomcount1 = slice_string(tempstr)
                atomcount2 = atomcount1[1]
    return(int(atomcount2))

def nicspizz_values():
    
    filetext_temp = []
    temp_line = []
    sum = 0
    initial = []
    nicspizz_orbs.append(0)
    found_value = False
    for j in range(0, len(nicspizz_orbs)-1):
 
        for (i, line) in enumerate(filetext):
            filetext_temp.append(line)
            if(i> int(nicspizz[-1]) + 5):
                temp_line = slice_string(line)
                try:
                    initial = int(temp_line[0][0:-1])
                except:
                    initial = temp_line[0][0:-1]
                if(initial != (int(nicspizz_orbs[j]) + 1) and found_value == True):

                    counter = 0
                    for k in range(0, len(temp_line)):
                        
                        try:
                            temp = float(temp_line[k])
                            
                            counter = counter + 1
                            if(counter == 4):
                                sum = sum + float(temp_line[k])

                                break
                        except:
                            pass
                if(initial == int(nicspizz_orbs[j])):
                    counter = 0
                    charcount = 0
                    chars = False
                    for k in range(0, len(temp_line)):
                        tempchar = temp_line[k][1]
                        teststring = "cr"
                        if(temp_line[k].isalpha()):
                            charcount = charcount + 1
                        if ("cr" in line or "lp" in line):
                            chars = True
                            try:
                                temp = float(temp_line[k])
                                counter = counter + 1
                                if(chars == True):
                                    sum = sum + float(temp_line[5])
                                    break;
    
                            except:
                                pass
                        else:
                            if(chars == False):
                                  print(temp_line[7])
                                  sum = sum + float(temp_line[7])
                                  break;                     
                    found_value = True

                if(initial == int(nicspizz_orbs[j]) + 1 or initial == "-"):
                    print("SUM(",nicspizz_orbs[j],")", sum)
                    found_value = False
                    sum = 0;
                    break;
                    


def print_CMO(orbitals):
    CMOlist = []
    CMOlist_arr_temp = []
    CMOlist_arr = []
    filetext_temp = []
    for (i, line) in enumerate(filetext): # append all lines of the perturbation list into perturbnrg_list
        filetext_temp.append(line)
        if (i > int(CMO_arr[-1]) + 3):
            if(len(line) > 3):
                if(line[2] == 'M' and line[3] == 'O'):
                    slicedstring = slice_string(line)
                    if(slicedstring[2] == ' (occ):'):
                        CMO_list.append(line) # finding the line number where the very last perturbation energy list is
                        CMOlist.append(i)
                        if(line[0] == '\n'):
                            break  
                            
    done = False
    while(done != True):
        for i in range(CMOlist[-1], CMOlist[-1] + 15):
            tempstring = slice_string(filetext_temp[i])

            if(tempstring[2] == ' (vir):'):
                CMOlist.append(i)
                done = True

    for i in range(0,len(CMOlist)-3):
        for j in range(CMOlist[i], CMOlist[i+1]):
            CMOlist_arr_temp.append(filetext_temp[j])
        CMOlist_arr.append(CMOlist_arr_temp)
        CMOlist_arr_temp = []
    
    if(orbitals == 1):
        for i in range(0, len(CMOlist_arr)):
            print(''.join(CMOlist_arr[i]))  
          
    if(orbitals == 2):
        CMOPrint = []
        CMOPrintsum = []
        printed = False
        for i in range(0, len(CMOlist_arr)):
            printed = False
            for j in range(0, len(CMOlist_arr[i])):
                tempstring = slice_string(CMOlist_arr[i][j])
                if(tempstring[4] == ' 2)' and printed == False and tempstring[2] != ' LP'):
                    CMOPrintsum.append(CMOlist_arr[i][0])
                    CMOPrint.append(i)
                    
                    printed = True
        
       # print(''.join(CMOPrint))   
        for k in range(0, len(CMOPrint)-1):
            for l in range(CMOPrint[k], CMOPrint[k+1]):
                print(''.join(CMOlist_arr[l]))
        print("Molecular orbitals likely containing pi-NBO's:\n-------------------------------------------------")
        print(''.join(CMOPrintsum)) 


def print_high_energy_perturb(prnt, energycap):
    alkynewarning = False
    if(int(prnt) == 1):
        if(NBO6 == 1):
            print("                                                           E(2)  E(j)-E(i) F(i,j)")
            print("         Donor NBO (i)           Acceptor NBO (j)        kcal/mol   a.u.    a.u. ")
            print(" ================================================================================")
            for i in range(0, len(perturbnrg_list) - 1):
                perturb_temp = slice_string(perturbnrg_list[i])
                if(perturb_temp[1] == ' BD' or perturb_temp[1] == ' LP'):
                    if(float(perturb_temp[-3]) > float(energycap)):
                                    print(perturbnrg_list[i])
    if(int(prnt) == 2):
        if(NBO6 == 1):
            print("                                                           E(2)  E(j)-E(i) F(i,j)")
            print("         Donor NBO (i)           Acceptor NBO (j)        kcal/mol   a.u.    a.u. ")
            print(" ================================================================================")
            for i in range(0, len(perturbnrg_list) - 1):
                              
                printvalue = False
                pi_NBO_temp = []
                perturb_temp = slice_string(perturbnrg_list[i])
                #print(" |", perturb_temp[i], "| ") 
                if(perturb_temp[1] == ' BD'):
                    if(perturb_temp[3] == ' 2)' or perturb_temp[3] ==  ' 3)'):
                        if(float(perturb_temp[-3]) > float(energycap)):
                                        print(perturbnrg_list[i])
                                        pi_NBO_temp.append(perturb_temp[0])
                                        pi_NBO_temp.append(perturb_temp[8])
                                        pi_NBO.append(pi_NBO_temp)
                if(perturb_temp[1] == ' LP'):
                    if(perturb_temp[4] == ' O' and perturb_temp[3] == ' 2)'):
                        printvalue = True
                    if(perturb_temp[4] == ' N'):
                        printvalue = True
                    if((float(perturb_temp[-3]) > float(energycap)) and (printvalue == True)):
                        print(perturbnrg_list[i])   
                        pi_NBO_temp.append(perturb_temp[0])
                        pi_NBO_temp.append(perturb_temp[6])
                        pi_NBO.append(pi_NBO_temp)
                if(perturb_temp[3] == ' 3)'):
                    alkynewarning = True
                    
                       
            print("Pi NBO's for NBODel method: \n", pi_NBO)                    
            if(alkynewarning):
                print("\nWarning: possible alkyne in molecule, potentially need to use BD ( 3) for conjugated alkyne!\n")
    
def print_perturbnrg(a, b, c):
    tempstr = []
    printcount = 0
    sumperturb = 0
    if(a>b):
        a,b = b, a
    if(NBO6 == 0):
        print("                                                                              E(2)  E(j)-E(i) F(i,j)")
        print("         Donor NBO (i)                     Acceptor NBO (j)                 kcal/mol   a.u.    a.u. ")
        print(" ===================================================================================================")
        for i in range(0, len(perturbnrg_list) - 1):
            perturb_temp = slice_string(perturbnrg_list[i])
            print(perturb_temp[-3])
            sumperturb =+ float(perturb_temp[-3])
            if(perturb_temp[1] != ' BD*('):
                if(perturb_temp[1] != ' CR'):
                    for (i,v) in enumerate(perturb_temp):
                        perturb_temp[i] = v.replace('-', '')
                    if (int(perturb_temp[5]) == a):
                        if (int(perturb_temp[8]) == b):
                            if(float(perturb_temp[-3]) > float(c)):
                                print(perturbnrg_list[i])

                                sumperturb = sumperturb + float(perturb_temp[-3])
                                printcount += 1
            if(perturb_temp[1] == " BD*("):
            #    print("numbers:" + int(perturb_temp[4]) + ' ' + int(perturb_temp[7]))
                if (int(perturb_temp[4]) is a):
                    if (int(perturb_temp[7]) is b):
                        if(float(perturb_temp[-3]) > float(c)):
                            print(perturbnrg_list[i])

                            sumperturb = sumperturb + float(perturb_temp[-3])
                            printcount += 1
        if(printcount == 0):
            print("Error: No bond between atoms or no 2nd order perturbation energies within specified range")
    if(NBO6 == 1):
        print("                                                           E(2)  E(j)-E(i) F(i,j)")
        print("         Donor NBO (i)           Acceptor NBO (j)        kcal/mol   a.u.    a.u. ")
        print(" ================================================================================")
        for i in range(0, len(perturbnrg_list) - 1):
            perturb_temp = slice_string(perturbnrg_list[i])
            if(perturb_temp[1] != ' BD*('):
                if(perturb_temp[1] != ' CR'):
                    if (int(perturb_temp[5]) == a):
                        if (int(perturb_temp[7]) == b):
                            if(float(perturb_temp[-3]) > float(c)):
                                print(perturbnrg_list[i])

                                sumperturb = sumperturb + float(perturb_temp[-3])
                                printcount += 1
            if(perturb_temp[1] == " BD*("):
            #    print("numbers:" + int(perturb_temp[4]) + ' ' + int(perturb_temp[7]))
                if (int(perturb_temp[4]) is a):
                    if (int(perturb_temp[6]) is b):
                        if(float(perturb_temp[-3]) > float(c)):
                            print(perturbnrg_list[i])

                            sumperturb = sumperturb + float(perturb_temp[-3])
                            printcount += 1
        if(printcount == 0):
            print("Error: No bond between atoms or no 2nd order perturbation energies within specified range")
        print("Total perturbation energies: " + str(sumperturb))

        
def print_idx():
    var=0
    tempstr = ' '
    npy_array = np.array([])
    n = 4
	# Find the number of matrices containing data
    if atoms/9 >= int(atoms/9):
        sections = int((atoms/9)+1)
    else:
        sections = atoms/9
    for j in range(atoms):
        for k in range(sections):
            for i in range(len(filetext)):
                q = 4+j
                start = bndidx_arr[-1] + q + (atoms+3)*k
                if i == start:
                    tempstr = filetext[i]
                    unmodstr = slice_string(tempstr)
                    unmodstr.remove(unmodstr[0])
                    unmodstr.remove(unmodstr[0])
                    test = np.array(unmodstr)
                    npy_array = np.append(npy_array,test,axis = 0)
    #print(npy_array)
    npy_array.shape = (atoms,atoms)
    return npy_array


# Find Wiberg bond index between two given atoms
def call_bndidx():
    inputstr = []
    # Run the script until user quits
    while(inputstr != 'q'):
        print("Input two numbers to get the Wiberg bond index ('q' to quit):")
        inputstr = input()
        sliceinput = slice_string(inputstr)
        if(sliceinput[0].isdigit()):
                print(npy_array_global[int(sliceinput[0])-1][int(sliceinput[1])-1])

				
def slice_string(string1):
    tempstr = '' 	 	
    finalarr = []
    whitespace = []
# Find locations of all whitespaces and record values
    whitespace.append(0) # add a white space for the start
    for i in range(0,len(string1)):
        if(string1[i] == ' '):
            whitespace.append(i)
    whitespace.append(len(string1))
# Append all values inbetween whitespaces into an array
# the number of whitespaces determines the number of iterations
    for j in range(0,len(whitespace)-1):
        for k in range(whitespace[j],whitespace[j+1]):
            tempstr += string1[k]
        finalarr.append(tempstr) # Append the string to the array
        tempstr = "" # Clear variable after use 
# Delete all the extra white spaces in the array

    while (' ' in finalarr):
        finalarr.remove(' ')
    #while ('-' in finalarr):
    #    finalarr.remove('-')    
    while ('' in finalarr):
        finalarr.remove('')
    #for (i,v) in enumerate(finalarr):
    #    finalarr[i] = v.replace('-', '')        
    return finalarr

def usage(): 
    print("\n\n")
    print("\t\t\t\t\t\t\t\t<<nbo.py>>")
    print("\t\t\t\t\t\t\t\t----------\n\n")
    print("USAGE:")
    print("------------------------------")
    print("filename [option] [parameters]")
    print("------------------------------\n")
    print("options:\n")
    print("1. Second order perturbation energies (E2PERT) : 2nd order perturbation energies\n ")
    print("\t\t a. Perturbation energies between bonding NBO's between two atoms and all other anti-bonding NBO's over a threshold energy")
    print("\t\t\t\t > nbo.py filename.log perturbation atom1 atom2 threshold")
    print("\t\t\t example > nbo.py file.log perturbation 1 2 5.0\n")
    print("\t\t b. ALL perturbation energies above a threshold (kcal/mol)")
    print("\t\t\t\t > nbo.py filename.log perturbation -h threshold (kcal/mol)")
    print("\t\t\t example > nbo.py file.log perturbation -h 10.0\n")
    print("\t\t c. All perturbation energies above a threshold likely to be related to pi-bonds")    
    print("\t\t\t\t > nbo.py filename.log perturbation -hb threshold (kcal/mol)")
    print("\t\t\t example > nbo.py file.log perturbation -hb 10.0\n")    
    print("\n2. Wiberg bond index    : Wiberg bond indices\n")
    print("\t\t a. Wiberg bond index for specified atoms")
    print("\t\t\t\t > nbo.py filename.log wiberg")
    print("\t\t\t requests user input (atom1, atom2) and produces the Wiberg bond order between indicated atoms\n")
    print("\n3. CMO   : Canonical molecular orbitals. Lists the NBO contributions to molecular orbitals.")
    print("\t\t a. NBO contributions to ALL molecular orbitals")
    print("\t\t\t\t > nbo.py filename.log CMO -a\n")
    print("\t\t b. All molecular orbitals with pi-NBO orbital contributions")
    print("\t\t\t\t > nbo.py filename.log CMO -pi\n")    
      
    print("Karnjit Parmar 2020")
    
# MAIN BODY

if (len(sys.argv) > 1):
    
    filename = os.getcwd() + '/' + sys.argv[1]

    read_filetext()
    atoms = atom_count()
    if(sys.argv[2] == 'CMO'):
        call_CMO()
        if(sys.argv[3] == '-pi'):
            print_CMO(2)
        if(sys.argv[3] == '-a'):
            print_CMO(1)
            
    elif(sys.argv[2] == 'nicspizz'):
        call_nicspizz()
        for i in range(3, len(sys.argv)):
            nicspizz_orbs.append(sys.argv[i])
        nicspizz_values()

  	  	  	
        
                
    elif(sys.argv[2] == 'perturbation'):
        if(sys.argv[3] == '-h'):
            call_perturbnrg()
            print_high_energy_perturb(1, sys.argv[4])
        if(sys.argv[3] == '-hb'):
            call_perturbnrg()
            print_high_energy_perturb(2, sys.argv[4])
        if(sys.argv[3] != '-h' and sys.argv[3] != '-hb'):
            call_perturbnrg()
            print_perturbnrg(int(sys.argv[3]), int(sys.argv[4]), float(sys.argv[5]))
            
        #gives the last value:
    elif(sys.argv[2] == 'wiberg'):
        call_wiberg()
        npy_array_global = np.array([])
        npy_array_global = print_idx()
        call_bndidx()
    else:
        usage()
else:
    usage()
    

print("TESTING")

        