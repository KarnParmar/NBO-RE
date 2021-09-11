#
#Karnjit Parmar
#2021
#
#!/usr/bin/env python
"""
2020
Karn Parmar
"""
import sys
import os
import numpy as np

# This script will use the nics.py bond searching algorithm to generate a sigma bond framework for a molecule
# the user can then specify the locations of all the double bonds and lone pairs

"""
THINGS LEFT TO DO:
[ ] Request user input for double bonds, lone pairs or lone valences
[ ] Request user input for gaussian input file where desired $CHOOSE keylist will go
[ ] Print formatted $CHOOSE keylist and insert into the gaussian input file
{ } Collect a $CHOOSE keylist from user specified gaussian output file

"""

atomlist = []
filetext = []
npycoord = np.array([])
carbonlist = []
trimmedcarbonlist = [] #contains only carbon atoms that are bonded to two others
trimmedcarbonbonds = []
bondassignment = []
tempbondlist = []
tempbondassignment = []
tempring = [] # holds carbon numbers of potential rings
counter2 = 0
counter = 0
ginputfile = sys.argv[2]

def call_choose(filename):
    tempstr = []
    chooselist = []
    phrase = " $CHOOSE"
    phrase_loc = 0
    phrase_found = False
    phrase2a = "LONE"
    phrase2a_loc = 0
    phrase2a_found = False
    phrase2_found = False
    phrase2 = "BOND"
    phrase2_loc = 0
    phrase3 = "$END"
    phrase3_loc = 0
    phrase_terminator_found = False
    line2 = ''
             
    with open(filename, 'r') as f:
        for (i, line) in enumerate(f):
            if ((phrase in line) and (phrase_found == False)):
                tempstr = '\n' + line
                phrase_loc = i
                phrase_found = True
                print("1")

                
            # FOUND LONE PAIR    
            if ((i == phrase_loc + 1) and (phrase2a in line) and (phrase_found == True)):
                phrase2a_loc = i
                tempstr = tempstr + line
                phrase2a_found = True
                phrase_found = False     
                print("2")
            # NO LONE PAIRS        
            if ((i == phrase_loc + 1) and (phrase2 in line) and (phrase_found == True) and phrase2a_found == False):
                phrase2_loc = i
                tempstr = tempstr + line
                phrase2_found = True
                phrase_found = False
                print("3:",tempstr,"\n\n")
                
            # CONTINUE AFTER LONE PAIR
            if ((i == phrase2a_loc + 1) and (phrase2 in line) and (phrase_found == False) and phrase2a_found == True):
                phrase2_loc = i
                tempstr = tempstr + line
                phrase2_found = True
                phrase_found = False
                print("4")  
                              
            if ( (phrase_found == True) and (i == phrase_loc+1) and  ((phrase2a_found == False) or (phrase2_found == False)) ):
                phrase_found = False
                #print("TEST")  
                              
            """if((i == phrase_loc + 1) and (phrase2 not in line) (phrase_found == True)):
                phrase_found = False
                tempstr = ''"""
            if ((i > phrase2_loc) and (phrase_terminator_found is False) and (phrase2_found == True)):
                tempstr = tempstr + line
                if(phrase3 in line):
                    phrase_terminator_found = True
                if((i>phrase2_loc+20) and (phrase_terminator_found is False)):
                    phrase_terminator_found = True
                    tempstr = ''
                print("5")
                    
                    
            if(phrase_terminator_found == True):
                phrase_found = False
                phrase_loc = 0
                phrase2_found = False
                phrase2_loc = 0
                phrase_terminator_found = False
        
        rw_ginputfile(tempstr)                    

def distance_finder(atom1,atom2) :
    return (((float(npycoord[atom2][0])-float(npycoord[atom1][0]))**2)+((float(npycoord[atom2][1])-float(npycoord[atom1][1]))**2)+((float(npycoord[atom2][2])-float(npycoord[atom1][2]))**2))**(1/2)

## From nics.py, used to collect atoms and xyz coordinates from each line
def slice_string(string1):
    tempstr = '' 	 	
    finalarr = []
    whitespace = []
	
# Find locations of all whitespaces and record values
    whitespace.append(0) # add a white space for the start
    for i in range(0,len(string1)):
        if(string1[i] == ' '):
            whitespace.append(i+1)  ##MODIFICATION (+1)
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
    if('' in finalarr):
        finalarr.remove('')

    counter=0
    for i in range(0, len(finalarr) - 1):
        finalarr[i] = finalarr[i].replace(" ", "")
        finalarr[i] = finalarr[i].replace('\n', "")
    atomlist.append(finalarr[0])
    return finalarr

# An abridged version of autosearchring() from nics.py. Locates all bonds.
def autosearchring():
    global carbonlist
    global trimmedcarbonlist#contains only carbon atoms that are bonded to two others
    global trimmedcarbonbonds
    global tempbondlist
    global tempbondassignment
    global bondassignment
    global tempring # holds carbon numbers of potential rings
    global counter2
    global counter
    #find first carbon atom from beginning
    for i in range(0, len(atomlist)):
        counter = 0
        if(atomlist[i] == 'C' or atomlist[i] == 'S' or atomlist[i] == 'O' or atomlist[i] == 'N' or atomlist[i] == 'B' or atomlist[i] == 'P' or atomlist[i] == 'H'):
            for j in range(0, len(atomlist)-1):
                if(atomlist[i] == 'S' or atomlist[j] == 'S' or atomlist[i] == 'P' or atomlist[j] == 'P'):
                    if(distance_finder(i, j) < 1.71):
                        counter = counter +  1
                else:
                    if(distance_finder(i, j) < 1.6):
                        counter = counter +  1                
            if(counter < 5):
                carbonlist.append(i)
    print("CARBONLIST", carbonlist)            # DEBUGGING
    #trim carbon list to include only carbons with two others bonded
    for i in range(0, len(carbonlist)):
        tempbondlist = []
        tempbondassignment = []
        counter = 0
        for j in range(0, len(carbonlist)):
            if(atomlist[i] == 'S' or atomlist[i] == 'N' or atomlist[i] == 'O' or atomlist[i] == 'P'):
                if (distance_finder(carbonlist[i], carbonlist[j]) < 1.75):
                    if(i is not j):
                        counter = counter + 1
                        tempbondlist.append(carbonlist[j])        #make an array containing all bonded atoms    
                        tempbondassignment.append('S')
            else:
                if(atomlist[j] == 'S' or atomlist[j] == 'P'):
                    if (distance_finder(carbonlist[i], carbonlist[j]) < 1.75):
                        if(i is not j):
                            counter = counter + 1
                            tempbondlist.append(carbonlist[j]) 
                            tempbondassignment.append('S')   
                else:             
                    if (distance_finder(carbonlist[i], carbonlist[j]) < 1.6):
                        if(i is not j):
                            counter = counter + 1
                            tempbondlist.append(carbonlist[j])  
                            tempbondassignment.append('S')           
        if (counter>0 and counter<4):
            trimmedcarbonlist.append(carbonlist[i])
            trimmedcarbonbonds.append(tempbondlist)
            bondassignment.append(tempbondassignment)

"""
 There are two modes to this scriipt:

1. The $CHOOSE list is generated from an xyz file and the double bonds are specified
2. The $CHOOSE list is generated from a gaussian output file

"""

def find_choose():

    filename = sys.argv[1] ## gaussian log file
    ginputfile = sys.argv[2] ## gaussian input file
    call_choose(sys.argv[1])
    print("in function")

def rw_ginputfile(value):

    with open(ginputfile, 'r') as f:
        content = f.read()

    
    f.close()
    
    
    #print(content)
    
    for i in range(0, len(content)):
        if(content[i] == '$'):
            if(content[i+1] == 'e' and content[i+2] == 'n' and content[i+3] == 'd'):
                linenum = i+4
                print("FOUND TERMINATOR AT ", i+3)
                break;                
    
    templist = str(content)
    print(templist[:linenum])
    tempval = linenum + 2
    noreturn = templist[tempval:]
    templist = templist[0:linenum] + value + templist[linenum+2:]  
    content = str(templist)
    
    
    with open(ginputfile, 'w') as f:
        f.write(content)
         
def choose_xyz():

    global carbonlist
    global trimmedcarbonlist#contains only carbon atoms that are bonded to two others
    global trimmedcarbonbonds
    global tempbondlist
    global tempbondassignment
    global bondassignment
    global tempring # holds carbon numbers of potential rings
    global counter2
    global counter
    global npycoord
    
    ## Open an xyz file
    filename = sys.argv[1] ## XYZ file
    ginputfile = sys.argv[2] ## gaussian input file
    
    with open(filename, 'r') as f:
        atomcount = 0
        start_at_2 = False
        for (i, line) in enumerate(f): 
            templine = line[0]
            if(templine.isdigit()):
                start_at_2 = True
            if(start_at_2 == False):
                filetext.append(line)
            else:
                if(i > 1):
                    filetext.append(line)
            #print("LINE[0]", lines[5])
        for i in range(0, len(filetext)):  
            currentline = slice_string(filetext[i])  
            currentline.remove(currentline[0])
            npycoord = np.append(npycoord, currentline, axis = 0)
        for i in range(0, len(npycoord) - 1):
            npycoord[i] = npycoord[i].replace(" ", "")
            npycoord[i] = npycoord[i].replace('\n', "")  
        atomcount=len(filetext) #2019-11-12
            
    
    npycoord.shape = (atomcount,3)    
    autosearchring()
    phrase1 =  "$nbo"
    phrase2  = "$end"
    linenum = -1
    value = "\nTESTING"
    content = []
    j = 0
    chooselist = '$CHOOSE'
    rw_ginputfile(value)
    """
    with open(ginputfile, 'r') as f:
        for (i,line) in enumerate(f):
            content = f.read()
    
    f.close()
    
    
    #print(content)
    
    for i in range(0, len(content)):
        if(content[i] == '$'):
            if(content[i+1] == 'e' and content[i+2] == 'n' and content[i+3] == 'd'):
                linenum = i+4
                print("FOUND TERMINATOR AT ", i+3)
    
    templist = str(content)
    print(templist[:linenum])
    templist = templist[:linenum] + value + templist[linenum:]  
    content = str(templist)
    
    
    with open(ginputfile, 'w') as f:
        f.write(content)
    """
    
    print(trimmedcarbonlist)
    print(trimmedcarbonbonds)
    print(bondassignment)
    
    doublebondlist = []
    tempval = []
    
    #find the double bonds
    for i in range(0, len(sys.argv)):
        if(sys.argv[i] == '-d'):
            for j in range(i+1, len(sys.argv)):
                tempval = sys.argv[j]
                if (tempval.startswith("-")):
                    break;
                if tempval[0].isdigit:
                    print(tempval[0])
                    doublebondlist.append(sys.argv[j])
    
    
    """
    for i in range(0, len(bondassignment):
        if
    """
    
    
    print(doublebondlist)   
    
    """
    with open(ginputfile, 'r') as f:
        for (i, line) in enumerate(f):
            print("test", j)
            for(j, line) in enumerate(f):
                if phrase2 in line:
                    linenum = j
    
    print(linenum)
    """
    
#choose_xyz()
find_choose()