import string
import numpy as np
import math
import operator
import os
import sys
import re #STEFANO: Tests regular expressions in string

def HELP():
	print """
DESCRIPTION:
-------------------------------------------------------------------------

BEPPE processes the non-bonded, intramolecular energetics of a protein
to predict putative antibody and MHC-II binding sites.

BEFORE USING:

In order to use the program you will need:

- a pdb file of atomic coordinates. Make sure it does NOT contain holes
or duplicates in the residue numbering. BEPPE will not fix it for you.

- a set of eigenvectors and eigenvalues after energy decomposition. You
may use the incuded package ISABEL for this task (a version of Amber is
required)

- optionally, you can run predictions based on a previously generated
energy matrix (i.e produced by ISABEL). Include the energy matrix file
with argument -e. If an energy matrix is provided, eigenvectors and
eigenvalues files are no longer mandatory.

ADJUSTING PREDICTION SENSITIVITY:

- The prediction sensitivity may be adjusted at the variation of a
parameter, expressed as the fraction of lowest coupling energies to be
considered for prediction.

- This fraction is set with -s argument. If not stated, it is set to 15
by default. Usually, 10 is used for a tight prediction, 20 or 25 for a
soft one.

PRINTING THE CONTACT MAP TO FILE:

- Add argument -t to the input line.

REFERENCE: Scarabelli et al. (2010) Biophys J

=========================================================================
Option    Filename(example) 	 Type       Description
-------------------------------------------------------------------------
  -f      structure.pdb         Input       Protein structure file: pdb
  -v      EIGENVECT.txt         Input       Energy eigenvectors file
  -a      EIGENVAL.txt          Input       Energy eigenvalues file
  -e      enematrix.dat         Input(opt)  Energy matrix
  -s      15                    Input(opt)  Prediction softness
  -t      xxx_topology.dat      Output      Exported Contact Map
  -m      mlce.dat              Output      MLCE printing
          xxx_beppe.log         Output      Prediction logfile
          xxx_beppe.pml         Ouptut      Visual prediction (Pymol)
-------------------------------------------------------------------------
	"""
	return

def HEADER(version,sversion,ssversion):
	print """
=========================================================================
		
                       (py)B    E    P    P    E
		
            Binding Epitope Prediction from Protein Energetics
		
                                ver. %s.%s.%s
		
-------------------------------------------------------------------------
  Copyright (c) 2010-2016, National Research Council, CNR-ICRM, MILANO
  Guido Scarabelli, Claudio Peri, Riccardo Capelli, Giorgio Colombo.
  
  This work is licensed under a Creative Commons
  Attribution-NonCommercial-ShareAlike 3.0 Unported License.
  
=========================================================================
""" % (version,sversion,ssversion)
	return

def ERROR(errstring,version,sversion,ssversion):
	print """
01 11 1  11
0    1101
|11111011
010000  01
0011100    10
001100001   10
10011100011010
100111    0
01  011
....11   1....
....110110.......
......................
"""
	HEADER(version,sversion,ssversion)
	print "  ===============";
	print "  | I AM ERROR: |";
	print "  ===============";
	print errstring
	print "========================================================================="
	return

def READPDBFILE(pdbfilename):
	INPUT = open(pdbfilename,'r')
	structure = []
			
	for line in INPUT.readlines():
		atom = []                                   #list for a single line
		temp = line
		atom.append(''.join(list(temp[0:6])).replace(' ',''))			#Record name
		if atom[0] != 'ATOM':
			pass
		else:
			atom.append(''.join(list(temp[6:11])).replace(' ',''))      #Atom serial number
			atom.append(''.join(list(temp[12:16])).replace(' ',''))     #Atom name
			atom.append(''.join(list(temp[16])).replace(' ',''))        #AltLoc
			atom.append(''.join(list(temp[17:20])).replace(' ',''))     #Residue name
			atom.append(''.join(list(temp[21])).replace(' ',''))        #Chain Identifier
			atom.append(''.join(list(temp[22:27])).replace(' ',''))     #Sequence number
			atom.append(''.join(list(temp[27])).replace(' ',''))        #iCode - OVER 1000aa fix
			atom.append(''.join(list(temp[30:38])).replace(' ',''))     #x coordinate
			atom.append(''.join(list(temp[38:46])).replace(' ',''))     #y coordinate
			atom.append(''.join(list(temp[46:54])).replace(' ',''))     #z coordinate
			atom.append(''.join(list(temp[54:60])).replace(' ',''))     #Occupancy
			atom.append(''.join(list(temp[60:66])).replace(' ',''))     #Temperature factor
			atom.append(''.join(list(temp[72:76])).replace(' ',''))     #SegID
			atom.append(''.join(list(temp[76:78])).replace(' ',''))     #Element symbol
			atom.append(''.join(list(temp[78:80])).replace(' ',''))     #Charge

		structure.append(atom)
	
	
	return structure


def PDB_CONVERTER(structure,orig_structure):
	protein_length = 0		#actual protein length
	renumbering = 0			#return control value:
							#0: correct numbering
							#1: read renumbered.pdb
	span = 0				#difference in numbering between
							#original PDB and renumbered one
	
	#Extract the ATOM lines
	temp = [line for line in structure if 'ATOM' in line ]
	orig_temp = [line for line in orig_structure if 'ATOM' in line ]
	#extract last residue index
	last = int(temp[len(temp) - 1][6])
	#extract first residue index
	first = int(temp[0][6])

	if first == 1: protein_length = last
	#If the first residue is not 1, generate a new pdbfile with fixed numbering
	else:
		renumbering = 1
		OUTPUT = open("renumbered.pdb",'w')
		span = first - 1
		protein_length = last - first + 1
		for i in range(0,len(temp)):
			num = temp[i][6]
			counter = len(num)
			num2 = str(int(num) - span)
			counternew = len(num2)
			#space correction in PDB
			if counter == counternew:
				string.replace(orig_temp[i][:],' '+num+' ',' '+num2+' ')
			if counter == 2 and counternew == 1:
				string.replace(orig_temp[i][:],'  '+num+'	 ','   '+num2+'	 ')
			if counter == 3 and counternew == 1:
				string.replace(orig_temp[i][:],' '+num+'	 ','   '+num2+'	 ')
			if counter == 4 and counternew == 1:
				string.replace(orig_temp[i][:],num+'	 ','   '+num2+'	 ')
			if counter == 3 and counternew == 2:
				string.replace(orig_temp[i][:],' '+num+'	 ','  '+num2+'	 ')
			if counter == 4 and counternew == 2:
				string.replace(orig_temp[i][:],num+'	 ','  '+num2+'	 ')
			if counter == 4 and counternew == 3:
				string.replace(orig_temp[i][:],num+'	 ',' '+num2+'	 ')
			temp[i][6] = num2
			OUTPUT.write(orig_temp[i][:]+'\n')
		OUTPUT.write('TER\n')
		OUTPUT.close()
	return temp, protein_length, span, renumbering


def DISTANCE_MATRIX(pdbfileaccess, topolcontrol, radius, structure):
	
        CBx = [None]*(int(structure[len(structure)-1][6]) - int(structure[0][6]) +1)		#Coordinates x of Cbeta atom
	CBy = [None]*(int(structure[len(structure)-1][6]) - int(structure[0][6]) +1)		#Coordinates y of Cbeta atom
	CBz = [None]*(int(structure[len(structure)-1][6]) - int(structure[0][6]) +1)		#Coordinates z of Cbeta atom
	
        distance_matrix = np.zeros(shape=(len(CBx),len(CBx)))   #Contact matrix (1/0)
	distance_CB_CB = np.zeros(shape=(len(CBx),len(CBx)))	#Distance matrix

	x_distance = np.zeros(shape=(len(CBx),len(CBx)))		#x distances
	y_distance = np.zeros(shape=(len(CBx),len(CBx)))		#y distances
	z_distance = np.zeros(shape=(len(CBx),len(CBx)))		#z distances
        
	# Check the presence of hydrogens
	control = 0
	
	#browse PDB file and looks for C-beta
	for line in structure:
		atom = line[2]
		resname = line[4]
		if atom == 'H' and resname == 'GLY':
			control = 1
			break
        
	if control == 0: hydrogen = 'CA'
	elif control == 1: hydrogen = 'H'

	#browse PDB file and looks for C-beta
	for line in structure:
		atom = line[2]
		resname = line[4]
                #STEFANO TRANSFORMS THIS INTO resname = all other ******KNOWN, AMBER-TRANSLATED****** AA variants. I know this would have better been done with an array or something, but I don't master Python and I'm in a hurry (not GLY!!!!
		if atom == 'CB' and (resname == 'ALA' or resname == 'CYS' or resname == 'CYX' or resname == 'CYM' or resname == 'ASP' or resname == 'ASH' or resname == 'GLU' or resname == 'GLH' or resname == 'PHE' or resname == 'HIP' or resname == 'HIE' or resname == 'HID' or resname == 'HIS' or resname == 'ILE' or resname == 'LYS' or resname == 'LYN' or resname == 'LEU' or resname == 'MET' or resname == 'ASN' or resname == 'NLN' or resname == 'PRO' or resname == 'HYP' or resname == 'OLP' or resname == 'GLN' or resname == 'ARG' or resname == 'SER' or resname == 'OLS' or resname == 'THR' or resname == 'OLT' or resname == 'VAL' or resname == 'TRP' or resname == 'TYR'):
			num = int(line[6])
			CBx[num-1] = float(line[8])
			CBy[num-1] = float(line[9])
			CBz[num-1] = float(line[10])
		#For glycine, it uses hydrogen coordinates
		elif atom == hydrogen and resname == 'GLY':
			num = int(line[6])
			CBx[num-1] = float(line[8])
			CBy[num-1] = float(line[9])
			CBz[num-1] = float(line[10])
                #Standard Termini when they are isolated residues
                elif atom == 'CH3' and (resname == 'NME' or resname == 'ACE'):
                        num = int(line[6])
                        CBx[num-1] = float(line[8])
                        CBy[num-1] = float(line[9])
                        CBz[num-1] = float(line[10])
                #Template for one ion. always check nomenclature
                elif atom == 'ZN':
                        num = int(line[6])
                        CBx[num-1] = float(line[8])
                        CBy[num-1] = float(line[9])
                        CBz[num-1] = float(line[10])
                #Template that should hopefully work for all sugars. First letter: 0-6,P-Z. Final letter: D,U,A,B
                elif atom == 'C1' and re.match(r"[0-6,P-Z].[A,B,D,U]",resname):
                        num = int(line[6])
                        CBx[num-1] = float(line[8])
                        CBy[num-1] = float(line[9])
                        CBz[num-1] = float(line[10])

                #STEFANO'S extra loops for unknown residues in which none of the atoms
                #Check out format for random ions (ZN), and wildcard criterion for sugars)
                #STEFANO DEBUG
                #print(CBx)
                #print("%s %f" % (resname,CBx))

	#builds the distance matrix. Distances between the same residues are set to 0, calculates the coordinates'difference, then the distance using Pitagora for CB-CB. For in-range distances puts the node distance_matrix[i,h] to 1. Else 0.
	length = len(CBx)
	for i in range(0,length):
		for h in range(0,length):
			if i == h:
				distance_matrix[i,h] = 0
			else:
				try:
					x_distance[i,h] = CBx[h] - CBx[i]
				except TypeError:
					print "ERROR! Distance matrix non complete!"
					if CBx[i] is None:
						print "CBx[%d] = None" % (i)
					elif CBx[h] is None:
						print "CBx[%d] = None" % (h)
					sys.exit()
				y_distance[i,h] = CBy[h] - CBy[i]
				z_distance[i,h] = CBz[h] - CBz[i]

				distance_CB_CB[i,h] = math.sqrt(x_distance[i,h] * x_distance[i,h] + y_distance[i,h] * y_distance[i,h] + z_distance[i,h] * z_distance[i,h])
				if distance_CB_CB[i,h] <= radius: distance_matrix[i,h] = 1
				else: distance_matrix[i,h] = 0

	if topolcontrol == 1:
		OUTPUT = open(pdbfileaccess.split("/")[-1]+"_topology.dat",'w')
		for i in range(0,length):
			for h in range(0,length):
				OUTPUT.write(str(i+1)+"\t"+str(h+1)+"\t"+str(int(distance_matrix[i,h]))+"\n")
			OUTPUT.write("\n")
		OUTPUT.close()

	return distance_matrix

#It computes the energy matrix
def ENERGY_MATRIX(eigenvect, eigenval):
	length = int(np.shape(eigenvect)[0])+1
	counter = 1
	
	energy = np.zeros(shape=(length,length))			#energy matrix
	energy_component = np.zeros(shape=(length,length))  #component of energy matrix
														#associated to an eigenvector
	eigenvector = np.zeros(length)					  #first eigenvector
	
	for i in range(0,np.shape(eigenvect)[0]):
		energy_values = eigenvect[i,:]
		eigenvector[i] = energy_values[0]
		counter += 1
	
	#it builds the energy matrix as sum of (energy coupling)*(eigenvalue)
	for i in range(0,length):
		for j in range(0,length):
			energy_component[i,j] = eigenvector[i] * eigenvector[j] * eigenval[0]
			energy[i,j] = energy[i,j] + energy_component[i,j]

	return energy

#It computes the MLCE matrix
def MLCE(protein_length,threshold,distance_matrix,energy,mlceprint):
	counter=0
	mlce = np.zeros(shape=(protein_length,protein_length))			  #Matrix of local coupling energies
	approx_matrix = []		  #MLCE rounded up
	epitopes_raw = []		   #Selection of energy couplings
	
	#Hadamard product of distance and energy matrices, obtaining MLCE matrix
	for i in range(0,protein_length):
		for h in range(0,protein_length):
			mlce[i,h] = distance_matrix[i,h] * energy[i,h]
			if mlce[i,h] != 0: counter += 1
			approx_matrix.append([i,h,round(mlce[i,h],6)])

	if mlceprint == 1:
		np.savetxt("mlce.dat", mlce)

	approx_matrix.sort(key=lambda sl: (sl[2]))

	#It defines the number of couplings to consider
	limit = counter - int(counter * threshold / 100)

	counter = 0

	for line in approx_matrix:
		if line[2] != 0. and counter > limit: epitopes_raw.append(line)
		elif line[2] != 0. and counter <= limit: counter += 1

	counter = 0
	for line in epitopes_raw:
		epitopes_raw[counter] = line[0]
		counter += 1
	

	return REDUND(epitopes_raw),mlce

#remove redundancies from a list
def REDUND(array):
	return list(set(array))


def PATCHES(epitopes_raw, distance_matrix,protein_length):
	h = 0
	patch_strings = [None] * protein_length

	#gather the selected residues in single patches, based on the neighbouring list.
	length = len(epitopes_raw)
	for i in range(0,length):
		for j in range(0,length):
			if distance_matrix[epitopes_raw[i]][epitopes_raw[j]] == 1:
				if patch_strings[h] is None: patch_strings[h] = [ epitopes_raw[j] ]
				else: patch_strings[h].append(epitopes_raw[j])
		
		if patch_strings[h] is None: patch_strings[h] = [ epitopes_raw[j] ]
		else: patch_strings[h].insert(0,epitopes_raw[i])
		h += 1
	
	patch_strings = patch_strings[:h]
	
	return patch_strings


def JOIN(patch_strings):
	
	merged_strings = []
	while len(patch_strings)>0:
		first, rest = patch_strings[0],patch_strings[1:]
		first = set(first)
			
		lf = -1
		while len(first)>lf:
			lf = len(first)
			rest2 = []
			for r in rest:
				if len(first.intersection(set(r)))>0:
					first |= set(r)
				else:
					rest2.append(r)
			rest = rest2

		merged_strings.append(first)
		patch_strings = rest
	
	for patch in merged_strings:
		patch_strings.append(list(patch))

	return patch_strings


def SCORE(patches, mlce):

	patch_average = [0]*len(patches)	#List of patches' average energy
	join = []						   #List of joined patches with scores

	x = 0
	
	#For every patch, retrieve residues' coupling energy from mlce matrix and
	#calculate average energy
	for residues in patches:
		n_residues = len(residues)
                average = (n_residues * n_residues - n_residues) / 2

		patch_energy = 0
		for i in range(0,n_residues):
			for j in range(0,n_residues):
				patch_energy += float(mlce[residues[i]][residues[j]])
	
		patch_average[x] =   patch_energy/average
		residues.sort()
		#append residues numbers
		join.append(residues)
		#append score to list of residues
		join[x].append(round(patch_average[x],9))
		
		x += 1
	
	#sort patches based on score
	join.sort(key=lambda sl: sl[-1])
	
	#remove the score
	for i in range(0,len(join)):
		join[i].pop()

	return join

def RENUMBER(span, scored):
	
	for i in range(0,len(scored)):
		scored[i] = [val + span for val in scored[i]]

	return scored

def WRITE(pdbfileaccess, renumbering, span, scored, structure,
		  version, sversion, ssversion):
	info = []				   #a list of 3 element lists with index, 3-letters
								#code and 1-letter code
	code = [None]*len(scored)   #list of string with 1-letter sequence of patches
	gaps = []				   #list of position of gaps in patches
        
	# reading PDB data
	for line in structure:
		num = int(line[6])
		if renumbering == 1: num += span
		nam = line[4]
		info.append(str(num)+" "+str(nam))
	
	# deleting duplicates
	info = REDUND(info)

	# cast residue id as int
	temp = []
	for line in info:
		temp.append(line.split())
	for i in range(0,len(temp)):
		temp[i][0] = int(temp[i][0])

	info = temp
	# sorting by residue id
	info.sort(key=lambda sl: (sl[0]))

	# adding 1 letter residue name
        #STEFANO ADDS NLN, LYN, ASH, GLH... requires O-glycosylated
	#X for all other residues/ions
        for i in range(0,len(info)):
		if info[i][1] == 'VAL': info[i].append('V')
		elif info[i][1] == 'LEU': info[i].append('L')
		elif info[i][1] == 'ILE': info[i].append('I')
		elif info[i][1] == 'MET': info[i].append('M')
		elif info[i][1] == 'PHE': info[i].append('F')
		elif info[i][1] == 'ASN': info[i].append('N')
                elif info[i][1] == 'NLN': info[i].append('N')
		elif info[i][1] == 'GLU': info[i].append('E')
                elif info[i][1] == 'GLH': info[i].append('E')
		elif info[i][1] == 'GLN': info[i].append('Q')
		elif info[i][1] == 'ASP': info[i].append('D')
                elif info[i][1] == 'ASH': info[i].append('D')
		elif info[i][1] == 'HIS': info[i].append('H')
		elif info[i][1] == 'HID': info[i].append('H')
		elif info[i][1] == 'HIE': info[i].append('H')
		elif info[i][1] == 'HIP': info[i].append('H')
		elif info[i][1] == 'LYS': info[i].append('K')
                elif info[i][1] == 'LYN': info[i].append('K')
		elif info[i][1] == 'ARG': info[i].append('R')
		elif info[i][1] == 'GLY': info[i].append('G')
		elif info[i][1] == 'ALA': info[i].append('A')
		elif info[i][1] == 'SER': info[i].append('S')
		elif info[i][1] == 'THR': info[i].append('T')
		elif info[i][1] == 'TYR': info[i].append('Y')
		elif info[i][1] == 'TRP': info[i].append('W')
		elif info[i][1] == 'CYS': info[i].append('C')
		elif info[i][1] == 'CYX': info[i].append('C')
		elif info[i][1] == 'CYP': info[i].append('C')
		elif info[i][1] == 'PRO': info[i].append('P')
                else: info[i].append('X')

	# opens log file and writes header
	LOG = open("beppe_"+pdbfileaccess.split("/")[-1]+".log",'w')

	LOG.write("""
=========================================================================
	
                       (py)B    E    P    P    E
	
            Binding Epitope Prediction from Protein Energetics
	
                                ver. %s.%s.%s
	
-------------------------------------------------------------------------
	Copyright (c) 2010-2015p, National Research Council, CNR-ICRM, MILANO
	Guido Scarabelli, Claudio Peri, Riccardo Capelli, Giorgio Colombo.
	
	This work is licensed under a Creative Commons
	Attribution-NonCommercial-ShareAlike 3.0 Unported License.

=========================================================================
""" % (version,sversion,ssversion))

	LOG.write("BEPPE has run successfully!\n")
	LOG.write("Your results will be displayed sorted by energy score.\n\n")

	#prints every patches' residue index, 3-letter and 1-letter code
	length = len(scored)


	#shifts by 1 the residue number (lists start from 0, pdb from 1)
	for i in range(0,length):
		scored[i] = [x+1 for x in scored[i]]

	for i in range(0,length):
		counter = i+1
		LOG.write("PATCH "+str(counter)+"\n")
		for res in scored[i]:
			index = zip(*info)[0].index(res)
			LOG.write(str(info[index][0]) + "\t" + info[index][1] + "\t" + info[index][2] + "\n")
			if code[i] is None:
				code[i] = [ info[index][2] ]
			else:
				code[i].append(info[index][2])

		LOG.write("\n")

	#Finds gaps and inserts - in 1-letter sequence
	for i in range(0,length):
		gaps = []
		for j in range(0,len(code[i])-1):
			if scored[i][j+1] > scored[i][j]+1:
				gaps.append(j+1)
		for j in range(0,len(gaps)):
			code[i].insert(gaps[j],'-')
			gaps = [x + 1 for x in gaps]

	#Writes patches sequences
	LOG.write("\n")
	for i in range(0,length):
		counter = i + 1
		LOG.write("PATCH "+str(counter)+"\n")
		for h in code[i]:
			LOG.write(h)
		LOG.write("\n")
		for h in scored[i]:

			LOG.write(str(h)+" ")
		LOG.write("\n")

	LOG.close()


def PYMOL(pdbfileaccess,pdbfilename,protein_length,scored):
	length = len(scored)
	colors = ['marine',
			  'orange',
			  'forest',
			  'firebrick',
			  'brown',
			  'deeppurple',
			  'deepteal',
			  'olive',
			  'deepsalmon',
			  'smudge',
			  'slate',
			  'acquamarine',
			  'yelloworange',
			  'green',
			  'blue',
			  'sand',
			  'red',
			  'lightmagenta',
			  'lightblue',
			  'lightorange',
			  'palegreen']

	OUTPUT = open('beppe_'+pdbfileaccess.split("/")[-1]+'.pml','w')
	OUTPUT.write("load "+pdbfilename+"\nhide lines\nshow cartoon\nset cartoon_fancy_helices\ncolor gray60, all\n\n")
	counter = 1
        #STEFANO MODIFIES COLOUR COUNTER TO GO BACK TO 0 ONCE IT GETS TO 20
        colorcounter = 0
	for i in range(0,length):
                if colorcounter == 21: colorcounter=0       
		OUTPUT.write("select patch" + str(counter) + ", (i;,")
		for h in range(0,len(scored[i])):
			if h == len(scored[i])-1: OUTPUT.write(str(scored[i][h]))
			else: OUTPUT.write(str(scored[i][h])+",")
                print("Called %d"%i)
		OUTPUT.write(")\ncolor "+colors[colorcounter]+", (i;,")
		for h in range(0,len(scored[i])):
			if h == len(scored[i])-1: OUTPUT.write(str(scored[i][h]))
			else: OUTPUT.write(str(scored[i][h])+",")
		OUTPUT.write(")\n")
		counter += 1
                colorcounter += 1

	OUTPUT.write("deselect")
	OUTPUT.close()
