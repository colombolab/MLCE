#!/usr/bin/python

#===========================================================================================================================
#pyBEPPE CHANGELOG ver 1.0.3
#---------------------------------------------------------------------------------------------------------------------------
#ver 1.0 - 07 apr 2014
#Porting of BEPPE 1.1.1 on Python
#
#ver 1.0.1 - 08 feb 2016
#Merging the fork without NumPy
#
#ver 1.0.2 - 01 mar 2016
#Bugfixing (over 1000 residues)
#
#ver 1.0.3 - 01 mar 2016
#Bugfixing (patches)
#
#ver 1.0.4 - 02 nov 2016
#Added a printing flag for MLCE
#===========================================================================================================================

import sys, getopt, copy, py_compile, os
import numpy as np


#Import functions from toolbox (change name!)
from lib import toolbox

#VERSION NUMBER (for header file)
version = '1'
sversion = '0'
ssversion = '3'


#PARAMETERS
threshold = 15  #It's the variable balancing the sensitivity of prediction. If not stated, it's set by default to 15
radius = 6.0	#neighbourhood cutoff for the generation of the contact map, defined in angstroms.


#===========================================================================================================================
#BEPPE 01: import data from arguments, manage input errors, provide help, clear output files.
#===========================================================================================================================

try:
	opts, args = getopt.getopt(sys.argv[1:],"f:v:a:s:e:tmh",['help'])
except getopt.GetoptError:
	toolbox.HELP()
	sys.exit(2)

for opt, arg in opts:
	if opt in ('-h','--help'):
		toolbox.HELP()
		sys.exit(2)
	elif opt in ('-f'):
		pdbfilename = arg
		pdbfileaccess = pdbfilename[:-4]
	elif opt in ('-v'):
		vectorsname = arg
	elif opt in ('-a'):
		valuesname = arg
	elif opt in ('-s'):
		if arg is None:
			threshold = 15. #if no number follows the argument, reset to 15.
		else:
			try:
				threshold = float(arg)
			except ValueError:
				toolbox.ERROR("Threshold is not a valid number",
							  version, sversion, ssversion)
				sys.exit(2)
	elif opt in ('-e'):
		enername = arg
	elif opt in ('-t'):
		topolcontrol = 1
	elif opt in ('-m'):
		mlceprint = 1
	else:
		print "Your argument is invalid"

#check on mandatory arguments: protein structure, eigenvectors and values OR protein structure and energy matrix
try:
	pdbfilename
except NameError:
	pdbfilename = None

try:
	vectorsname
except NameError:
	vectorsname = None

try:
	valuesname
except NameError:
	valuesname = None

if not ((pdbfilename is None and vectorsname is None and valuesname is None ) or (pdbfilename is None and enername is None)):
	pass
else:
	toolbox.ERROR("missing input argument(s). Type -h for help",
				  version, sversion, ssversion)
	sys.exit(2)

try:
	vectorsname, valuesname
except NameError:
	try:
		enername
	except NameError:
		toolbox.ERROR("Missing input argument(s). Type -h for help",
					  version, sversion, ssversion)
		sys.exit(2)

#If arguments are ok, let's start with the header.
toolbox.HEADER(version,sversion,ssversion)

#saves the path of the working directory
bash_call = os.popen('pwd')
path = bash_call.readline().replace('\n','')
bash_call.close()


#Open files I/O
try:
	PDBFILE = open(pdbfilename,'r')
except IOError:
	toolbox.ERROR("Structure file %s not found" % pdbfilename,
				  version, sversion, ssversion)
	sys.exit(2)

#raw copy line per line of PDB file
orig_structure = [i for i in PDBFILE.read().split('\n')]
PDBFILE.close()
#formatted list of PDB data
structure = toolbox.READPDBFILE(pdbfilename)

if vectorsname is not None:
	try:
		#EIGENVECT = open(vectorsname,'r')
		print vectorsname
		eigenvect = np.genfromtxt(vectorsname)
	except IOError:
		toolbox.ERROR("Eigenvectors file %s not found" % vectorsname,
					  version, sversion, ssversion)
		sys.exit(2)


	#eigenvect = []
	#for line in EIGENVECT:
	#	eigenvect.append([i for i in line.split()])
	#EIGENVECT.close()

if valuesname is not None:
	try:
		#EIGENVAL = open(valuesname,'r')
		
		eigenval = np.genfromtxt(valuesname)
	except IOError:
		toolbox.ERROR("Eigenvalues file %s not found" % valuesname,
					  version, sversion, ssversion)
		sys.exit(2)
		
	#eigenval = [i for i in EIGENVAL.read().split('\n')]
	#EIGENVAL.close()

if enername is not None:
	try:
		ENERMAT = open(enername,'r')
	except IOError:
		toolbox.ERROR("Energy matrix file %s not found" % enername,
					  version, sversion, ssversion)
		sys.exit(2)
		
	enermat = [i for i in ENERMAT.read().split('\n')]
	for i in range(0,len(enermat)):
		enermat[i] = enermat[i].split()
	
	enermat.pop(); enermat.pop() #blank lines


if 'topolcontrol' not in vars() or topolcontrol != 1: topolcontrol = 0
if 'mlceprint' not in vars() or mlceprint != 1: mlceprint = 0


#===========================================================================================================================
#BEPPE 02: define the protein length, generate a standard pdbfile indexing the residues from n=1
#===========================================================================================================================

structure, protein_length, span, renumbering = toolbox.PDB_CONVERTER(structure,orig_structure)


#===========================================================================================================================
#BEPPE 03: builds the distance matrix from the coordinates of beta carbon atoms.
#===========================================================================================================================
# outname is created to redirect output file to the wd (ISABEL)
outname = os.getcwd()+"/"+pdbfilename.split('/')[-1][:-4]
distance_matrix = toolbox.DISTANCE_MATRIX(pdbfileaccess,topolcontrol,radius,structure)

#===========================================================================================================================
#BEPPE 04: creates the energy matrix from the first eigenvector, if not provided with argument -e
#===========================================================================================================================
if enername is None:
	energy = toolbox.ENERGY_MATRIX(eigenvect,eigenval)
else:
	energy = np.zeros(shape=(protein_length,protein_length))
	for line in enermat:
		try:
			energy[int(line[0])-1][int(line[1])-1] = float(line[2])
		except IndexError:
			pass

#===========================================================================================================================
#BEPPE 05: generation of the MLCE matrix and first selection of residues to be processed into an epitope prediction
#===========================================================================================================================
epitopes_raw,mlce = toolbox.MLCE(protein_length, threshold, distance_matrix, energy, mlceprint);

#===========================================================================================================================
#BEPPE 06: first selection of the epitopes to be processed into a patch prediction
#===========================================================================================================================
patch_strings = toolbox.PATCHES(epitopes_raw, distance_matrix, protein_length)

#===========================================================================================================================
#BEPPE 07: union of patch strings with residues in common. Generation of full localized patches
#===========================================================================================================================
patches = toolbox.JOIN(patch_strings);

#===========================================================================================================================
#BEPPE 08: calculation of the patches' energy score and sorting of the patches based on the score
#===========================================================================================================================
scored = toolbox.SCORE(patches,mlce);

#===========================================================================================================================
#BEPPE 09: Renumber patches in case of pdb starting from res != 1
#===========================================================================================================================
if(renumbering == 1): scored = toolbox.RENUMBER(span,scored)

#===========================================================================================================================
#BEPPE 10: write output logfile
#===========================================================================================================================
toolbox.WRITE(pdbfileaccess, renumbering, span, scored, structure,
			  version,sversion,ssversion)

#===========================================================================================================================
#BEPPE 11: write output pymol file
#===========================================================================================================================
toolbox.PYMOL(outname,pdbfilename,protein_length,scored)


print " Job finished.\n Output files have been written on:"
print "- beppe_"+pdbfileaccess+".log"
print "- beppe_"+pdbfileaccess+".pml";
if topolcontrol == 1:
	print "- " + pdbfileaccess + "_topology.dat"
print "=========================================================================";
