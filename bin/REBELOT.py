#!/usr/bin/env python3.6
# REBELOT CHANGELOG ver 1.3.2
#
# 1.3.2 - Corrected bugs hampering the python (--py) branch, made pyBEPPE.py "777"
# 1.3.1 - Glycan option optimized, --keep_h option added
# 1.3a  - Glycan option -g added, all the previous modes rolled back to 1.2.3, ported to python3
# 1.2.31- Glycan support being added progressively by Stefano Serapian
# 1.2.3 - Membrane proteins now supported (kinda)
# 1.2.2 - Added support to AmberTools16 (new inputfiles in data)
# 1.2.1 - Small bugfixing and pyBEPPE interface modified
# 1.2   - Added coevolutionary matrix reading (CoCaInE)
# 1.1   - Added MPI support for minimization/Cluster mode
# 1.0   - Initial Release (ex pyISABEL-M)
#
# pyISABEL-M CHANGELOG
#
# 0.6 - Added FF14
# 0.5 - MMPBSA.py (AMBER)
# 0.4 - Multiple eigenvector (ex Domain decomposition)
# 0.3 - First eigenvector
# 0.2 - Migrating to pyBEPPE
# 0.1 - plain porting from ISABEL-M.pl
#
# REBELOT by Riccardo Capelli
# ISABEL-M.pl updated by Claudio Peri
# ISABEL.pl is copyright of Dario Corrada
# ============================================================

# Modules
import sys
import os
import re
import subprocess
import shutil		# Shell interaction
import random					    # Quotes :)
import getopt					    # Command line argument
import copy
import math							# Math
import numpy as np				    # Math
from time import strftime as date   # Time stuff


def main():
    # Global parameters
    workdir = os.getcwd()	   # Working directory

    # Path for REBELOT
    REB_path = subprocess.Popen('which REBELOT.py',
                                stdout=subprocess.PIPE,
                                shell=True).communicate()[0].decode('utf-8').replace('\n', '')


    # Removes the /bin/REBELOT.py
    REB_path = REB_path[:-15]

    # Path for REBELOT
    # STEFANO ADDS CHANGES TO PATH AND NAME OF FORCEFIELD
    GP_path = subprocess.Popen('which gnuplot',
                               stdout=subprocess.PIPE,
                               shell=True).communicate()[0].decode('utf-8').replace('\n', '')
    address = {
        'gnuplot_bin': GP_path,
        'leaprc_rc': REB_path + '/data/leaprc.GLYCAM_06j-1_ff14SB.v18_water.tip3p',
        'minin_rc': REB_path + '/data/min16.in',
        'mmpbsain_rc': REB_path + '/data/mm_pbsa.in',
        'diel_mmpbsain_rc': REB_path + '/data/diel_mm_pbsa.in',
        'mmpbsapyin_rc': REB_path + '/data/mm_pbsa_py.in',
        'tleap_bin': '/home/administrator/amber18/bin/tleap',
        'sander_bin': '/home/administrator/amber18/bin/sander',
        'mmpbsa_bin': '/home/administrator/amber18/bin/mm_pbsa.pl',
        'mmpbsapy_bin': '/home/administrator/amber18/bin/MMPBSA.py',
        #'mmpbsapympi_bin': '/home/administrator/amber18/bin/MMPBSA.py.MPI',
        'sandermpi_bin': '/home/administrator/amber18/bin/sander.MPI',
        'mpirun_bin': '/opt/local/bin/mpirun'
    }
    REB_log = '/REBELOT.log'
    beppe_bin = REB_path + '/bin/pyBEPPE.py'

    ######
    # 01: check arguments, manage input errors, provide help, initialize input files and output folders.
    ######

    # Set optional parameters to zero (default)
    domain = 0
    python = 0
    mpi = 0
    nproc = 1
    coevo = 0
    ext_mat = 0
    dielectric = 0.
    glycans = 0
    keep_h = 0

    try:
        options, remainder = getopt.getopt(sys.argv[1:],
                                           'hm:f:dc:tMD:g',
                                           ['help', 'mode=', 'file=', 'minrange=',
                                            'maxrange=', 'domain', 'cutoff',
                                            'topology', 'mlce', 'py', 'mpi=',
                                            'cluster=', 'coevo=', 'matrix=',
                                            'dielectric=','glycans','keep_h'])
    except getopt.GetoptError as err:
        print(str(err))
        print('Please use the correct arguments, for usage type --help|-h')
        Fine(2)

    for opt, arg in options:
        if opt in ('-h', '--help'):
            Header()
            Help()
            sys.exit(2)
        elif opt in ('-m', '--mode'):
            mode = arg
        elif opt in ('-f', '--file'):
            refpdb = arg
        elif opt in ('-r', '--minrange'):
            minrange = int(arg)
        elif opt in ('-R', '--maxrange'):
            maxrange = int(arg)
        elif opt in ('-d', '--domain'):
            domain = 1
        elif opt in ('-c', '--cutoff'):
            beppe = int(arg)
        elif opt in ('-t', '--topology'):
            topology = '-t'
        elif opt in ('-M', '--mlce'):
            mlceprint = '-m'
        elif opt in ('-D', '--dielectric'):
            print('Dielectric constant' + arg)
            dielectric = float(arg)
        elif opt in ('--py'):
            python = 1
        elif opt in ('--cluster'):
            clustfile = arg
        elif opt in ('--mpi'):
            mpi = 1
            try:
                nproc = int(arg)
            except ValueError:
                Error('Please specify the number of MPI threads.\n')
        elif opt in ('--coevo'):
            coevo = 1
            coevo_path = arg
        elif opt in ('--matrix'):
            ext_mat = 1
            ext_mat_path = arg
        elif opt in ('-g', '--glycans'):
            glycans = 1
        elif opt in ('--keep_h'):
            keep_h  = 1

    # mandatory arguments
    try:
        refpdb
    except NameError:
        Error('Please specify the PDB filename.')

    try:
        mode
    except NameError:
        Error('Please specify the REBELOT mode (SINGLE, MULTIFRAME, CLUSTER or BEPPE).\n')

    if (mode == 'b') and 'res_range' in vars():
        Header()
        Error(
            'With [--range] option enabled it is not possible to start BEPPE.\n')

    # merging maxrange and minrange in res_range
    if 'minrange' in vars() and 'maxrange' in vars():
        res_range = [minrange, maxrange]
    elif 'minrange' in vars() or 'maxrange' in vars():
        Header()
        print('Please specify both minrange and maxrange\n')
        Fine(2)
    else:
        res_range = [0, 0]			# needed to use it in Plots()

    if (mode == 'b') and not 'topology' in vars():
        topology = ' '
    if (mode == 'b') and not 'mlceprint' in vars():
        mlceprint = ' '
    if (mode == 'b') and not 'beppe' in vars():
        beppe = 15	  # Default cutoff
    if (mode == 's' or mode == 'm' or mode == 'c') and 'beppe' in vars():
        Header()
        Error('Cannot use BEPPE cutoffs with mode SINGLE, MULTIFRAME or CLUSTER.\n')
    if (mode == 's' or mode == 'm' or mode == 'c') and 'topology' in vars():
        Header()
        Error('Cannot calculate topology with mode SINGLE, MULTIFRAME or CLUSTER.\n')
    if (mode == 's' or mode == 'm' or mode == 'c') and 'mlceprint' in vars():
        Header()
        Error('Cannot calculate MLCE with mode SINGLE, MULTIFRAME or CLUSTER.\n')
    if (mode == 'c') and not 'clustfile' in vars():
        Error('Cluster file needed for CLUSTER mode')
    if (mode == 's' or mode == 'm' or mode == 'b') and 'clustfile' in vars():
        Error('Cluster file cannot be used in SINGLE, MULTI or BEPPE mode.')
    if (mode == 'm' or mode == 'c') and coevo == 1:
        Error('Coevolutionary matrix cannot be used in MULTI or CLUSTER mode.')
    if (mode == 'm' or mode == 'c') and ext_mat == 1:
        Error('External matrix cannot be used in MULTI or CLUSTER mode.')
    if ext_mat == 1 and coevo == 1:
        Error('Only one external matrix can be specified.')

    # If all the arguments are ok, REBELOT prints the header
    Header()

    # Control of the required files
    for element in address:
        if mpi == 0 and (element != 'mpirun_bin' or element != 'sandermpi_bin'):
            pass
        else:
            if not os.path.exists(address[element]):
                Error('Path ' + address[element] +
                      ' for ' + element + ' not found\n')

    # Single PDB mode
    if mode == 's' or mode == 'b':

        # Reference structure check
        try:
            PDBFILE = open(refpdb, 'r')
        except IOError:
            Error('Structure file ' + refpdb + ' not found\n')

        PDBFILE.close()
        refpdb = workdir + '/' + refpdb
        beppepdb = refpdb

        # WorkDir check
        if not os.path.exists(workdir):
            Error('Path ' + workdir + ' not found\n')
        os.chdir(workdir)
        workdir = os.getcwd()
        folder = workdir + "/REBELOT"
        # check if the REBELOT folder exists and moves the old one
        if os.path.exists(folder):
            mtime = os.stat(folder)[9]
            subprocess.call(["mv", folder, folder + "." + str(mtime)])

        os.mkdir(folder)

        # Sort the min/maxrange
        res_range.sort()

        # Enter the REBELOT folder and set it as destination path
        workdir = folder
        os.chdir(workdir)

        # Log files
        REB_log = workdir + REB_log
        logfile = open(REB_log, 'w')
        logfile.write("REBELOT log file")

        # PDB preprocessing
        if glycans == 0:
            AmberPDB(refpdb, workdir, logfile,keep_h)
            string_for_bonds=""
        currentpdb = workdir + '/snapshot.AMBER.pdb'

        if glycans == 1:
            print('\n/!\\ WARNING /!\\\n')
            print('GLYCANS OPTION ACTIVATED!')
            print('THE -g OPTION IS RECOMMENDED ONLY FOR THE ANALYSIS OF GLYCANS\n')
            string_for_bonds=AmberPDB_Glycans(refpdb, workdir, logfile,keep_h)

        #print(currentpdb)

        beppepdb = currentpdb

        # Prepare AMBER input
        resnumber = AmberConfig(
            currentpdb, folder, logfile, address, dielectric,string_for_bonds)

        if coevo == 0 and ext_mat == 0:
            # Run minimization and MM_GBSA
            AmberJobs(currentpdb, folder, logfile,
                      address, python, resnumber, nproc)
        elif coevo == 1:
            # Check coevolutionary matrix existance
            if not os.path.exists('../' + coevo_path):
                Error('Path ' + coevo_path +
                      ' for coevolutionary matrix not found.')
            # Convert coevolutionary matrix to a readable format for REBELOT
            ConvertCoevoMatrix('../' + coevo_path)
        elif ext_mat == 1:
            if not os.path.exists('../' + ext_mat_path):
                Error('Path ' + ext_mat_path +
                      ' for coevolutionary matrix not found.')
            shutil.copyfile('../' + ext_mat_path, 'ele-vdw.dat')

        # Diagonalize energy matrix
        Diagonalize(currentpdb, folder, logfile,
                    resnumber, 'ele-vdw.dat', address)

        # Produce energy matrix and (optionally) energy profile
        profile_1d, profile_2d = Automan(
            folder, logfile, 'EIGENVECT.txt', 'EIGENVAL.txt', mode, domain)

        # Generate plots with gnuplot
        Plots(folder, logfile, address, res_range, resnumber, mode)
        Summa(folder, logfile, res_range)

    # Remove temporary files
        # Cleansweep(folder)
        if mode == 'b':
            os.chdir(folder)
            Beppe(folder, logfile, beppe, beppe_bin,
                  beppepdb, topology, mlceprint)
        Fine(0)

    # MultiPDB
    if mode == 'm' or mode == 'c':

        # Reference structure check
        try:
            PDBFILE = open(refpdb, 'r')
        except IOError:
            Error('Trajectory file ' + refpdb + ' not found\n')

        # Split the PDB file
        traj = re.split('TER\s*\nENDMDL\s*\n',PDBFILE.read())
        # Remove the last empty list element
        traj.pop()

        PDBFILE.close()

        if mode == 'c':
            try:
                cluster_length = np.genfromtxt(clustfile)
            except IOError:
                Error('Cluster file ' + clustfile + ' not found\n')

        workdir = os.getcwd()

        # WorkDir check
        if not os.path.exists(workdir):
            Error('Path ' + workdir + ' not found\n')
        os.chdir(workdir)
        workdir = os.getcwd()
        workdir += "/REBELOT"
        # check if the REBELOT folder exists and moves the old one
        if os.path.exists(workdir):
            mtime = os.stat(workdir)[9]
            subprocess.call(["mv", workdir, workdir + "." + str(mtime)])

        os.mkdir(workdir)
        
        if glycans == 1:
            print('\n/!\\ WARNING /!\\\n')
            print('GLYCANS OPTION ACTIVATED!')
            print('THE -g OPTION IS RECOMMENDED ONLY FOR THE ANALYSIS OF GLYCANS\n')

        # Set the frame order
        frame_number = 0
        for frame in traj:

            frame_number += 1

            folder = workdir + "/FRAME_" + str(frame_number)
            os.mkdir(folder)
            with open(folder + "/frame_" + str(frame_number) + ".pdb", 'w') as f:
                f.write("%sENDMDL\n" % frame)

        # List of numpy matrices
        multienergy = []

        # Log files
        REB_log = workdir + REB_log
        logfile = open(REB_log, 'w')
        logfile.write("REBELOT log file")

        # Sort the min/maxrange
        res_range.sort()

        # Begin from 1
        for frame in range(1, frame_number+1):
            folder = workdir + "/FRAME_" + str(frame)
            refpdb = folder + "/frame_" + str(frame) + ".pdb"

            # Enter the frame folder and set it as destination path
            os.chdir(folder)

            print("\n\nI- %s Processing frame %d/%d..." %
                  (str(date('%Y-%m-%d %H:%M:%S')), frame, frame_number))

            # Analysis STEFANO OMITS (I removed this for now - RC)
            if glycans == 0:
                AmberPDB(refpdb, folder, logfile,keep_h)
                string_for_bonds=""
            currentpdb = folder + '/snapshot.AMBER.pdb'

            # STEFANO ADDS (I removed this for now - RC)
            if glycans == 1:
                string_for_bonds=AmberPDB_Glycans(refpdb, folder, logfile,keep_h)

            # Prepare AMBER input
            resnumber = AmberConfig(
                currentpdb, folder, logfile, address, dielectric,string_for_bonds)
            print("RESNUMBER: %d" % resnumber)
            # Run minimization and MM_GBSA
            AmberJobs(currentpdb, folder, logfile,
                      address, python, resnumber, nproc)

            # Store energetic information
            multienergy.append(StoreEnergy(folder, resnumber))

        print("\nI- %s averaging electrostatics and vdw...\n" %
              str(date('%Y-%m-%d %H:%M:%S')))

        # Go back to parent directory
        os.chdir(workdir)

        # Compute the average and print it to a file
        if mode == 'm':
            AverageEnergy(workdir, frame_number, resnumber, multienergy)
        elif mode == 'c':
            ClustersEnergy(workdir, frame_number, resnumber,
                           multienergy, cluster_length)
        # Set the first frame as reference structure
        currentpdb = workdir + '/FRAME_1/snapshot.AMBER.pdb'

        # Diagonalize the average matrix
        Diagonalize(currentpdb, workdir, logfile, resnumber,
                    'ele-vdw_average.dat', address)

        profile_1d, profile_2d = Automan(
            workdir, logfile, 'EIGENVECT.txt', 'EIGENVAL.txt', mode, domain)

        # Generate plots with gnuplot
        Plots(workdir, logfile, address, res_range, resnumber, mode)
        Summa(workdir, logfile, res_range)

        # Cleansweep(workdir)

        Fine(0)


def AmberPDB(pdb, workdir, logfile,keep_h):

    print("\nI-", date('%Y-%m-%d %H:%M:%S'),
          "converting pdb format for AMBER...\n")

    logfile.write("\n\n*** AMBER PDB CONVERSION %s" %
                  str(date('%Y-%m-%d %H:%M:%S')))

    # Open the original pdb file and the destination pdb file
    pdb_infile = open(pdb, 'r')
    pdb_outfile = open(workdir + '/snapshot.AMBER.pdb', 'w')

    # Store original pdb data in a list of lines
    pdb_data = pdb_infile.read().split('\n')

    # Close the input file
    pdb_infile.close()

    # Create an empty list
    pdb_content = []

    # Fill the chain if it does not exist
    chainid = 65

    for line in pdb_data:
        # Check if the line contains atom informations
        if 'TER' in line:
            chainid += 1
        elif len(pdb_content) > 0:
            if 'OXT' in pdb_content[-1][3]:
                chainid += 1
        if line[:4] == 'ATOM':
            if line[21] == ' ':
                pdb_content.append([line[0:5],          # [0] Datatype
                                                        # [1] Atom serial number
                                                        line[5:11],
                                                        # [2] (NOT DEFINED IN STANDARD)
                                                        line[11],
                                                        # [3] Atom name
                                                        line[12:16],
                                                        # [4] Alternate location
                                                        line[16],
                                                        # [5] Residue name
                                                        line[17:20],
                                                        # [6] (NOT DEFINED IN STANDARD)
                                                        line[20],
                                                        # [7] Chain ID
                                                        chr(chainid),
                                                        # [8] Residue seq. number
                                                        line[22:26],
                                                        # [9] Code for residues insertion
                                                        line[26],
                                                        line[27:30],		# [10]
                                                        # [11] x coordinate
                                                        line[30:38],
                                                        # [12] y coordinate
                                                        line[38:46],
                                                        # [13] z coordinate
                                                        line[46:54],
                                                        # [14] Occupancy volume
                                                        line[54:60],
                                                        # [15] T-factor
                                                        line[60:66],
                                                        line[66:len(line)]])  # [16] charge, atom type, etc.
            else:
                pdb_content.append([line[0:5],          # [0] Datatype
                                                        # [1] Atom serial number
                                                        line[5:11],
                                                        # [2] (NOT DEFINED IN STANDARD)
                                                        line[11],
                                                        # [3] Atom name
                                                        line[12:16],
                                                        # [4] Alternate location
                                                        line[16],
                                                        # [5] Residue name
                                                        line[17:20],
                                                        # [6] (NOT DEFINED IN STANDARD)
                                                        line[20],
                                                        # [7] Chain ID
                                                        line[21],
                                                        # [8] Residue seq. number
                                                        line[22:26],
                                                        # [9] Code for residues insertion
                                                        line[26],
                                                        line[27:30],		# [10]
                                                        # [11] x coordinate
                                                        line[30:38],
                                                        # [12] y coordinate
                                                        line[38:46],
                                                        # [13] z coordinate
                                                        line[46:54],
                                                        # [14] Occupancy volume
                                                        line[54:60],
                                                        # [15] T-factor
                                                        line[60:66],
                                                        line[66:len(line)]])  # [16] charge, atom type, etc.

    Nterms = []
    Cterms = []

    for i in range(0, len(pdb_content)):
        # Find N-term residues
        if pdb_content[i][3] == ' N  ' or pdb_content[i][3] == 'N   ':
            if pdb_content[i+1][3] == ' H1 ':
                Nterms.append([pdb_content[i][7], pdb_content[i][8]])
        # Find C-term residues
        elif pdb_content[i][3] == ' C  ' or pdb_content[i][3] == 'C   ':
            if pdb_content[i+1][3] == ' O1 ' or pdb_content[i+1][3] == ' OC1' \
                    or pdb_content[i+1][3] == ' O2 ' or pdb_content[i+1][3] == ' OC2':
                Cterms.append([pdb_content[i][7], pdb_content[i][8]])

    # Histidine fix
    histidine_list = []
    prot_list = []

    # Create a list of all the histidines in the system
    for i in range(0, len(pdb_content)):
        if pdb_content[i][5] == 'HIS':
            id_tag = [pdb_content[i][5], pdb_content[i][7], pdb_content[i][8]]
            if not id_tag in histidine_list:
                histidine_list.append(id_tag)
                prot_list.append('')
            if pdb_content[i][3] == ' HD1':
                hist_index = histidine_list.index(id_tag)
                prot_list[hist_index] += 'HD1'
            elif pdb_content[i][3] == ' HE2':
                hist_index = histidine_list.index(id_tag)
                prot_list[hist_index] += 'HE2'

    for i in range(0, len(pdb_content)):
        if pdb_content[i][5] == 'HIS':
            id_tag = [pdb_content[i][5], pdb_content[i][7], pdb_content[i][8]]
            if id_tag in histidine_list:
                hist_index = histidine_list.index(id_tag)
                if prot_list[hist_index] == 'HD1HE2':
                    pdb_content[i][5] = 'HIP'
                elif prot_list[hist_index] == 'HE2':
                    pdb_content[i][5] = 'HIE'
                else:
                    pdb_content[i][5] = 'HID'

    # Isoleucine fix
    for i in range(0, len(pdb_content)):
        if pdb_content[i][5] == 'ILE' and pdb_content[i][3] == ' CD ':
            pdb_content[i][3] = ' CD1'

    # Sulfur atoms
    SG = []
    for i in range(0, len(pdb_content)):
        if pdb_content[i][5] == 'CYS' and pdb_content[i][3] == ' SG ':
            SG.append(pdb_content[i])

    # disulphide bonds
    bridge = []
    S_keys = []
    for i in range(0, len(SG)):
        for j in range(0, len(SG)):
            distance = math.sqrt((float(SG[i][11]) - float(SG[j][11]))**2 +
                                 (float(SG[i][12]) - float(SG[j][12]))**2 +
                                 (float(SG[i][13]) - float(SG[j][13]))**2)
            if distance < 2.15:
                if not [SG[i][7], SG[i][8]] in S_keys and [SG[j][7], SG[j][8]] in S_keys:
                    S_keys.append([SG[i][7], SG[i][8]])
                    S_keys.append([SG[j][7], SG[j][8]])
                    bridge.append(1)
                    bridge.append(1)

    # Cysteine fix (old, FF03 dont have support for N/C-terminal CYX)
    for i in range(0, len(pdb_content)):
        id_tag = [pdb_content[i][7], pdb_content[i][8]]
        if id_tag in S_keys:
            bridge_index = S_keys.index(id_tag)
            if bridge[bridge_index] == 1:
                pdb_content[i][5] = 'CYX'

    # Cysteine fix (CYX removal - see comment above)
    # for i in range(0,len(pdb_content)):
    #	if pdb_content[i][5] == 'CYX':
    #		pdb_content[i][5] = 'CYS'

    # C-term fix
    for i in range(0, len(pdb_content)):
        pattern = [pdb_content[i][7], pdb_content[i][8]]
        if pattern in Cterms and (pdb_content[i][3] == ' O1 ' or pdb_content[i][3] == ' OC1'):
            pdb_content[i][3] = ' O  '
        if pattern in Cterms and (pdb_content[i][3] == ' O2 ' or pdb_content[i][3] == ' OC2'):
            pdb_content[i][3] = ' OXT'

    if keep_h:
        pdb_noH = pdb_content
    else:
    # Hydrogens remove
        pdb_noH = []
        for line in pdb_content:
            if line[3][1] == 'H' or line[3][0] == 'H':
                pass
            else:
                pdb_noH.append(line)

    # atom number, residue number and chain ID sorting
    pdb_sorted = []
    letters = list(map(chr, list(range(65, 91))))
    inc_chain = letters.pop(0)
    inc_atom = 1
    inc_res = 1
    current_chain = pdb_noH[0][7]
    current_res = pdb_noH[0][8]
    prev_resn = 1
    prev_ins = pdb_noH[0][9]

    for line in pdb_noH:
        deep_line = copy.deepcopy(line)

        # Residue number
        newer_res = line[8]
        if newer_res != current_res:
            inc_res += 1
            current_res = newer_res

        # Chain ID
        newer_chain = line[7]
        # Change of chain or last atom
        if newer_chain != current_chain:
            TER = "TER   % 5d     % 4s %s% 4d%s" % \
                (inc_atom, prev_resn, inc_chain, inc_res, prev_ins)
            inc_atom += 1
            inc_chain = letters.pop(0)
            current_chain = newer_chain
            pdb_sorted.append(TER)

        # Right-justified atom and residue number
        deep_line[1] = "% 5d" % inc_atom
        deep_line[8] = "%4d" % inc_res
        deep_line[7] = "%s" % inc_chain

        #if inc_res >= 1000:
        #    deep_line[9] = ''

        # Add new line to pdb
        pdb_sorted.append(deep_line)
        prev_resn = deep_line[5]
        prev_ins = deep_line[9]

        # Increase atom counter
        inc_atom += 1

    # Writing pdb file (snapshot.AMBER.pdb)
    pdb_outfile.write("TITLE     REBELOT\nMODEL	     1\n")
    for line in pdb_sorted:
        outline = ''
        # Formatting for > 10000 atoms
        if 'TER' not in line:
            if int(line[1]) >= 10000:
                line[0] += ''
            else:
                line[0] += ' '
        for col in line:
            outline += col
        pdb_outfile.write("%s\n" % outline)

    # Add a TER and a ENDMDL at the end of the structure
    pdb_outfile.write("TER   % 5d     % 4s %s% 4d%s\n" %
                      (inc_atom, prev_resn, inc_chain, inc_res-1, prev_ins))
    pdb_outfile.write("ENDMDL\n")

    # Close pdb file
    pdb_outfile.close()

    # Update log
    logfile.write("\nI- file <%s> converted in AMBER format\n" % pdb)

def AmberPDB_Glycans(pdb, workdir, logfile,keep_h):

    print("\nI-", date('%Y-%m-%d %H:%M:%S'),
          "converting pdb format for AMBER...\n")

    logfile.write("\n\n*** AMBER PDB CONVERSION %s" %
                  str(date('%Y-%m-%d %H:%M:%S')))


    # Open the original pdb file and the destination pdb file
    pdb_infile = open(pdb, 'r')
    pdb_outfile = open(workdir + '/snapshot.AMBER.pdb', 'w')

    # Store original pdb data in a list of lines
    pdb_data = pdb_infile.read().split('\n')

    # Close the input file
    pdb_infile.close()

    # Create an empty list
    pdb_content = []

    # Fill the chain if it does not exist
    chainid = 65

    for line in pdb_data:
        # Check if the line contains atom informations
        #if 'TER' in line:
         #   chainid += 1
        if len(pdb_content) > 0:
            if 'OXT' in pdb_content[-1][3]:
                chainid += 1
        if line[:4] == 'ATOM':
            if line[21] == ' ':
                pdb_content.append([line[0:5],          # [0] Datatype
                                                        # [1] Atom serial number
                                                        line[5:11],
                                                        # [2] (NOT DEFINED IN STANDARD)
                                                        line[11],
                                                        # [3] Atom name
                                                        line[12:16],
                                                        # [4] Alternate location
                                                        line[16],
                                                        # [5] Residue name
                                                        line[17:20],
                                                        # [6] (NOT DEFINED IN STANDARD)
                                                        line[20],
                                                        # [7] Chain ID
                                                        chr(chainid),
                                                        # [8] Residue seq. number
                                                        line[22:26],
                                                        # [9] Code for residues insertion
                                                        line[26],
                                                        line[27:30],		# [10]
                                                        # [11] x coordinate
                                                        line[30:38],
                                                        # [12] y coordinate
                                                        line[38:46],
                                                        # [13] z coordinate
                                                        line[46:54],
                                                        # [14] Occupancy volume
                                                        line[54:60],
                                                        # [15] T-factor
                                                        line[60:66],
                                                        line[66:len(line)]])  # [16] charge, atom type, etc.
            else:
                pdb_content.append([line[0:5],          # [0] Datatype
                                                        # [1] Atom serial number
                                                        line[5:11],
                                                        # [2] (NOT DEFINED IN STANDARD)
                                                        line[11],
                                                        # [3] Atom name
                                                        line[12:16],
                                                        # [4] Alternate location
                                                        line[16],
                                                        # [5] Residue name
                                                        line[17:20],
                                                        # [6] (NOT DEFINED IN STANDARD)
                                                        line[20],
                                                        # [7] Chain ID
                                                        line[21],
                                                        # [8] Residue seq. number
                                                        line[22:26],
                                                        # [9] Code for residues insertion
                                                        line[26],
                                                        line[27:30],		# [10]
                                                        # [11] x coordinate
                                                        line[30:38],
                                                        # [12] y coordinate
                                                        line[38:46],
                                                        # [13] z coordinate
                                                        line[46:54],
                                                        # [14] Occupancy volume
                                                        line[54:60],
                                                        # [15] T-factor
                                                        line[60:66],
                                                        line[66:len(line)]])  # [16] charge, atom type, etc.

    
   


    Nterms = []
    Cterms = []

    for i in range(0, len(pdb_content)):
        # Find N-term residues
        if pdb_content[i][3] == ' N  ' or pdb_content[i][3] == 'N   ':
            if pdb_content[i+1][3] == ' H1 ':
                Nterms.append([pdb_content[i][7], pdb_content[i][8]])
        # Find C-term residues
        elif pdb_content[i][3] == ' C  ' or pdb_content[i][3] == 'C   ':
            if pdb_content[i+1][3] == ' O1 ' or pdb_content[i+1][3] == ' OC1' \
                    or pdb_content[i+1][3] == ' O2 ' or pdb_content[i+1][3] == ' OC2':
                Cterms.append([pdb_content[i][7], pdb_content[i][8]])

    # Histidine fix
    histidine_list = []
    prot_list = []

    # Create a list of all the histidines in the system
    for i in range(0, len(pdb_content)):
        if pdb_content[i][5] == 'HIS':
            id_tag = [pdb_content[i][5], pdb_content[i][7], pdb_content[i][8]]
            if not id_tag in histidine_list:
                histidine_list.append(id_tag)
                prot_list.append('')
            if pdb_content[i][3] == ' HD1':
                hist_index = histidine_list.index(id_tag)
                prot_list[hist_index] += 'HD1'
            elif pdb_content[i][3] == ' HE2':
                hist_index = histidine_list.index(id_tag)
                prot_list[hist_index] += 'HE2'

    for i in range(0, len(pdb_content)):
        if pdb_content[i][5] == 'HIS':
            id_tag = [pdb_content[i][5], pdb_content[i][7], pdb_content[i][8]]
            if id_tag in histidine_list:
                hist_index = histidine_list.index(id_tag)
                if prot_list[hist_index] == 'HD1HE2':
                    pdb_content[i][5] = 'HIP'
                elif prot_list[hist_index] == 'HE2':
                    pdb_content[i][5] = 'HIE'
                else:
                    pdb_content[i][5] = 'HID'

    # Isoleucine fix
    for i in range(0, len(pdb_content)):
        if pdb_content[i][5] == 'ILE' and pdb_content[i][3] == ' CD ':
            pdb_content[i][3] = ' CD1'

    if keep_h:
        pdb_noH = pdb_content
    else:
    # Hydrogens remove
        pdb_noH = []
        for line in pdb_content:
            if line[3][1] == 'H' or line[3][0] == 'H':
                pass
            else:
                pdb_noH.append(line) 
    # Sulfur atoms
    SG = []
    for i in range(0, len(pdb_content)):
        if ( pdb_content[i][5] == 'CYS' or pdb_content[i][5] == 'CYX') and pdb_content[i][3] == ' SG ':
            SG.append(pdb_content[i])
            

    # disulphide bonds
    bridge = []
    S_keys = []
    for i in range(0, len(SG)):
        for j in range(0,i):
            distance = math.sqrt((float(SG[i][11]) - float(SG[j][11]))**2 +
                                 (float(SG[i][12]) - float(SG[j][12]))**2 +
                                 (float(SG[i][13]) - float(SG[j][13]))**2)
            if distance < 2.15:
                if (not [SG[i][3],SG[i][7], SG[i][8]] in S_keys ) and  (not [SG[j][3],SG[j][7], SG[j][8]] in S_keys ):
                    S_keys.append([SG[i][3],SG[i][7], SG[i][8]])
                    S_keys.append([SG[j][3],SG[j][7], SG[j][8]])
                    bridge.append(1)
                    bridge.append(1)

    # Cysteine fix (old, FF03 dont have support for N/C-terminal CYX)
    for i in range(0, len(pdb_content)):
        id_tag = [pdb_content[i][7], pdb_content[i][8]]
        if id_tag in S_keys:
            bridge_index = S_keys.index(id_tag)
            if bridge[bridge_index] == 1:
                pdb_content[i][5] = 'CYX'
    

    # atom number, residue number and chain ID sorting
    pdb_sorted = []
    letters = list(map(chr, list(range(65, 91))))
    inc_chain = letters.pop(0)
    inc_atom = 1
    inc_res = 1
    current_chain = pdb_noH[0][7]
    current_res = pdb_noH[0][8]
    current_res_type = pdb_noH[0][5]
    prev_resn = 1
    prev_ins = pdb_noH[0][9]
    gly_code   = re.compile('[0-6,P-Z].[D,U,A,B]')

    for line in pdb_noH:
        deep_line = copy.deepcopy(line)

        # Residue number
        newer_res = line[8]
        if newer_res != current_res:
            inc_res += 1
            current_res = newer_res
           
            if gly_code.match(current_res_type):     # Add TER after every glycan
                TER = ["TER","","","","","","","","","","","","","","",""]
                #inc_atom += 1
                pdb_sorted.append(TER)
            current_res_type = line[5]

        # Chain ID
        newer_chain = line[7]
        # Change of chain or last atom
        if newer_chain != current_chain:
            TER = ["TER","","","","","","","","","","","","","","",""]
            #inc_atom += 1
            inc_chain = letters.pop(0)
            current_chain = newer_chain
            pdb_sorted.append(TER)

        # Right-justified atom and residue number
        deep_line[1] = "% 5d" % inc_atom
        deep_line[8] = "%4d" % inc_res
        deep_line[7] = "%s" % inc_chain
        
        #if inc_res >= 1000:
        #    deep_line[8] = deep_line[8].replace(" ","")
        
        # Add new line to pdb
        pdb_sorted.append(deep_line)
        prev_resn = deep_line[5]
        prev_ins = deep_line[9]

        # Increase atom counter
        inc_atom += 1

    # Sulfur atoms with sorted numbering
    SG = []
    for i in range(0, len(pdb_sorted)):
        if ( pdb_sorted[i][5] == 'CYS' or pdb_sorted[i][5] == 'CYX') and pdb_sorted[i][3] == ' SG ':
            SG.append(pdb_sorted[i])
            

    # disulphide bonds
    bridge = []
    S_keys = []
    for i in range(0, len(SG)):
        for j in range(0,i):
            distance = math.sqrt((float(SG[i][11]) - float(SG[j][11]))**2 +
                                 (float(SG[i][12]) - float(SG[j][12]))**2 +
                                 (float(SG[i][13]) - float(SG[j][13]))**2)
            if distance < 2.15:
                if (not [SG[i][3],SG[i][7], SG[i][8]] in S_keys ) and  (not [SG[j][3],SG[j][7], SG[j][8]] in S_keys ):
                    S_keys.append([SG[i][3],SG[i][7], SG[i][8]])
                    S_keys.append([SG[j][3],SG[j][7], SG[j][8]])
                    bridge.append(1)
                    bridge.append(1)

    # C1 atoms
    C1 = []
    for i in range(0, len(pdb_sorted)):
        if (pdb_sorted[i][3] == ' C1 '):
            C1.append(pdb_sorted[i])

    # NLN , O atoms
    oxigen_code= re.compile('\s*O[0-6]\s*')
    others_to_c1 = []
    for i in range(0, len(pdb_sorted)):
        if (pdb_sorted[i][5] == 'NLN' or ( gly_code.match(pdb_sorted[i][5]) and oxigen_code.match(pdb_sorted[i][3]) )):
            others_to_c1.append(pdb_sorted[i])

   
    # Other Bonds
    Bonds_keys=[]
    for c1 in C1:
        old_distance = 1000
        old_atom = []
        for i in range(0, len(others_to_c1)):
            new_distance = math.sqrt((float(c1[11]) - float(others_to_c1[i][11]))**2 +
                                 (float(c1[12]) - float(others_to_c1[i][12]))**2 +
                                 (float(c1[13]) - float(others_to_c1[i][13]))**2)
            if (old_distance > new_distance) and ( c1[8] != others_to_c1[i][8]):
                old_atom = others_to_c1[i]
                old_distance = new_distance
        
        Bonds_keys.append([c1[3],c1[7],c1[8],old_atom[3],old_atom[7],old_atom[8]])    

    # String for bonds
    all_bonds=""
    for i in range(0,len(S_keys),2):
        all_bonds += "bond mol." + str(S_keys[i][2]).replace(" ","") + "." + str(S_keys[i][0]).replace(" ","") + " mol." + str(S_keys[i+1][2]).replace(" ","") + "." + str(S_keys[i+1][0]).replace(" ","") + "\n"
    
    for i in range(0,len(Bonds_keys)):
        all_bonds += "bond mol." +  str(Bonds_keys[i][2]).replace(" ","") + "." + str(Bonds_keys[i][0]).replace(" ","") + " mol." + str(Bonds_keys[i][5]).replace(" ","") + "." + str(Bonds_keys[i][3]).replace(" ","") + "\n"

    # Writing pdb file (snapshot.AMBER.pdb)
    pdb_outfile.write("TITLE     REBELOT\nMODEL	     1\n")
    for line in pdb_sorted:
        outline = ''
        # Formatting for > 10000 atoms
        if 'TER' not in line:
            if int(line[1]) >= 10000:
                line[0] += ''
            else:
                line[0] += ' '
        for col in line:
            outline += col
        pdb_outfile.write("%s\n" % outline)

    # Add a TER and a ENDMDL at the end of the structure
    pdb_outfile.write("TER\n")
    pdb_outfile.write("ENDMDL\n")

    # Close pdb file
    pdb_outfile.close()

    # Update log
    logfile.write("\nI- file <%s> converted in AMBER format\n" % pdb)

    return all_bonds



def AmberConfig(pdb, folder, logfile, address, dielectric, bondstring=""):

    print("I- %s configuring AMBER input files..." %
          str(date('%Y-%m-%d %H:%M:%S')))

    logfile.write("\n\n*** AMBER SETUP %s" % str(date('%Y-%m-%d %H:%M:%S')))

    # Open, read and close the leaprc file
    with open(address['leaprc_rc'], 'r') as f:
        content = f.read()

    # Add the tleap file writing part
    content += '''
mol = loadpdb %s 
%s
saveamberparm  mol prot.prmtop prot.inpcrd
quit
''' % (pdb, bondstring)

    # Write the tleap file
    with open(folder + '/leaprc.ff14', 'w') as f:
        f.write("%s" % content)

    # Open, read and close the pdb file
    with open(pdb, 'r') as f:
        pdb_data = f.read().split('\n')

    res = []

    for line in pdb_data:
        # Check if the line contains atom informations
        if line[:4] == 'ATOM':
            if not line[21:27] in res:		# Chain ID + Residue seq. number
                res.append(line[21:27])

    tot_res = len(res)

    # Editing min.in
    with open(address['minin_rc'], 'r') as f:
        minfile = f.read()

    with open(folder + '/min.in', 'w') as f:
        f.write('%s' % minfile.replace(
            "RES start_number end_number", "RES 1 %d" % tot_res))

    # Editing mm_pbsa.in (no dielectric change)
    if dielectric == 0:
        with open(address['mmpbsain_rc'], 'r') as f:
            mmpbsafile = f.read()

        with open(folder + '/mm_pbsa.in', 'w') as f:
            f.write('%s' % mmpbsafile.replace("start-end", "1-%d" % tot_res))

    # Editing mm_pbsa.in (with dielectric change)
    else:
        with open(address['diel_mmpbsain_rc'], 'r') as f:
            mmpbsafile = f.read()

        with open(folder + '/mm_pbsa.in', 'w') as f:
            f.write('%s' % mmpbsafile.replace("start-end", "1-%d" %
                                              tot_res).replace("_FILL_", str(dielectric)))

    with open(address['mmpbsapyin_rc'], 'r') as f:
        mmpbsafile = f.read()

    with open(folder + '/mm_pbsa_py.in', 'w') as f:
        f.write('%s' % mmpbsafile)

    # Write logfile
    logfile.write("\nI- AMBER config files imported")

    return tot_res


def AmberJobs(pdb, folder, logfile, address, python, resnumber, nproc):

    print("\nI- launching AMBER components...")
    print("\t%s generating topology" % str(date('%Y-%m-%d %H:%M:%S')))

    logfile.write("\n\n*** TOPOLOGY RUN %s\n" % str(date('%Y-%m-%d %H:%M:%S')))

    # Generate topology
    amberlog = subprocess.Popen('%s -s -f leaprc.ff14 2>&1' % address['tleap_bin'],
                                stdout=subprocess.PIPE,
                                shell=True).communicate()[0].decode('utf-8')

    logfile.write("%s" % amberlog)

    # Structure minimization
    # Singlecore
    if (nproc == 1):
        logfile.write("\n\n*** MINIMIZATION RUN %s\n" %
                      str(date('%Y-%m-%d %H:%M:%S')))
        print("\t%s minimization" % str(date('%Y-%m-%d %H:%M:%S')))
        amberlog = subprocess.Popen("%s -O -i min.in -p prot.prmtop -c prot.inpcrd -r prot-min.restrt -o prot.mdout 2>&1"
                                    % address['sander_bin'],
                                    stdout=subprocess.PIPE,
                                    shell=True).communicate()[0].decode('utf-8')
    # Multicore
    else:
        logfile.write("\n\n*** MINIMIZATION RUN (MPI) %s\n" %
                      str(date('%Y-%m-%d %H:%M:%S')))
        print("\t%s minimization (MPI)" % str(date('%Y-%m-%d %H:%M:%S')))
        amberlog = subprocess.Popen("%s -np %d %s -O -i min.in -p prot.prmtop -c prot.inpcrd -r prot-min.restrt -o prot.mdout 2>&1"
                                    % (address['mpirun_bin'], nproc, address['sandermpi_bin']),
                                    stdout=subprocess.PIPE,
                                    shell=True).communicate()[0].decode('utf-8')
    logfile.write("%s" % amberlog)

    # Perl routine
    if python == 0:
        print("\t%s MM-GBSA with energy decomposition" %
              str(date('%Y-%m-%d %H:%M:%S')))

        logfile.write("\n\n*** MM-GBSA RUN %s\n" %
                      str(date('%Y-%m-%d %H:%M:%S')))

        # MM-GBSA
        subprocess.call(["cp", folder + '/prot-min.restrt',
                         folder + '/prot-min_rec.crd.1'])
        amberlog = subprocess.Popen("%s mm_pbsa.in 2>&1" % address['mmpbsa_bin'],
                                    stdout=subprocess.PIPE,
                                    shell=True).communicate()[0].decode('utf-8')

        logfile.write("%s" % amberlog)

        # Read MM-GBSA output
        try:
            with open(folder + '/prot-min_statistics.out', 'r') as f:
                in_content = f.read().split('\n')
        except IOError:
            Error('MM-GBSA output file not found. Check REBELOT.log for details.')

        out_content = ''

        # Parse output file
        for line in in_content:
            record = line.split()
            if len(record) == 52:
                out_content += '%4d  %4d  %10.3f\n' % \
                    (int(record[1]) - 1, int(record[3]) - 1,
                     float(record[14]) + float(record[20]) + float(record[44]))

        with open(folder + '/ele-vdw.dat', 'w') as f:
            f.write('%s' % out_content.replace('-0.000', ' 0.000'))

    # Python routine
    elif python == 1:
        print("\t%s MM-GBSA with energy decomposition (python)" %
              str(date('%Y-%m-%d %H:%M:%S')))

        logfile.write("\n\n*** MM-GBSA RUN (.py) %s\n" %
                      str(date('%Y-%m-%d %H:%M:%S')))
                      
        #STEFANO removes presumable small bug which read unminimised structure
        #(prot.inpcrd)

        #amberlog = subprocess.Popen("%s -O -i mm_pbsa_py.in -cp %s -y %s 2>&1" %
        #                            (address['mmpbsapy_bin'], folder + "/prot.prmtop",
        #                             folder + "/prot.inpcrd"),
        #                            stdout=subprocess.PIPE,
        #                            shell=True).communicate()[0].decode('utf-8')

        #STEFANO FIXES THE ABOVE COMMENTED OUT LINE
        #IF MMPBSA.py was parallelised for decomp too, here is where one would invoke its MPI version with an appropriate branch
        #However, the dirtier workaround is that sander.MPI should be used on the serial MMPBSA.py-generated input.
        amberlog = subprocess.Popen("%s -O -i mm_pbsa_py.in -cp %s -y %s 2>&1" %
                                    (address['mmpbsapy_bin'], folder + "/prot.prmtop",
                                     folder + "/prot-min.restrt"),
                                    stdout=subprocess.PIPE,
                                    shell=True).communicate()[0].decode('utf-8')
        
        logfile.write("%s" % amberlog)

        # Read MM-GBSA output
        try:
            #STEFANO changes to split by newline only, maybe format has changed between versions
            with open(folder + '/FINAL_DECOMP_MMPBSA.dat', 'r') as f:
                in_content = f.read().split('\n')
        except IOError:
            Error('MM-GBSA output file not found. Check REBELOT.log for details.')

        out_content = ''
        # Create an empty matrix
        en_matrix = np.zeros((resnumber, resnumber))
        cont = 0
        # Variable to control output
        ctrl = 0
        for line in in_content:
            record = line.split(',')
            if record[0] == 'Sidechain Energy Decomposition:':
                break
            if record[0] == 'Backbone Energy Decomposition:':
                break
            elif record[0] == 'Resid 1' or record[0] == '':
                pass
            elif len(record) == 20 and ctrl == 0:
                en_matrix[int(record[0][3:]) - 1, int(record[1][3:]) - 1] += \
                    float("%.2f" % float(record[5])) + float("%.2f" % float(record[8])) + float(
                        "%.2f" % float(record[11])) + float("%.2f" % float(record[14]))
            else:
                pass

        for i in range(0, en_matrix.shape[0]):
            for j in range(0, en_matrix.shape[1]):
                out_content += '%4d  %4d  %10.3f\n' % \
                    (i, j, en_matrix[i, j])

        with open(folder + '/ele-vdw.dat', 'w') as f:
            f.write('%s' % out_content.replace('-0.000', ' 0.000'))


def ConvertCoevoMatrix(coevo_path):
    with open(coevo_path, 'r') as f:
        raw_mat = f.read().split('\n')
        # Remove the header
        raw_mat.pop(0)
        length = int(raw_mat[-1].split('\t')[2])

    en_matrix = np.zeros((length, length))

    # Reconstruct everything in a numpy matrix.
    # The -1 at the end is for shift the starting aa to zero.
    for row in raw_mat:
        i = int(row.split('\t')[0]) - 1
        j = int(row.split('\t')[2]) - 1
        en_matrix[i, j] = float(row.split('\t')[4])
        en_matrix[j, i] = float(row.split('\t')[4])

    # Print the matrix in a 'ele-vdw.dat' file
    with open('ele-vdw.dat', 'w') as outf:
        for i in range(0, length):
            for j in range(0, length):
                outf.write("%d\t%d\t%f\n" % (i, j, en_matrix[i, j]))


def Diagonalize(pdb, folder, logfile, resnumber, datafile, address):

    print("\nI- %s diagonalize energy matrix...\n" %
          str(date('%Y-%m-%d %H:%M:%S')))

    logfile.write("\n\n*** DIAGONALIZATION %s\n" %
                  str(date('%Y-%m-%d %H:%M:%S')))

    # Read matrix in list form
    mat_list = np.genfromtxt(datafile)

    # Create a squared matrix nres x nres
    mat_squared = np.zeros((resnumber, resnumber))

    # Convert the matrix in squared form, deleting autointeraction terms
    for el in mat_list:
        if el[0] == el[1]:
            pass
        else:
            mat_squared[int(el[0]), int(el[1])] = el[2]

    # Diagonalize the matrix
    try:
        eigenval, eigenvect = np.linalg.eigh(mat_squared)
    except LinAlgError:
        logfile.write('Eigenvalue computation did not converge!!!\n')
        Error('Eigenvalue computation did not converge!!!')

    # Sort the eigenvalues
    idx = eigenval.argsort()
    eigenval = eigenval[idx]
    eigenvect = eigenvect[:, idx]

    f = open(folder + "/EIGENVAL.txt", 'w')
    g = open(folder + "/EIGENVECT.txt", 'w')

    for i in range(0, resnumber):
        for j in range(0, resnumber):
            g.write("%7s\t" % (str(eigenvect[i, j])))
        f.write("%7s\n" % (str(eigenval[i])))
        g.write("\n")


def Automan(folder, logfile, evecfile, evalfile, mode, domain):

    print("\nI- %s energy eigendecomposition..." %
          str(date('%Y-%m-%d %H:%M:%S')))

    logfile.write("\n*** EIGEN DECOMPOSITION %s\n" %
                  str(date('%Y-%m-%d %H:%M:%S')))

    # Import eigenvectors and eigenvalues as numpy arrays
    evecs = np.genfromtxt(evecfile)
    evals = np.genfromtxt(evalfile)

    resnumber = evals.shape[0]

    sel_eig = []

    if domain == 1:
        logfile.write(
            "\nW- selecting the significant eigenvectors for energy reconstruction\n")

        # Compute a T/F matrix with the significant eigenvector components (absolute values!!!)
        for col_id in range(0, resnumber):
            boolcolumn = np.array(
                [x > np.median(np.abs(evecs[:, col_id])) for x in np.abs(evecs[:, col_id])])
            try:
                # Add new eigenvector
                boolmatrix = np.c_[boolmatrix, boolcolumn]
            except NameError:
                # Otherwise, define matrix as first eigenvector
                boolmatrix = boolcolumn

        # Convert boolmatrix to a 1/0 (int) matrix
        #boolmatrix = boolmatrix.astype(int)

        # Add the first eigenvector to the selection
        sel_eig.append(0)

        # Create a new array that shows which residues are covered (bool) and one
        coverage = boolmatrix[:, 0]
        newcoverage = boolmatrix[:, 0].astype(int)

        # Create a list of remaining eigenvectors
        remaining = list(range(1, resnumber))

        while len(remaining) > 0:
            # If more than 50% of the protein is covered by at least 3 eigenvectors, exit
            # print np.sum(coverage)
            if np.sum(newcoverage >= 3) > resnumber/2.:
                print("\t%s More than 50%% of the sequence is covered by at least 3 eigenvectors!" % str(
                    date('%Y-%m-%d %H:%M:%S')))
                break

            redundancy = []
            for eig in remaining:
                # Compute the redundancy as ((A XOR B) AND A)
                axorb = np.logical_xor(coverage, boolmatrix[:, eig])
                redundancy.append(np.sum(np.logical_and(axorb, coverage)))

            sel_id = redundancy.index(max(redundancy))

            # Compute the added information as ((A XOR B) AND B)
            axorb = np.logical_xor(coverage, boolmatrix[:, remaining[sel_id]])
            newinfo = np.sum(np.logical_and(axorb, coverage))

            if newinfo/float(resnumber) > 0.01:
                coverage = np.logical_or(
                    coverage, boolmatrix[:, remaining[sel_id]])
                newcoverage += axorb
                print("\t%s New eigenvector selected!" %
                      str(date('%Y-%m-%d %H:%M:%S')))
                sel_eig.append(remaining[sel_id])

            remaining.remove(remaining[sel_id])

        selected = ""
        for el in sel_eig:
            selected += "%d " % el

        logfile.write("\n*** SELECTED EIGENS [ %s ]\n" % selected)
        #logfile.write("Coverage: %2.2f\n" % np.sum(coverage)/float(resnumber) )

    else:
        logfile.write(
            "\nW- selecting the first eigenvector for energy reconstruction\n")
        sel_eig.append(0)

    if not mode == 'b':
        profile = Distrib(folder, resnumber, mode, evecs, sel_eig)

    content = EneMatrix(folder, resnumber, evecs, evals, mode, sel_eig)

    if mode == 'b':
        return 'dummy', content
    else:
        return profile, content


def Distrib(folder, resnumber, mode, eigenvectors, sel_eig):

    # Define a threshold for hotspot highlighting
    threshold = math.sqrt(1./float(resnumber))

    components = np.zeros(resnumber)

    for index in sel_eig:
        for i in range(0, resnumber):
            if np.abs(eigenvectors[i, index]) > np.abs(components[i]):
                components[i] = np.abs(eigenvectors[i, index])

    profile = ''

    for i in range(0, components.shape[0]):
        profile += "%d  %.6f\n" % (i, components[i])

    with open(folder + '/enedist.dat', 'w') as f:
        f.write("%s" % profile)

    return profile


def EneMatrix(folder, resnumber, eigenvectors, eigenvalues, mode, sel_eig):

    # Generate a zeros matrix
    enematrix = np.zeros((resnumber, resnumber))

    for index in sel_eig:
        for i in range(0, resnumber):
            for j in range(0, resnumber):
                enematrix[i, j] += eigenvalues[index] * \
                    eigenvectors[i, index] * eigenvectors[j, index]

    out_data = ''

    for i in range(0, resnumber):
        for j in range(0, resnumber):
            out_data += "%d  %d  %.6f\n" % (i+1, j+1, enematrix[i, j])
        out_data += "\n"

    with open(folder + '/enematrix.dat', 'w') as f:
        f.write("%s" % out_data)

    return out_data


def Plots(folder, logfile, address, res_range, resnumber, mode):

    print("\nI- %s generating graphics..." % str(date('%Y-%m-%d %H:%M:%S')))

    logfile.write("\n\n*** PLOTTING %s\n" % str(date('%Y-%m-%d %H:%M:%S')))

    if res_range == [0, 0]:
        # res_range = [1,resnumber]	#Ma tipo fare cosi' no?
        res_range[0] = 1
        if not mode == 'b':
            log = np.genfromtxt(folder + "/enedist.dat").shape[0]
        else:
            log = np.genfromtxt(folder + "/EIGENVAL.txt").shape[0]
        res_range[1] = log

    if not mode == 'b':
        Profilo(folder, logfile, address, res_range, mode)

    Matrice(folder, logfile, address, res_range, mode)


def Profilo(folder, logfile, address, res_range, mode):

    content = np.genfromtxt(folder + "/enedist.dat")[:, 1]
    threshold = np.sqrt(1./float(content.shape[0]))

    # Separate and write the values under/over the threshold
    up = open(folder + "/up", 'w')
    down = open(folder + "/down", 'w')

    i = 1
    for value in content:
        if np.abs(value) < threshold:
            up.write("%d  %.6f\n" % (i, 0.))
            down.write("%d  %.6f\n" % (i, value))
        else:
            up.write("%d  %.6f\n" % (i, value))
            down.write("%d  %.6f\n" % (i, 0.))
        i += 1

    up.close()
    down.close()

    gnuplot = '''#
# *** REBELOT - Renewed and Extendend BEppe Layer On Trajectories ***
#
# Gnuplot script to draw energy distribution
#
set terminal png size 2400, 800
set size ratio 0.33
set output "enedist.png"
set key
set tics out
set xrange[%d:%d]
set xtics rotate
set xtics nomirror
set xtics 10
set mxtics 10
set xlabel "resID"
set yrange [%f:%f+0.1]
set ytics 0.1
set mytics 10
plot    "up" with impulses lw 3 lt 1, \\
"down" with impulses lw 3 lt 9''' % \
        (res_range[0], res_range[1], np.min(content), np.max(content)+0.1)

    with open(folder + '/enedist.gp', 'w') as f:
        f.write("%s" % gnuplot)

    gplot_log = subprocess.Popen("%s enedist.gp 2>&1"
                                 % address['gnuplot_bin'],
                                 stdout=subprocess.PIPE,
                                 shell=True).communicate()[0].decode('utf-8')

    logfile.write("\n%s\n" % gplot_log)

    # Remove datafiles
    os.unlink(folder + "/up")
    os.unlink(folder + "/down")


def Matrice(folder, logfile, address, res_range, mode):

    content = np.genfromtxt(folder + "/enematrix.dat")
    ave = np.sum(content[:, 2])/(float((res_range[1] -
                                        res_range[0])*(res_range[1]-res_range[0])))

    if mode == 'b':
        colorpattern = '35,35,36'
        cbrangemin = ave - 0.5
        cbrangemax = ave
    else:
        colorpattern = '-35,-35,-36'
        cbrangemin = ave + 0.5
        cbrangemax = ave - 0.5

    gnuplot = '''
#
# *** REBELOT - Renewed and Extendend BEppe Layer On Trajectories ***
#
# Gnuplot script to draw energy matrix
#
set terminal png size 2400, 2400
set output "enematrix.png"
set size square
set pm3d map
set palette rgbformulae %s
set cbrange[%f to %f]
set tics out
set xrange[%d:%d]
set yrange[%d:%d]
set xtics 10
set xtics rotate
set ytics 10
set mxtics 10
set mytics 10
splot "enematrix.dat" ''' % \
        (colorpattern, cbrangemin, cbrangemax, res_range[0], res_range[1],
         res_range[0], res_range[1])

    with open(folder + '/enematrix.gp', 'w') as f:
        f.write("%s" % gnuplot)

    gplot_log = subprocess.Popen("%s enematrix.gp 2>&1"
                                 % address['gnuplot_bin'],
                                 stdout=subprocess.PIPE,
                                 shell=True).communicate()[0].decode('utf-8')

    logfile.write("\n%s\n" % gplot_log)


def Summa(folder, logfile, res_range):

    print("\nI- %s calculating total energy...\n" %
          str(date('%Y-%m-%d %H:%M:%S')))

    content = np.genfromtxt(folder + '/enematrix.dat')

    energy_sum = 0.

    for element in content:
        if not element[0] == element[1]:
            if element[0] <= res_range[1] and element[0] >= res_range[0] and \
                    element[1] <= res_range[1] and element[1] >= res_range[0]:
                energy_sum += element[2]

    energy_sum *= 4.184		# kJ/mol -> kcal/mol

    logfile.write('\n\n*** TOTAL ENERGY %.3f kJ/mol\n' % energy_sum)


def Cleansweep(folder):

    print('\nI- %s cleaning temporary files...\n' %
          str(date('%Y-%m-%d %H:%M:%S')))
    try:
        os.unlink(folder + '/blocks.error')
    except OSError:
        pass
    try:
        os.unlink(folder + '/ele-vdw.dat')
    except OSError:
        pass
    try:
        os.unlink(folder + '/enedist.gnuplot')
    except OSError:
        pass
    try:
        os.unlink(folder + '/enematrix.gnuplot')
    except OSError:
        pass
    try:
        os.unlink(folder + '/energy.txt')
    except OSError:
        pass
    try:
        os.unlink(folder + '/energy_blocks1.txt')
    except OSError:
        pass
    try:
        os.unlink(folder + '/leap.log')
    except OSError:
        pass
    try:
        os.unlink(folder + '/leaprc.ff14')
    except OSError:
        pass
    try:
        os.unlink(folder + '/matrice')  # raw energy matrix
    except OSError:
        pass
    try:
        os.unlink(folder + '/mdinfo')
    except OSError:
        pass
    try:
        os.unlink(folder + '/min.in')
    except OSError:
        pass
    try:
        os.unlink(folder + '/mm_pbsa.in')
    except OSError:
        pass
    try:
        os.unlink(folder + '/prot.inpcrd')
    except OSError:
        pass
    try:
        os.unlink(folder + '/prot.mdout')
    except OSError:
        pass
    try:
        os.unlink(folder + '/prot-min.restrt')
    except OSError:
        pass
    try:
        os.unlink(folder + '/prot-min_rec.crd.1')
    except OSError:
        pass
    try:
        os.unlink(folder + '/prot-min_rec.all.out')
    except OSError:
        pass
    try:
        os.unlink(folder + '/prot-min_statistics.in')
    except OSError:
        pass
    try:
        os.unlink(folder + '/SELECTED_EIGENVECTORS.txt')
    except OSError:
        pass


def Beppe(folder, logfile, beppe, beppe_bin, beppepdb, topology, mlceprint):

    print('\nI- %s begin BEPPE job...\n' % str(date('%Y-%m-%d %H:%M:%S')))

    os.chdir(folder)
    beppelog = subprocess.Popen('%s -f %s -e %s -s %s %s %s 2>&1' %
                                (beppe_bin, beppepdb, folder + "/enematrix.dat",
                                 beppe, topology, mlceprint),
                                stdout=subprocess.PIPE,
                                shell=True).communicate()[0].decode('utf-8')

    logfile.write("\n%s\n" % beppelog)


def StoreEnergy(folder, resnumber):

    energymatrix = np.zeros((resnumber, resnumber))

    try:
        raw_data = np.genfromtxt(folder + "/ele-vdw.dat")
    except IOError:
        Error('%s/ele-vdw.dat not found\n' % folder)

    for element in raw_data:
        energymatrix[int(element[0]), int(element[1])] = element[2]

    return energymatrix


def AverageEnergy(workdir, frame_number, resnumber, multienergy):

    ave_mat = np.zeros((resnumber, resnumber))

    for mat in multienergy:
        ave_mat += mat

    ave_mat /= frame_number

    f = open(workdir + "/ele-vdw_average.dat", 'w')

    for i in range(0, resnumber):
        for j in range(0, resnumber):
            avg = "%.3f" % ave_mat[i, j]
            f.write("%5s%5s%12s\n" % (str(i), str(j), avg))

    f.close()


def ClustersEnergy(workdir, frame_number, resnumber, multienergy, cluster_length):

    ave_mat = np.zeros((resnumber, resnumber))

    for i in range(0, len(multienergy)):
        ave_mat += multienergy[i] * cluster_length[i]

    ave_mat /= cluster_length.sum()

    f = open(workdir + "/ele-vdw_average.dat", 'w')

    for i in range(0, resnumber):
        for j in range(0, resnumber):
            avg = "%.3f" % ave_mat[i, j]
            f.write("%5s%5s%12s\n" % (str(i), str(j), avg))

    f.close()


def Error(err):
    error = '''
 _    _   ___  _____  ___   _
| |  | | / _ \|_   _||__ \ | |
| |  | |/ /_\ \ | |     ) || |
| |/\| ||  _  | | |    / / | |
\  /\  /| | | | | |   |_|  |_|
 \/  \/ |_| |_| |_|   (_)  (_)
'''
    Header()
    print("  ===============")
    print("  | I AM ERROR: |")
    print("  ===============")
    print(error)
    print(err)
    sys.exit(2)


def Header():
    header = '''
================================================================================

				   R E B E L O T
			
		Renewed and Extended BEppe Layer On Trajectories
				
				    ver. 1.3.1

--------------------------------------------------------------------------------
                             Copyright (c) 2011-2021, 
        Riccardo Capelli, Claudio Peri, Dario Corrada, Stefano A. Serapian
		
                  This work is licensed under a Creative Commons
            Attribution-NonCommercial-ShareAlike 3.0 Unported License.

================================================================================
'''
    print(header)


def Help():
    help = '''DESCRIPTION:
-------------------------------------------------------------------------

REBELOT is a comprehensive script for the study of the energetic footprint
of your favourite protein, designed to study the determinants of protein
stability as well as the location of putative immunogenic epitopes.

The script operates in three different modes, according to the purpose.
MODE = input argument -m (s|m|b|c)

NOTE: to speed up the minimization part (and if AMBERTools is compiled in
MPI) it is possible to add the flag --mpi followed by the number of cores.
	
==================
| REBELOT SINGLE | (-m s)
==================

The original script: your input PDB is converted in AMBER format and the
energy of nonbonded interactions are decomposed via MM-GBSA.pl from
Ambertools. The interaction energy matrix values include the energy terms
derived from pairwise electrostatic and van der
Waals interactions. Intramolecular type 1-4 interaction are excluded from
the calculation. The energy matrix is reconstructed according to the
chosen method, displaying the determinants of protein stability.
It is also possible to use a pre-calculated raw matrix and skip the GBSA 
evaluation with the flag "--matrix". Furthermore, it is possible to use
a coevolutionary matrix from the CoCaInE[4] analysis using the flag 
"--coevo".
ARGUMENTS FOR MODE REBELOT SINGLE
=========================================================================
Option     example         Type           Description
-------------------------------------------------------------------------

-f(file)   structure.pdb Input          PDB file to be processed.
--matrix   filename      Input(opt)     Raw "ele-vdw.dat" matrix
--coevo    filename      Input(opt)     Coevolutionary potential matrix
-d(domain)               Input(opt)     Alternative energy approximation
                                        via domain decomposition [1].
                                        Default is first eigenvector.
--minrange 5             Input(opt)     Range of residues for which
--maxrange 30                           interaction energy matrix will be
                                        calculated (default: all)
-g                       Input(opt)     Add terms for glycans treatment
-------------------------------------------------------------------------
SYNOPSYS: REBELOT.py  -m s -f structure.pdb  --minrange 1 --maxrange 100
	
======================
| REBELOT MULTIFRAME | (-m m)
======================
The multiframe-variant of REBELOT is designed to operate on a small PDB 
trajectory producing a regular MM-PBSA analysis with AMBER for every frame 
in the trajectory (an equal number of folders will be created).
The electrostatics and vdw are averaged for each frame, and the set of
eigenvectors and values after diagonalization are produced on the average
matrix.
NOTICE: mm_pbsa.pl is very slow, and the execution time grows fast with
the number of amino acids inside the structure. Keep it in mind when
choosing the number of frames to process.
ARGUMENTS FOR MODE REBELOT MULTIFRAME
=========================================================================
Option     example      Type        Description
-------------------------------------------------------------------------

-f(file)   traj.pdb     Input       PDB trajectory to be processed.
-d(domain)              Input(opt)  Alternative energy approximation
                                    via domain decomposition [1].
                                    Default is first eigenvector.
--minrange 5            Input(opt)  Range of residues for which
--maxrange 30                       interaction energy matrix will be
                                    calculated (default: all)
-g                      Input(opt)  Add terms for glycans treatment
-------------------------------------------------------------------------
SYNOPSYS: REBELOT.py  -m m  -f traj.pdb  --minrange 1 --maxrange 100

======================
|   REBELOT CLUSTER  | (-m c)
======================
The cluster-variant of REBELOT is designed to operate on a clustering 
analysis result producing a regular MM-PBSA analysis with AMBER for every 
frame in the trajectory (an equal number of folders will be created).
The only difference with respect to multiframe mode is the averaging of
the energy matrix, which is weighted with the dimension of every cluster.
NOTICE: mm_pbsa.pl is very slow, and the execution time grows fast with
the number of amino acids inside the structure. Keep it in mind when
choosing the number of frames to process.
ARGUMENTS FOR MODE REBELOT CLUSTER
=========================================================================
Option     example      Type        Description
-------------------------------------------------------------------------

-f(file)   traj.pdb     Input       PDB trajectory to be processed.
--cluster  filename     Input       File with the dimension of clusters.
-d(domain)              Input(opt)  Alternative energy approximation
                                    via domain decomposition [1].
                                    Default is first eigenvector.
--minrange 5            Input(opt)  Range of residues for which
--maxrange 30                       interaction energy matrix will be
                                    calculated (default: all)
-g                      Input(opt)  Add terms for glycans treatment
-------------------------------------------------------------------------
SYNOPSYS: REBELOT.py  -m m  -f traj.pdb  --minrange 1 --maxrange 100

=========
| BEPPE | (-m b)
=========
The last mode is for epitope prediction out of a single PDB structure.
The generation of the energy matrix is all in all similar to the other
modes, but highlighting in the final image the areas of energetic
uncoupling instead of stability.
In this case a list of putative epitope patches is returned along with
a pymol script to visualize the results. The script requires a new input
parameter, to define the prediction softness (ratio between sensitivity
and specificity of prediction). This parameter is a number usually
comprised between 5 and 25 (default 15). Use 5/10 for maximum sensitivity
and 20/25 for maximum specificity. Going further is pointless in the vast
majority of cases (although it is possibile).
For the task of epitope prediction, an energy matrix composed only of the
first eigenvector is best suited for small or single domain proteins.
In case of very large multi-domain proteins, or when in presence of
different chains, the domain decomposition energy matrix may be used.
It is possible to export the contact map used for prediction (based on
beta carbon positions and 6 angstroms cutoff) with argument -t.
For more details on the method see [2][3].
Like in single frame mode, it is also possible to use a pre-calculated 
raw matrix and skip the GBSA evaluation with the flag "--matrix". 
Furthermore, it is possible to use a coevolutionary matrix from the 
CoCaInE[4] analysis using the flag "--coevo".
ARGUMENTS FOR MODE BEPPE
=========================================================================
Option      example         Type        Description
-------------------------------------------------------------------------

-f(file)    structure.pdb   Input       PDB structure to be processed.
--matrix    filename        Input(opt)  Raw "ele-vdw.dat" matrix
--coevo     filename        Input(opt)  Coevolutionary potential matrix
-d(domain)                  Input(opt)  Alternative energy approximation
                                        via domain decomposition [1].
                                        Default is first eigenvector.
-c(cutoff)  15              Input(opt)  Prediction softness (default 15)
-t(contact) topology.dat    Output(opt) Exported Contact Map
-M          mlce.dat        Output(opt) Exported MLCE Map
-g                          Input(opt)  Add terms for glycans treatment
-------------------------------------------------------------------------
SYNOPSYS: REBELOT.py  -m b  -f structure.pdb  -c 15 -t


CHANGELOG VERSIONS:
REBELOT ver 1.3 is based on:
ISABEL-M ver 1.0 by Claudio Peri
ISABEL   ver 14.4 by Dario Corrada
pyBEPPE	ver 1.0.1 by Claudio Peri, Riccardo Capelli

[1] Genoni A, Morra G, Colombo G. J Phys Chem B. 2012
[2] Scarabelli G, Morra G, Colombo G. Biophys J. 2010
[3] Peri C, Sole' O.C., Corrada D., [...], Colombo G. Method Mol Biol 2015
[4] Contini A, Tiana G. J Chem Phys 2015
=========================================================================
'''
    print(help)
    sys.exit(0)


def Fine(ret):
    RandomQuote()
    sys.exit(ret)


def RandomQuote():
    quotes = ["<<La uno, la due, olatreee!>> M. Bongiorno",
              "<<Life is Pain>> A. Galli",
              "<<Lo scienziato non puo` pretendere di capire un fenomeno. Puo` solo abituarsi a prendere dimestichezza con esso>> A. Bernacchi",
              "<<They're eating her...and then they're going to eat me.. ... oh my GOOOOOOOOOOOOOO>> Troll 2",
              "<<Lei e` una nanerottola, e la sua pettinatura fa schifo, PORTATELA VIA!>> Star Whores",
              "<<Come dico sempre... non l'hai messo dentro finche` non l'hai messo dentro>> Il Signore dell'Agrizia",
              "<<Bazinga!>> Sheldon Cooper",
              "knock-knock-knock <<Penny!>> knock-knock-knock <<Penny!>> knock-knock-knock <<Penny!>> Sheldon Cooper",
              "<<Tutto quello che ho e` questo pollo di gomma con una carrucola in mezzo>> Guybrush Treepwood",
              "<<We do what we must, because we can. For the good of all of us, except the ones who are dead>> Portal credit song",
              "<<All science is either physics or stamp collecting>> Ernest Rutherford",
              "<<I'm going to be a boooss>> <<Am I hallucinating?>> <<You like meeeen!>> <<You can see into my mind?!?>> <<No>> <<Fuck>> Metal Gear Awesome.",
              "<<I puffi sanno che un tesoro c'e`>> C. D'Avena",
              "<<Cowabunga? Cowafuckinpieceo'dogshit!>> Angry VideoGame Nerd",
              "<<Tutto disintegra quando gli girano le lame boomerang>> Daltanius opening",
              "<<Hai mai fatto un sogno tanto realistico da sembrarti vero? E se da un sogno cosi` non ti dovessi piu` svegliare... come potresti distinguere il mondo dei sogni da quello della realta`?>> The Matrix",
              "<<E' inevitabile, signor Anderson>> The Matrix",
              "<<Mi consenta>> S. Berlusconi",
              "<<Werewolf!>> <<Werewolf?>> <<There>> <<What?>> <<There, wolf. There, castle>>",
              "<<I personally believe, that U.S. Americans, are unable to do so, because uh, some, people out there, in our nation don't have maps and uh... I believe that our education like such as in South Africa, and the Iraq, everywhere like such as... and, I believe they should uh, our education over here, in the U.S. should help the U.S. or should help South Africa, and should help the Iraq and Asian countries so we will be able to build up our future, for us.>> Miss Teen South Carolina 2007",
              "<<Io... vado a pisciare>> <<...erano tre anni che non andava a pisciare>> Evanghelion LSD",
              "<<Quando il gioco si fa duro... io vorrei essere da un'altra parte>> PK",
              "<<The cake is a lie>>",
              "<<Thank you, Mario, but our princess is in another castle!>> Toad",
              "<<How are you gentlemen! All your base are belong to us>> Cats",
              "<<Venendo qui, hai per caso visto un cartello davanti a casa mia con scritto ''deposito di negri morti''?>> Pulp Fiction",
              "<<Non ci vuole un pennello grande, ma un grande pennello>>",
              "<<Aaaah, ma e` Ronco>>",
              ">:3   <<Jesus Christ it's a lion get in the car!>>",
              "<<It's over nine thousaaaaaaaaaaaaaand!!>> Vegeta ",
              "<<Put a banana in your ear>> Charlie the Unicorn",
              "<<Balls to you!>> Liza Minnelli",
              "<<Sono bat-terie, e io sono Bat-man, il paladino della giustizia>> <<Con quelle tette? Io pensavo fossi Bat-tona, la paladina del ribaltabile!>> Evanghelion LSD",
              "<<Solo Puffin ti dara` forza e grinta a volonta`!>> Chi trova un amico trova un tesoro",
              "<<Dacci il mezzuomo, figlia degli Aerosmith>> <<Venite a prenderlo, figli dei Mordorhead!>> Il Signore dell'Agrizia",
              "<<Maremma impestata! Millenni di esperienza per farci fregare ancora col trucco dei cuscini!>> Il Signore dell'Agrizia",
              "<<Contro coloro che vi inseguono servira` ben piu` dei miei occhi verdi penetranti>> <<ma non sono verdi...>> Il Signore dell'Agrizia",
              "<<Urca Urca tirulero, oggi splende il sol>> Cantagallo",
              "<<Si - puo` - fare!>> Frankenstein Junior",
              "<<Per questo e` cosi` illogico essere uno dei Puffi, perche`... Che cazzo vivi a fare, se non hai il pisello?>> Donnie Darko",
              "<<Per me, nnumero uno!>> Dan Peterson",
              "<<Come ti chiami?>> <<Luna!>> <<Ma che bel nome! E perche` non Le due, o Le tre!>> Laime Nailo",
              "<<...>>  Wile E. Coyote",
              "<<Don't panic>> Hitchickers guide to the galaxy",
              "<<Mi hai fatto ingoiare la gomma... Mi rimarra` nell'intestino per sette anni!>> Scott Pilgrim vs The World",
              "<<I'm in lesbians with you>> Scott Pilgrim vs The World",
              "<<Developers Developers Developers Developers>> Steve Ballmer",
              "<<Sei tu un DIO?!?>> Gozer il gozeriano",
              "<<Procol harum zovirax platinum collection sandeibladisa`ndei!>> <<eh?!?>> <<Puppa!!>> Don Zauker",
              "BUDDA BUDDA BUDDA <<Muori, anima dannata di Farinata degli Uberti! Muori! Muori!!!>> <<Oorghlll--->> Dante II",
              "<<badger badger badger badger badger badger badger badger badger badger badger badger mushroom mushroooom>>",
              "<<UN POLLOOOOOOOOOOOO!!!>> Richard Benson",
              "<<Una scrivania ordinata e` sintomo di una mente malata>> O. Wilde",
              "<<My body is ready>> Reginald Fils Aime",
              "<<This is the best handgun ever made: the Colt Single Action Army. Six bullets... more than enough to kill anything that moves>> Revolver Ocelot",
              "<<Vendero` care le penne! ...anche se non so quale sia la loro quotazione sul mercato... >> PK",
              "<<Winter is coming>> House Stark motto",
              "<<Oppa gangnam style>>",
              "<<FUS RO DAH>>",
              "<<Ser Gregor, Dunsen, Polliver, Chiswyck, Weese , Raff Dolcecuore, Messer Sottile , il Mastino. Ser Amory, ser Ilyn, ser Meryn, re Joffrey, regina Cersei>>",
              "<<Feeeeenomenali poteri cosmici... in un minuscolo spazio vitale>>",
              "<<Siamo in missione per conto di Dio>> Elwood Blues",
              "<<Fletto i muscoli e sono nel vuoto>> Ratman",
              "<<Deep into that darkness peering, long I stood there, wondering, fearing, doubting, dreaming dreams no mortal ever dared to dream before>> E.A. Poe",
              "<<When I get sad, I stop being sad and be awesome instead. True Story>> Barney Stinson",
              "<<Ci scusiamo per il disagio>> Trenitalia",
              "<<E' vietato attraversare i binari. Servirsi del sottopassaggio>>",
              "<<Could you please volunteer?>>",
              "<<Che bello, due amici, una chitarra e uno spinello>> Stefano Rosso",
              "<<Fettine panate! Dai amore, ti prego, fai le fettine panate, tu le fai da impazzire di bonta`>> Mirko",
              "<<Prima ero in collera con voi, Locksley, ma adesso sono addirittura incazzato!>> Smerdino, Sceriffo di Ruttingham",
              "<<Ora basta! Il re, delle foreste illegali maialino uccisamente sapete selvatici in un a del?>> Smerdino, Sceriffo di Ruttingham",
              "<<Ecco le chiavi di tutte le porte di tutti i forzieri. Ecco le chiavi di tutte di tal, ecco le chiavi di tutti di tutte. Tutte di tutte di tutte di tu. Tutte di tutte di tutti di tal, tutte di tutti di tutte di chiavi, tutte di tutti di tutte di VADO.>> Ulbabrab",
              "<<Vieni giu` subito che ti voglio ammazzare adesso!>> Ulbabrab",
              "<<Eammecheccacchiomenefregamme`, io c'ho il diesel!>> Italiano Medio",
              "<<Vedo la gente scema>> Il Sesto Scemo",
              "<<beam me up, Scotty>>",
              "<<Guarda dietro di te! Una scimmia a tre teste!>> Guybrush Threepwood",
              "<<Svegliate Rei>> <<Ma non e` ancora morta?>> <<Pero` e` viva>> <<... non ho capito>> Evanghelion LSD",
              "<<The Door>> <<The door?>> <<The Door!>> <<What is the door?>> <<The Door is everything! All that once was and all that will be! The Door controls Time and Space! Love and Death! The Door can see into your mind! The Door can see into your SOUL!>> <<Really, th-the door can do all that?>> <<eh, no>> Charlie the Unicorn.",
              "<<suca>>",
              "<<42>>",
              "<<The true object of all human life is play. Earth is a task garden; heaven is a playground>> G.K. Chesterton",
              "<<The poets have been mysteriously silent on the subject of cheese>> G.K. Chesterton",
              "<<The object of opening the mind, as of opening the mouth, is to shut it again on something solid>> G.K. Chesterton",
              "<<Usa la forza, Harry>> J.R.R. Kirk",
              "<<The fewer the facts, the stronger the opinion>> R.A. Heinlein",
              "<<Research is what I'm doing when I don't know what I'm doing>> Wernher von Braun",
              "<<The saddest aspect of life right now is that science gathers knowledge faster than society gathers wisdom>> Isaac Asimov",
              "<<If an elderly but distinguished scientist says that something is possible, he is almost certainly right; but if he says that it is impossible, he is very probably wrong>> Arthur C. Clarke",
              "<<Science may never come up with a better office communication system than the coffee break>> Earl Wilson",
              "<<Science is a differential equation. Religion is a boundary condition>> Alan Turing",
              "<<When I die, I'm leaving my body to science fiction>> Steven Wright"]

    print(random.choice(quotes))


if __name__ == '__main__':
    main()
