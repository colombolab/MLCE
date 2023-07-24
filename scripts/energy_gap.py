import os
import sys
import numpy as np

###########################################################
# usage: python3 energy_gap.py <eigval_file> <output_file>
###########################################################


# Read in the eigenvalues
eigval_file = sys.argv[1]

# Output
output_file = sys.argv[2]

# Make a list of the eigenvalues
eigval_list = []
with open(eigval_file, 'r') as f:
    for line in f:
        eigval_list.append(float(line.strip()))

# Calculate the energy gap between the first in position and the second in position (absolute value)
energy_gap = abs(eigval_list[0] - eigval_list[1])

# Calculate the average value of each energy gap between consecutive eigenvalues
energy_gap_avg = np.mean([abs(eigval_list[i] - eigval_list[i+1]) for i in range(0, len(eigval_list)-1)])

# Result
result = energy_gap / energy_gap_avg

# Write the result to the output file
with open(output_file, 'w') as f:
    f.write('Energy gap: {}\n'.format(energy_gap))
    f.write('Average energy gap: {}\n'.format(energy_gap_avg))
    f.write('Result: {}\n'.format(str(result)))
