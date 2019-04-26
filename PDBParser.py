#!/usr/bin/python
import sys
import os.path
#from collections import OrderedDict

scriptpath = os.path.dirname(__file__)
infile_name = input("Enter in a PDB filename with extension (.pdb): ")
outfile_name = input("Enter in an output file name with extension (.pdb): ")

try:
	inFile = open(infile_name, "r")
except IOError:
	print("There was an error with reading from the input file.")
	sys.exit()

try:
	outFile = open(outfile_name, "a")
except IOError:
	print("There was an error with writing to the output file.")
	sys.exit()

print("Input File:", infile_name)
print("Output File:", outfile_name)

"""
Read all lines containing HELIX or SHEET, and then extract residue sequence number data from each line
to determine associated atoms which comprise the secondary structures. Print to an output PDB file only
those lines containing HELIX or SHEET, and the atoms which share a residue sequence number within one of
the specified ranges.

As each line is read in and stored as a string variable, this means that you should begin indexing (counting columns) from 0.
HELIX residue sequence numbers can be found in columns 21-24 (initial) - entry 6 on line, and columns 33-36 (terminal), entry 9 on line.
SHEET residue sequence numbers can be found in columns 22-25 (initial) - entry 7 on line, and columns 33-36 (terminal), entry 10 on line.
ATOM singular residue sequence number can be found in columns 22-25 - entry 6 on line.
"""

helixCount = 0
sheetCount = 0
atomCount  = 0
lineNumber = 1

# List atoms contains the residue sequence number ranges of all atoms which comprise an alpha helix or beta sheet
atoms = []
# List atom contains the individual sequence number of every initial or terminal that has no partner not already found in atoms
# atom = []
line = inFile.readline()

# Store all lines containing string "ATOM" into a list for use in search and comparison of atom residue sequence number to alpha helix/beta sheet [initial, terminal] residue sequence range
totalAtoms = []
for line in inFile:
	if "ATOM" in line and "REMARK" not in line:
		totalAtoms.append(line)
		#print(line, end="")

# Before search for and printing of atoms within protein secondary structure residue sequence range, search inFile for all instances of "HELIX" and "SHEET" and print those lines first to standard output and the outFile
inFile.seek(0)
for line in inFile:
	if "HELIX" in line or "SHEET" in line:
		print(line, sep = '', end = "")
		print(line, sep = '', end = "", file = outFile)

inFile.seek(0)
for line in inFile:
	if "HELIX" in line:
		lineNumber += 1
		helixCount += 1
		
		# Access characters from string line that refer to initial and terminal residue sequence numbers
		# Update atoms to have sequence numbers listed in order of initial then terminal after conversion to integer types
		
		# print(line[21:25], line[33:37])
		initial = int(line[21:25])
		terminal = int(line[33:37])
		# print(int(initial), int(terminal))

		# print(line, sep = '', end = "")
		# print(line, sep = '', end = "", file = outFile)

		# Print atoms to the output file which have a residue number within the range [initial, terminal]
		for atom in totalAtoms:
			residueSequenceNumber = int(atom[22:26])
			if residueSequenceNumber >= initial and residueSequenceNumber <= terminal and atom not in atoms:
				# print(residueSequenceNumber, initial, terminal)
				atomCount += 1
				print(atom, sep = '', end = "")
				print(atom, sep = '', end = "", file = outFile)
				atoms.append(atom)
	elif "SHEET" in line:
		lineNumber += 1
		sheetCount += 1

		# Access character from string line that refer to initial and terminal residue ssequence numbers
		# Update atoms to have sequence numbers listed in order of initial then terminal after conversion to integer type

		# print(line[22:26], line[33:37])
		initial = int(line[22:26])
		terminal = int(line[33:37])
		# print(int(initial), int(terminal))

		# print(line, sep = '', end = "")
		# print(line, sep = '', end = "", file = outFile)

		# Print atoms to the output file which have a residue number within the range [initial, terminal]
		# Print atoms to the output file which have a residue number within the range [initial, terminal]
		for atom in totalAtoms:
			residueSequenceNumber = int(atom[22:26])
			if residueSequenceNumber >= initial and residueSequenceNumber <= terminal and atom not in atoms:
				# print(residueSequenceNumber, initial, terminal)
				atomCount += 1
				print(atom, sep = '', end = "")
				print(atom, sep = '', end = "", file = outFile)
				atoms.append(atom)
	else:
		lineNumber += 1

structureCount = helixCount + sheetCount
print("")
# print("Base Atoms List:", atoms)
# print("Total Atoms: ", len(atoms))
# print("Individual Atom List:", atom)

"""
Does not produce a list containing n*2 elements, making it impossible to dynamically assign initial <-> terminal sequence number pairs
Will remove duplicate terminal residue sequence numbers for having already occurred, but leave the initial value for being unique, and vice-versa
This leaves incorrect pairs when reading pairss of values from the list in sequential order. (Example: ..., 2, 99, ...)

pruned = OrderedDict((x, True) for x in atoms).keys()
print("Pruned Atoms List:", pruned)
"""

print("There are", structureCount, "protein secondary structures described in the parsed PDB file!")
print("Total Alpha Helices:", helixCount)
print("Total Beta Sheets:", sheetCount)
print("Total Atoms:", atomCount)
inFile.close()
outFile.close()