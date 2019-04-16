#!/usr/bin/python
import sys
import os.path

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

helixCount = 0
sheetCount = 0
atomCount  = 0
lineNumber = 1
line = inFile.readline()
for line in inFile:
	if "HELIX" in line:
		lineNumber += 1
		helixCount += 1
		print(lineNumber, ": ", line, sep = '', end = "")
		print(line, sep = '', end = "", file = outFile)
	elif "SHEET" in line:
		lineNumber += 1
		sheetCount += 1
		print(lineNumber, ": ", line, sep = '', end = "")
		print(line, sep = '', end = "", file = outFile)
	elif "ATOM" in line:
		lineNumber += 1
		atomCount += 1
		print(lineNumber, ": ", line, sep = '', end = "")
		print(line, sep = '', end = "", file = outFile)
	else:
		lineNumber += 1

structureCount = helixCount + sheetCount
print("")
print("There are", structureCount, "protein secondary structures described in the parsed PDB file!")
print("Total Alpha Helices:", helixCount)
print("Total Beta Sheets:", sheetCount)
print("Total Atoms:", atomCount)
inFile.close()
outFile.close()