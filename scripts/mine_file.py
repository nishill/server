# mine_file.py

"""
Takes as input a column and filename. The file must include columns (tsv or csv).
The column name is processed into a json file with the columns as keys in the .json.

run:
python mine_file.py -c [COLUMN_NAME] [filename] 

"""

import os
import sys
import argparse 
import json 
import csv

# parses file to into json document specified by the column name
def parse_file (filename, column):
	vacant_dict = {}
	print("Loading biodata from {}".format(filename))
	empty_dict = {}
	with open(filename, 'rw') as tsvfile:
		reader = csv.DictReader(tsvfile,delimiter=str("\t"), quoting=csv.QUOTE_NONE)
		for row in reader:
			type = row[column]  
			vacant_dict[type] = empty_dict
	data = json.dumps(vacant_dict)
	out = 'data.json'
	with open(out, 'w') as file:
		file.write(data)
	return out

# main()
def main():
	parser = argparse.ArgumentParser(prog='mine_file',description=
							"Processes a specific column in a file into a .json " 
							"file with that column's cells as its keys.")

	parser.add_argument('filename', metavar='FILE_NAME',
                    help='specifies the name of file to be processed')
	parser.add_argument('-c', metavar='COLUMN_NAME', type=str,
					help='species the column to be processed into the json file')
	args = parser.parse_args()
	column = args.c
	filename = args.filename
	out = parse_file(filename, column)		
	print ("Done! Output written to {}".format(out))

# run main ----------------->
if __name__ == "__main__":
	main()
