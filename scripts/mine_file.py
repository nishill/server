# mine_file.py

"""
Takes as input a column and filename. The file must include columns (tsv or csv).
The column name is processed into a json file with the columns as keys in the .json.

run:
python mine_file.py -c [COLUMN_NAME] -o [OUTPUT_FILE] [filename] 


"""

import os
import sys
import argparse 
import json 
import csv

# parses file to into json document specified by the column name
def parse_file (filename, column, out):
	vacant_dict = {}
	file_name,file_extension = os.path.splitext(filename)
	file_type = ''
	if ( file_extension == '.tsv' ): file_type = "\t"
	elif (file_extension == '.csv'): file_type = ","

	print("Loading biodata from {}".format(filename))
	empty_dict = {}
	with open(filename, 'rw') as tsvfile:
		reader = csv.DictReader(tsvfile,delimiter=str(file_type), quoting=csv.QUOTE_NONE)
		for row in reader:
			type = row[column]  
			vacant_dict[type] = empty_dict
	data = json.dumps(vacant_dict)
	with open(out, 'w') as file:
		file.write(data)

# main()
def main():
	parser = argparse.ArgumentParser(prog='mine_file',description=
							"Processes a specific column in a file into a .json " 
							"file with that column's cells as its keys.")

	parser.add_argument('filename', metavar='FILE_NAME',
                    help='specifies the name of file to be processed')
	parser.add_argument('-c', metavar='COLUMN_NAME', type=str,
					help='specifies the column to be processed into the json file')
	parser.add_argument('-o', metavar= 'OUTPUT_FILE', type=str, default='data.json',
					help= 'specifies the name of the output file; it will be '
						  'defaulted to data.json if none is specified')

	args = parser.parse_args()
	column = args.c
	filename = args.filename
	output = args.o
	parse_file(filename, column, output)		
	print ("Done! Output written to {}".format(output))

# run main ----------------->
if __name__ == "__main__":
	main()
