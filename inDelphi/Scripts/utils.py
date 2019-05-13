import os, pickle, random, re, sys
import pandas as pd
import numpy as np

def get_oligo_dir(oligo, maxlen = 5):
	# a function that returns the correct directory
	# extension of an oligonucleotide. This is accomplished
	# by regexing the the digits and extending it with zeros.
	oligo = list(re.findall(r'\d+', oligo)[0])
	oligo = ((maxlen-len(oligo))*['0'] + oligo)
	if oligo[0] == "0":
		return oligo[1]
	else:
		return ''.join(oligo[0:2])

def id_interval(file, oligo):
	# a function that returns the numeric interval
	# of the directory files.
	return [int(x) for x in re.findall(r'\d+', file)]
	
def collect_ground_truth(directory, oligo):
	
	# data = collect_ground_truth(gt_target, "Oligo80492")
	# df =  pd.DataFrame.from_records(data)
	# print(df)

	# find the oligos and oligo directory
	oligo_id = int(re.findall(r'\d+', oligo)[0])
	oligo_dir = directory + "Oligos_" + get_oligo_dir(oligo)
	print(oligo_dir)
	if os.path.exists(oligo_dir):
		oligo_found = False
		# go through the directories to open correct file
		for dirpath,_,filenames in os.walk(oligo_dir):
			for f in filenames:
	
				# enumeate the interval in the oligo id in order
				# to find the corresponding file, wihtout searching all
				interval =  id_interval(f, oligo)
				if oligo_id >= interval[0] and oligo_id <= interval[1]:
					fpath = os.path.abspath(os.path.join(dirpath, f))
					# Ready statefull parsing
					flag = False
					lines = []
					with open(fpath, 'r') as infile:
						for line in infile:

							# if the correct oligo starter has been
							# found, set the flag to TRUE
							if line.startswith("@@@"):
								oligo_line = "".join(list(line)[3:-1])
								if oligo_line == oligo:
									flag = True
									oligo_found = True
								else:
									flag = False

							# Collect data if flag is Truue
							if flag == True:
								# only add to file if length is greater
								# than 1 and that the oligo id is not there
								rline = line.strip("\n").split("\t")
								if len(rline) > 1:
									rline[1] = int(rline[1])
									lines.append(rline)

		# Check if the oligo was actually
		# present in the file
		if oligo_found == True:
			return lines
		else:
			raise TypeError(str(oligo) + " was not found in file!")
	else:
		raise TypeError("'" + str(oligo_dir) + "' is not a valid path!")

