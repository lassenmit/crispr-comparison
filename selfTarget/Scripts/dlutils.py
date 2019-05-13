import os, pickle, random, re, sys
import pandas as pd
import numpy as np
from dlpredict import predict_single
from dlmetrics import accuracy_type_agreement, accuracy_type_mutation, accuracy_type_disrupt_reading_frame, stats_frameshift
from dlmetrics import compile_mutations, stats_frequency_mutations

def sample_oligos(file, n = 100, seed = 0):
	# sample oligos randomly
	random.seed(seed)
	data = pd.read_csv(file, sep ="\t", header =(0)) 
	samples = random.sample(range(data.shape[0]), n)
	return data.iloc[samples]

def normalize(input, decimals=5):
	return np.round((input/np.sum(np.array(input))).tolist(),decimals)

def KL(P,Q):
	### Epsilon is used here to avoid conditional code for
	### checking that neither P nor Q is equal to 0. """
	epsilon = 0.00001
	P = np.asarray(P, dtype=np.float)
	Q = np.asarray(Q, dtype=np.float)
	P = P+epsilon
	Q = Q+epsilon

	divergence = np.sum(P*np.log(P/Q))
	return divergence

def available_oligos(file, gtdir, wfile = "sampled_oligos.txt", n = 100, seed = 42):
	# A function that returns the oligos from 'file' and checks whether they
	# can be found in the gtdir. The 'n' oligos that are found are written to
	# wfile. Uses random subsampling to look for oligos. 
	oligos = sample_oligos(file, n=n, seed = seed)
	oligo_count = 0; iteration = 0
	with open(wfile, 'w') as outfile:
		while oligo_count < n:

			oligo = oligos.iloc[int(iteration)].values.tolist()[0:8]
			oligo_id = oligo[0]; direction = oligo[7]
			try:
				collect_ground_truth(gtdir, oligo_id)
			except Exception as e:
				continue
			print(direction)
			if direction == "FORWARD":
				write_oligo = '\t'.join(map(str,oligo)) + "\n"
				outfile.write(write_oligo)
				oligo_count += 1

def predict_and_collect_ground_truth(df,  gt_target, n=1000, wfile = 'losses_top_3.txt'):
	
	#' @description opens up a pandas dataframe of oligos.
	#' Takes the 2nd column as oligo ID and makes a prediction
	#' using the selfTarget framework
	#' @param df a pandas dataframe with the columns in the format
	#' of <ID> <GUIDE> <SEQUENCE> <?> <?> <?> <PAM ID>
    #' @param n how many lossess should be gathered? 

	mutations_reference_actual = compile_mutations()
	mutations_reference_pred = compile_mutations()
	outfile = open(wfile, 'w')
	outfile2 = open("{0}.freq_pred".format(wfile), 'w')
	outfile3 = open("{0}.freq_actual".format(wfile), 'w')

	outfile2.write("\t".join([v for v in mutations_reference_pred.keys()]) + "\n")
	outfile3.write("\t".join([v for v in mutations_reference_pred.keys()]) + "\n")
	#outfile.write("oligo_id KL acc01 acc02 acc03 acc04 acc05 acc06 stat01 stat02")
	oligo_count = 0; iteration = -1
	acc_count = list()


	

	print(df.head)


	for oligo_line in df.values.tolist():


		#oligo_line = df.iloc[int(iteration)].values.tolist()
		oligo_id = oligo_line[0]
		seq = oligo_line[2]
		pam_idx = int(oligo_line[6])


		# This is ugly, but terminates the run
		# when we have enough data.
		if oligo_count >= n:
			break

		# get ID and Oligo

		print("\n>{0} OLIGO_ID='{1}'".format(oligo_count, oligo_id))

		# Prectict using selfTarget
		try:

			# Get the actual mutations by reading from the ground truth
			# Furthermore, discard the sequence, as this is collected in type and n_actual
			actual = pd.DataFrame.from_records(collect_ground_truth(gt_target, oligo_id))
			actual.columns = ['type', 'n_actual', 'seq']
			actual = actual[['type', 'n_actual']]

			# Use the inputted platform to predict profiles. for SelfTarget the output 
			# is a .txt file that will have to be read by pandas
			print("2. PREDICTING OLIGO..")
			predict_single(target_seq = seq, pam_idx = pam_idx, output_prefix = "single_data")
			prediction = pd.read_csv("single_data_predictedindelsummary.txt", sep ="\t", header =(0)) 
			prediction.columns = ['type', 'rmme1', 'n_pred']
			prediction = prediction[['type', 'n_pred']]


			# Merge the predicted and actual data.frames through an outer-join.
			# This ensures that all mutations types are kept
			df = prediction.merge(actual,how='outer', on='type' )
			df = df.fillna(0) # fill nas to zero, as they are just unseen mutations
			df['pred'] = normalize(list(df['n_pred'])); df['actual'] = normalize(list(df['n_actual']))

			print(df)

			# Summarize the results and write the to files, so they can be
			# anaylyzed later.
			print("4. ACC AND LOSS..")
			loss = KL(list(df['pred']),list(df['actual']))
			stat01 = stats_frameshift(df, 'pred')
			stat02 = stats_frameshift(df, 'actual')
			acc01 = accuracy_type_agreement(list(df['n_pred']),list(df['n_actual']), top = 1)
			acc02 = accuracy_type_agreement(list(df['n_pred']),list(df['n_actual']), top = 2)
			acc03 = accuracy_type_agreement(list(df['n_pred']),list(df['n_actual']), top = 3)
			acc04 = accuracy_type_mutation(df)
			acc05 = accuracy_type_mutation(df, loose = True)
			acc06 = accuracy_type_disrupt_reading_frame(df)


			stat_freq_pred = stats_frequency_mutations(df, 'n_pred', mutations_reference_pred.copy())
			stat_freq_actual = stats_frequency_mutations(df, 'n_actual', mutations_reference_actual.copy())
			outfile2.write("\t".join([str(v) for v in stat_freq_pred.values()]) +"\n")
			outfile3.write("\t".join([str(v) for v in stat_freq_actual.values()]) +"\n")
			print(stat_freq_pred)
			print(stat_freq_actual)

			print("deletions, actual=", stat02, "pred=", stat01)

			acc_count.append(acc06)


			outfile.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\t{9}\n".format(oligo_id, loss, acc01, acc02, acc03, acc04, acc05, acc06, stat01, stat02))
			oligo_count += 1
			
		except Exception as e:
			print("ERROR: " + str(e))
			continue

		iteration += 1
		print("Current running accuracy:" + str(round(100*sum(acc_count)/len(acc_count),3)) +" %")

	outfile.close()
	outfile2.close()
	outfile3.close()

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
		raise TypeError(str(oligo_dir) + " is not a valid path!")