import os, pickle, random, re, sys, time
import pandas as pd
import numpy as np
sys.path.append("/home/heymann/Desktop/deeplearn2/inDelphi-model") #<-- https://github.com/maxwshen/inDelphi-model
import sklearn
import inDelphi
from utils import collect_ground_truth
from metrics import accuracy_type_agreement, accuracy_type_mutation, accuracy_type_disrupt_reading_frame, stats_frameshift
from metrics import compile_mutations, stats_frequency_mutations


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

def normalize(input, decimals=5):
	return np.round((input/np.sum(np.array(input))).tolist(),decimals)

def format_ground_truth(gt_target, oligo_id):
	# A wrapper function that formats the output of 
	# collect_ground_truth into a proper dataframe.
	#print(gt_target)
	actual = pd.DataFrame.from_records(collect_ground_truth(gt_target, oligo_id))
	actual.columns = ['type', 'n_actual', 'seq']; 
	actual = actual[['type', 'n_actual']]
	actual['type'] = actual['type'].str.extract(r'^(\w\d?\d)')
	actual['actual'] = normalize(list(actual['n_actual']))
	return actual


def compile_mutations():
	# A function that creates a dictionary of all the
	# needed mutations in indel-predictions
	mut_types = dict()
	for d in reversed(range(1,200)):
		mut_types["D{0}".format(d)] = float(0)
	for i in range(1,100):
		mut_types["I{0}".format(i)] = float(0)
	return mut_types

def format_predction(seq, pam_idx):
	# A wrapper function that formats the output of
	# the inDelphi prediction method
	cutsite = pam_idx - 3
	seqA = seq[0:cutsite]
	seqB = seq[cutsite::]
	pred_df, stats = inDelphi.predict(seq, cutsite)
	pred_df = pred_df.groupby(['Category','Length'],as_index=False).agg({'Predicted frequency': 'sum'})
	pred_df['New category'] = pred_df['Category'].astype(str).str[0] + pred_df['Length'].map(str)
	pred_df = pred_df[['New category', 'Predicted frequency']]
	pred_df.columns = ['type','pred']
	pred_df['type'] = pred_df['type'].apply(lambda x: x.upper())
	pred_df['pred'] = pred_df['pred'].astype(float)/100 #.round(decimals = 4)/100

	return pred_df, stats


def runDelphi(dataframe, truth_reference, max_oligos = 1, cell_type = 'mESC', file_prefix = "inDelphi", pathout = ""):
	# init model
	inDelphi.init_model(celltype = cell_type)

	print("init runDelphi.."); time.sleep(1)

	# Open files so that data can be written
	outfile1 = open(pathout + file_prefix + "_" + cell_type + "_statistics.txt", 'w')
	outfile2 = open(pathout + file_prefix + "_" + cell_type + "_indels_frequency_predicted.txt", 'w')
	outfile3 = open(pathout + file_prefix + "_" + cell_type + "_indels_frequency_actual.txt", 'w')
	outfile4 = open(pathout + file_prefix + "_" + cell_type + "_modelstats.txt", 'w')

	# prepare dictionaries for storing relevant data
	mutations_reference_actual = compile_mutations()
	mutations_reference_pred = compile_mutations()

	# Save some info to the outfiles
	outfile2.write("oligo_id" + "\t" + "\t".join([v for v in mutations_reference_pred.keys()]) + "\n")
	outfile3.write("oligo_id" + "\t" + "\t".join([v for v in mutations_reference_pred.keys()]) + "\n")
	errs = list()
	klavg = 0
	# Variables for loop
	elapsed_time = 0
	iteration = 0; n_oligos = 0
	while n_oligos < max_oligos:	
		start = time.time()
		oligo = dataframe.iloc[[iteration]].values.tolist()[0]
		oligo_id = oligo[0]; seq = oligo[2]; pam_idx = int(oligo[6])
		sys.stdout.write("Query [{0}/{1}]: {2}\n".format(n_oligos, max_oligos, oligo_id))
		try: 
			# Get a dataframe from the ground truth and
			# the the prediction of in-delphi. Merge afterwards.
			#print("pam=",pam_idx)
			actual = format_ground_truth(truth_reference, oligo_id)
			prediction, stats = format_predction(seq, pam_idx) # cutsite is modified in func
			df = actual.merge(prediction,how='outer', on='type' ).fillna(0)

			## Run statistical analysis
			loss = L(list(df['pred']),list(df['actual']))
			stat01 = str(1-stats['Frameshift frequency']/100) #stats_frameshift(df, 'pred')
			stat02 = stats_frameshift(df, 'actual')
			acc01 = accuracy_type_agreement(x,y, top = 1)
			acc02 = accuracy_type_agreement(x, y, top = 2)
			acc03 = accuracy_type_agreement(x, y, top = 3)
			acc04 = None #accuracy_type_disrupt_reading_frame(df)
			acc05 = None #accuracy_type_disrupt_reading_frame(df)
			acc06 = None #accuracy_type_disrupt_reading_frame(df)
			stat_freq_pred = stats_frequency_mutations(df, 'pred', mutations_reference_pred.copy())
			stat_freq_actual = stats_frequency_mutations(actual, 'actual', mutations_reference_actual.copy())

			#time.sleep(5)
			outfile2.write(oligo_id + "\t" + "\t".join([str(v) for v in x]) +"\n") #predicted
			outfile3.write(oligo_id + "\t" + "\t".join([str(v) for v in y]) +"\n")
			outfile1.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\t{9}\n".format(oligo_id, loss, acc01, acc02, acc03, acc04, acc05, acc06, stat01, stat02))
			outfile4.write("\t".join([str(v) for k, v in stats.items()]))


			# Gather statistics about the code, i.e. timing
			klavg += loss
			n_oligos += 1
			
			end = time.time(); elapsed_time += (end-start)

			if n_oligos % 20 == 0:
				print("KL=",np.round(klavg/n_oligos,5))
				time.sleep(1)

			sys.stdout.write("SElapsed Time: " + str(round(elapsed_time,0)) + "s\n")
		except Exception as e:
			sys.stdout.write("\n>>>ERROR: " + str(e) + "\n")
			errs.append(oligo_id)
			print("CUR ERRORS: " + str(errs))
			time.sleep(3)
			pass
		iteration += 1
	# close all files
	outfile1.close(); outfile2.close()
	outfile3.close(); outfile4.close()
