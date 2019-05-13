import numpy as np

def accuracy_type_top(x,y, top=3):
	# compares the top three mutations in the predicted and
	# the actual profiles. Returns a discrete value depending
	# on how many of X is actually present in Y. 
	x = np.asarray(x, dtype=np.int) #.argsort()[-top:][::-1]
	y = np.asarray(y, dtype=np.int) #.argsort()[-top:][::-1]

	xind = sorted(range(len(x)), key=lambda i: x[i])[-top:]
	yind = sorted(range(len(y)), key=lambda i: y[i])[-top:]

	rsum = 0
	for i in range(top):
		if xind[i] in yind:
			rsum += 1
	return rsum/top

def accuracy_type_agreement(x,y, top = 3):
	# compares the top three mutations in the predicted and
	# the actual profiles. If just a single mutation in 'top' of X
	# is present in 'top' of Y, then returns 1.
	x = np.asarray(x, dtype=np.int) #.argsort()[-top:][::-1]
	y = np.asarray(y, dtype=np.int) #.argsort()[-top:][::-1]
	xind = sorted(range(len(x)), key=lambda i: x[i])[-top:]
	yind = sorted(range(len(y)), key=lambda i: y[i])[-top:]

	rsum = 0
	for i in range(top):
		if xind[i] in yind:
			rsum = 1
	return rsum


def argmax_mutation(df, col):
	# A funnction that returns the
	# the mutation type of the best col
	best = df.iloc[[df[col].idxmax()]]
	mut = best['type'].iloc[0]
	return mut[0:mut.find('_')]

def check_reading_frame(string):
	count = int(string[1:])
	return True if count % 3 == 0 else False

def accuracy_type_disrupt_reading_frame(df):

	pred = check_reading_frame(argmax_mutation(df, 'n_pred'))
	actual = check_reading_frame(argmax_mutation(df, 'n_actual'))
	if pred == actual:
		print("Reading frame correctly predicted")
		return 1
	else:
		print("reading frame inccorectly predicted")
		return 0

def accuracy_type_mutation(df, loose = False):

	# Find best rows
	best_pred = df.iloc[[df['n_pred'].idxmax()]]
	best_actual = df.iloc[[df['n_actual'].idxmax()]]

	# Localize type of mutation
	mut_pred = best_pred['type'].iloc[0]
	mut_actual = best_actual['type'].iloc[0]
	mut_pred = mut_pred[0:mut_pred.find('_')]
	mut_actual = mut_actual[0:mut_actual.find('_')]
	
	# loose criteria, i.e if both are insertions
	# or both are deletions, then it will be acceepted
	if loose == True:
		mut_actual = mut_actual[0]
		mut_pred = mut_pred[0]

	print("PRED={0}, ACTU={1}".format(mut_pred, mut_actual))
	if mut_actual == mut_pred:
		return 1
	else:
		return 0


def stats_frameshift(df, col = 'pred'):
	# find the overall percentages of frameshifts by counting
	# amount of deletions and insertions.
	inframe_mutations = 0
	for index, row in df.iterrows():
		mut = row['type']
		mut = mut[0:mut.find('_')]
		mut_type = mut[0]; mut_delta = int(mut[1:])
		if mut_delta % 3 == 0:
			inframe_mutations += row[col]
	return inframe_mutations 



def stats_compare_delins(df, col = 'pred'):
	# find the overall percentages of frameshifts by counting
	# amount of deletions and insertions.
	pins = 0; pdel = 0 
	for index, row in df.iterrows():
		mut = row['type']
		mut = mut[0:mut.find('_')]
		mut_type = mut[0]; mut_delta = int(mut[1:])
		if mut_type == "D":
			pdel += row[col]
		elif mut_type == "I":
			pins += row[col]
		else:
			raise ValueError("An unexpected mutation was encountered!")

	return pdel # returns proportion of deletions


def compile_mutations():
	# A function that creates a dictionary of all the
	# needed mutations in indel-predictions
	mut_types = dict()
	for d in reversed(range(1,31)):
		mut_types["D{0}".format(d)] = 0
	for i in range(1,11):
		mut_types["I{0}".format(i)] = 0
	return mut_types
			

def stats_frequency_mutations(df, col, reference, mutcol = 'type'):
	# Assume that df is a pandas dataframe containing at least one 
	# column with mutations (mutcol) and at least one with the amount
	# of observed mutations (col). The type of mutation and its prevalence
	# is appended to a reference dictionary (generated with compile_mutations) 
	for index, row in df.iterrows():
		mut = row[mutcol]
		mut = mut[0:mut.find('_')]
		if mut in reference:
			reference[mut] += int(row[col])
	return reference


	





