import io, os, csv, sys

from predictor.model import computePredictedProfile, readTheta, setFeaturesDir, setReadsDir
from predictor.features import calculateFeaturesForGenIndelFile, readFeaturesData
from predictor.predict import predictMutationsBulk, predictMutationsSingle

def predict_single(target_seq, pam_idx, output_prefix):
    # uses selfTarget to predict the outcome of a single
    # target_seq: "ATGCTTCATTCGAAAACTTGCAATAAGAGCGCACGATCCAGGCGGCTTGCAATAACCCAGTCCGCTCGACGACCTCTGC"
	# pam_idx = 42 
    # output_prefix = "test"
    
    if sum([x not in 'ATGC' for x in target_seq]) > 0:
        raise Exception('Invalid target sequence, expecting string containing only A,T,G,C:\n%s' % target_seq)
    predictMutationsSingle(target_seq, pam_idx, output_prefix)
    
def predict_batch(batch_file, output_prefix):
	# A function that predicts a batch file of PAMS
    if not os.path.isfile(batch_file):
        raise Exception('Count not find batch file ' + batch_file)    
    predictMutationsBulk(batch_file, output_prefix)

def find_pam(seq):
	# A function that returns the index(s)
	# for different Pams
    pams = ['AGG','TGG','CGG','GGG']
    all_indexes = []
    for index in range(len(seq)):
        codon = seq[index:(index+3)]
        if codon in pams:
            all_indexes.append(index)
    return all_indexes
