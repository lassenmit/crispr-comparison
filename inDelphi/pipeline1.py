
from runmodel import runDelphi
import pandas as pd

target = "/home/heymann/Dropbox/Msc/MIT/Courses/20.490 Foundations of Computational and Systems Biology/Project/Data/sampled_oligos.txt"
gt_target = "/home/heymann/Dropbox/Msc/MIT/Courses/20.490 Foundations of Computational and Systems Biology/Project/Data/SelfTarget Mutational Profiles/ST_June_2017_CHO_LV7A_DPI7/ST_June_2017_CHO_LV7A_DPI7/"
oligos =  pd.read_csv(target, sep ="\t", header = None) 

for cType in ['mESC']:
	runDelphi(oligos, gt_target, 3000, cell_type = cType, file_prefix = "run09_test_CHO", pathout = "../Derived/")
 
