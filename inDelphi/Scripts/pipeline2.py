
from runmodel import runDelphi
import pandas as pd

target = "/home/heymann/Dropbox/Msc/MIT/Courses/20.490 Foundations of Computational and Systems Biology/Project/Data/sampled_oligos_RPE.txt"
gt_target = "/home/heymann/Dropbox/Msc/MIT/Courses/20.490 Foundations of Computational and Systems Biology/Project/Data/ST_Feb_2018_RPE1_500x_7B_DPI7_dec/ST_Feb_2018_RPE1_500x_7B_DPI7_dec/"
oligos =  pd.read_csv(target, sep ="\t", header = None) 


for cType in ['mESC']:
	runDelphi(oligos, gt_target, 3000, cell_type = cType, file_prefix = "run10_RPE", pathout = "../Derived/")

