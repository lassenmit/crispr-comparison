import pandas as pd
from dlutils import sample_oligos, predict_and_collect_ground_truth, available_oligos


pd.set_option('display.max_rows', 10)
# Run on CHO cells
## Sample the oligos
target = "/home/heymann/Dropbox/Msc/MIT/Courses/20.490 Foundations of Computational and Systems Biology/Project/Data/sampled_oligos.txt"
gt_target = "/home/heymann/Dropbox/Msc/MIT/Courses/20.490 Foundations of Computational and Systems Biology/Project/Data/SelfTarget Mutational Profiles/ST_Feb_2018_RPE1_500x_7B_DPI7_dec/ST_Feb_2018_RPE1_500x_7B_DPI7_dec/"
available_oligos(target, gt_target, n=3000, wfile = "rmme.txt")


#oligos =  pd.read_csv(target, sep ="\t", header = None) 
#predict_and_collect_ground_truth(oligos, gt_target, n = 3000, wfile= "result_RPE_trial01_3000.txt")

# Run on 


