
import pandas as pd
from dlutils import sample_oligos, predict_and_collect_ground_truth, available_oligos

# Setup variables and targets
pd.set_option('display.max_rows', 10)
target = "/home/heymann/Dropbox/Msc/MIT/Courses/20.490 Foundations of Computational and Systems Biology/Project/Data/sampled_oligos.txt"
gt_target = "/home/heymann/Dropbox/Msc/MIT/Courses/20.490 Foundations of Computational and Systems Biology/Project/Data/SelfTarget Mutational Profiles/ST_June_2017_CHO_LV7A_DPI7/ST_June_2017_CHO_LV7A_DPI7/"

# Run simple analysis
oligos =  pd.read_csv(target, sep ="\t", header = None)
predict_and_collect_ground_truth(oligos, gt_target, n = 100, wfile= "losses_trial02.txt")


#target = "/home/heymann/Dropbox/Msc/MIT/Courses/20.490 Foundations of Computational and Systems Biology/Project/Data/test.txt"
#jide_available_oligos(target, gt_target, n = 100, wfile = "jide_sampled_oligos.txt")
