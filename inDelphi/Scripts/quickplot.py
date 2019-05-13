import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import collections  as mc
import scipy.stats
import pylab

def abline(slope, intercept):
    """Plot a line from slope and intercept"""
    axes = plt.gca()
    x_vals = np.array(axes.get_xlim())
    y_vals = intercept + slope * x_vals
    plt.plot(x_vals, y_vals, '--')

df = pd.read_csv("../Derived/run001_mESC_statistics.txt", header = None, sep = "\t")
df.columns = ["oligo_id", \
			"KL-Divergence", \
			"acc01_agr_1", \
			"acc02_agr_2", \
			"acc03_agr_3", \
			"acc04_mut", \
			"acc05_mut_loose", \
			"acc06_frameshift", \
			"stat01_frameshift_pred", \
			"stat02_frameshift_actual"]

x = df['stat01_frameshift_pred'].values.tolist()
y = df['stat02_frameshift_actual'].values.tolist()
corr = scipy.stats.pearsonr(x, y)[0]
print(corr)

plt.scatter(x, y, alpha = 0.5)
plt.title('2017 CHO CELLs (LV7A_DPI17)\n(InDelphi)')
plt.xlabel('Percent in-frame mutations (Predicted)')
plt.ylabel('Percent in-frame mutations (Measured')
abline(slope = 1, intercept = 0)

plt.text(0.54, 0.02, "n=3000".format(len(x)), style='italic')
plt.text(0.54, -0.025, "Pearson Correlation: {0}".format(round(corr,3)), style='italic')
plt.savefig("in_delphi_mutations_correlation.png")
#plt.show()
#plt.close()



xrng = [x for x in range(0,40)]
print(xrng)
data_actual = pd.read_csv("../Derived/run001_mESC_indels_frequency_actual.txt", header = 0, sep = "\t")
data_pred = pd.read_csv("../Derived/run001_mESC_indels_frequency_predicted.txt", header = 0, sep = "\t")


#fig = plt.figure(figsize=(10,30))
fig, ax = plt.subplots(figsize=(18,6))

fig.canvas.draw()
labels = list(data_actual)


## Ready data for plotting
y1 = data_actual.mean(0).values.tolist()
y2 = data_pred.mean(0).values.tolist()
y1 = [float(y)/sum(y1) for y in y1]
y2 = [float(y)/sum(y2) for y in y2]


## plot data (note: the order of this is important to properly depict legend)
plt.plot(xrng[0:30], y1[0:30], color='g', alpha = 0.8, marker = "*")
plt.plot(xrng[0:30], y2[0:30], color='orange', alpha = 0.8, marker = "*")
plt.plot(xrng[30:], y1[30:], color='g', alpha = 0.8, marker = "*")
plt.plot(xrng[30:], y2[30:], color='orange', alpha = 0.8, marker = "*")

#plt.plot(xrng, y2, color='orange', alpha = 0.8, marker = "*")

## Add title, axes and legends
plt.xlabel('Mutation')
plt.ylabel('Frequency')
plt.title('2017 CHO CELLs (LV7A_DPI17)\nMutations (InDelphi Assuming mESC)')
plt.gca().legend(('Measured','Predicted (SelfTarget)'))

## Axis ticks and interval
plt.xticks(xrng)
ax.set_xticklabels(labels)

## Add information about deletions and insertions
lines = [[(-0.5, -0.01 ), (29.5, -0.01)], [(30, -0.01 ), (39.5, -0.01)]]
lc = mc.LineCollection(lines, colors= ["black","grey"], linewidths=2)
ax.add_collection(lc)
ax.autoscale()
ax.margins(0.1)

plt.text(13.5, -0.0215, "Deletions", color = "black")
plt.text(33.5, -0.0215, "Insertions", color = "grey")


plt.savefig("inDelphi_mutation_frequencies.png")
plt.show()