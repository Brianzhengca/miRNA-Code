#Linear Regression to Find miRNAs that are expressed at different levels across stages
import pandas as pd
import numpy as np
import scipy.stats as stats
from statsmodels.sandbox.stats.multicomp import multipletests


df = pd.read_csv("colorectaltoalldesign", sep="\t", index_col=0)
exp = pd.read_csv("colorectaltoallexp", sep="\t", index_col=0)
def finddiff(mirna):
    running = []
    colorectal = []
    for index, item in enumerate(df["colorectal cancer"]):
        if int(item) == 0:
            running.append(exp.at[mirna, df.index[index]])
        else:
            colorectal.append(exp.at[mirna, df.index[index]])
    running = np.array(running)
    colorectal = np.array(colorectal)
    if np.var(running) > np.var(colorectal):
        if np.var(running) > 4*np.var(colorectal):
            print("hi")
            return stats.ttest_ind(a=colorectal, b=running, equal_var=False)
        else:
            return stats.ttest_ind(a=colorectal, b=running, equal_var=True)
    elif np.var(running) < np.var(colorectal):
        if np.var(running) < 0.25*np.var(colorectal):
            print("hii")
            return stats.ttest_ind(a=colorectal, b=running, equal_var=False)
        else:
            return stats.ttest_ind(a=colorectal, b=running, equal_var=True)
cnt = 0
dt = {}
p_list = []
gene_list = []
for gene in [line.rstrip() for line in open('sig_genes.txt')]:
    stat, p = finddiff(gene)
    p_list.append(p)
    gene_list.append(gene)
p_list = multipletests(p_list, method="bonferroni")
p_list = p_list[1].tolist()
for index, value in enumerate(p_list):
    if value > 0.01:
        print(gene_list[index])
"""
print("-----")
keys = list(dt.keys())
keys.sort()
for i in range(20):
    print(dt[keys[i]])
"""