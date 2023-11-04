#Finds miRNAs that are expressed at unique levels between colorectal cancer and all other cancers. 
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import statsmodels.formula.api as smf
from statsmodels.sandbox.stats.multicomp import multipletests

plt.figure(dpi=1200) 

def plot(mirna):

    df = pd.read_csv("seriesmatrix.txt", sep="\t")
    exp = pd.read_csv("stagesexp", sep="\t")
    exp.index = exp["ID_REF"]
    df = df.T
    df.columns = df.iloc[0]
    df = df.drop(df.index[0])
    df["state"] = [' '.join(item.split()[2:]) for item in df["state"].tolist()]
    df["sex"] = [' '.join(item.split()[1:]) for item in df["sex"].tolist()]
    df["age"] = [' '.join(item.split()[1:]) for item in df["age"].tolist()]
    df['stage'] = [' '.join(item.split()[1:]) for item in df["stage"].tolist()]
    coordinate_list = []
    coordinate0 = []
    coordinate1 = []
    coordinate2 = []
    coordinate3 = []
    coordinate4 = []
    x = []
    y = []

    for index, value in enumerate(df['state'].tolist()):
        line = df.iloc[index].tolist()
        if value == "colorectal cancer":
            if 'NA' not in line:
                if line[4] == "0":
                    x.append(0)
                    y.append(exp.at[mirna, df.index[index]])
                    coordinate0.append(exp.at[mirna, df.index[index]])
                elif line[4] == "1":
                    x.append(1)
                    y.append(exp.at[mirna, df.index[index]])
                    coordinate1.append(exp.at[mirna, df.index[index]])
                elif line[4] == "2":
                    x.append(2)
                    y.append(exp.at[mirna, df.index[index]])
                    coordinate2.append(exp.at[mirna, df.index[index]])                  
                elif line[4] == "3":
                    x.append(3)
                    y.append(exp.at[mirna, df.index[index]])
                    coordinate3.append(exp.at[mirna, df.index[index]])
                elif line[4] == "4":
                    x.append(4)
                    y.append(exp.at[mirna, df.index[index]])
                    coordinate4.append(exp.at[mirna, df.index[index]])
    
    coordinates = [[0, sum(coordinate0)/len(coordinate0)], [1, sum(coordinate1)/len(coordinate1)], [2, sum(coordinate2)/len(coordinate2)], [3, sum(coordinate3)/len(coordinate3)], [4, sum(coordinate4)/len(coordinate4)]]
    
    #plt.scatter(*zip(*coordinates))
    data = [coordinate0, coordinate1, coordinate2, coordinate3, coordinate4]
    ax = sns.violinplot(x=x, y=y)
    ax.set_xlabel("Colorectal Cancer Stage", fontsize=18)
    ax.set_ylabel("Normalized Expression Level", fontsize=18)
    #sns.regplot(x=x, y=y)
    df = pd.DataFrame([x, y])
    df = df.T
    df.columns = ["stage", "expression"]
    model = smf.ols(formula='stage ~ expression', data=df).fit()
    #plt.boxplot(data)
    fig = ax.get_figure()
    fig.savefig(mirna+'.png',bbox_inches='tight')
    plt.show()
    return (model.pvalues.loc['expression'])
plot("hsa-miR-4435")
"""
p_list = []
gene_list = []
for gene in [line.rstrip() for line in open('sig_genes_unique.txt')]:
    p_list.append(plot(gene))
    gene_list.append(gene)
p_list = multipletests(p_list, method="bonferroni")
p_list = p_list[1].tolist()
for index, value in enumerate(p_list):
    if value < 0.01:
        print(gene_list[index], value)
"""