import numpy as np
import matplotlib.pyplot as plt
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import os
import random

def qqR2_python(corvec, nn, title="QQ Plot", path_to_save=None, color="blue"):

    mm = len(corvec)
    nullcorvec = np.tanh(np.random.normal(size=mm)/np.sqrt(nn-3))   
    
    null_r2_vec = [cor**2 for cor in nullcorvec]
    r2_vec = [cor**2 for cor in corvec]

    plt.plot(np.sort(nullcorvec),np.sort(corvec),"o", color=color)
    plt.plot([-1,1],[-1,1],"-")
    plt.xlabel("Null R")
    plt.ylabel("Empirical R")
    plt.title(title+" R vs. Null")
    plt.savefig(path_to_save+"_pearsonR.png")
    plt.clf()

    plt.plot(np.sort(null_r2_vec),np.sort(r2_vec),"o", color=color)
    plt.plot([0,1],[0,1],"-")
    plt.xlabel("Null R2")
    plt.ylabel("Empirical R2")
    plt.title(title+" R2 vs. Null")
    plt.savefig(path_to_save+"_pearsonR2.png")
    plt.clf()



def hist_pearson_R_plot(df, nn, color_mapping=None, title="Histogram of R", path_to_save=None):

    df["Null"] = np.tanh(np.random.normal(size=len(df))/np.sqrt(nn-3))
    melted_df = df.melt(value_vars=df.columns)
    melted_df

    if color_mapping is not None:
        sns.histplot(data=melted_df, x="value", hue="variable", palette=color_mapping, bins=100)
    else:
        sns.histplot(data=melted_df, x="value", hue="variable", bins=100)
    plt.grid(False)
    plt.xlabel("Pearson Correlation")
    plt.ylabel("Count")
    plt.title("Correlation distribution of EnPACT, PredictDB cv, and Null")
    plt.savefig(path_to_save)
    plt.clf()


def empirical_null_correlation(empirical, ground_truth, corr_type="pearson", axis=1):

    empirical_col_permuted = np.random.permutation(empirical.columns)

    empirical_permuted = empirical[empirical_col_permuted]

    return empirical_permuted.corrwith(ground_truth, axis=axis, method=corr_type)

def analytical_null_correlation(empirical):

    mm = empirical.shape[0]
    nn = empirical.shape[1]
    nullcorvec = np.tanh(np.random.normal(size=mm)/np.sqrt(nn-3))   
    return nullcorvec
