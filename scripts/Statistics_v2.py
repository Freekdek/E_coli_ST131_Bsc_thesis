"""
Credits:
- Freek de Kreek, 2022-2023

Description:
A pairwise-comparison table is taken as input with the column names query and reference, distance columns and an
epi_link_id column.
- The Spearman Rank Correlation coefficient is calculated for each distance pair.
- The Mann-Whitney U test between epi_link_id related and unrelated is calculated for each distance.
- Thresholds (distance) are iterated to determine the percentage of misclassified as related or unrelated for a certain
distance. The cutoff where the percentage of samples is misclassified as unrelated, is zero is marked as well.
- The specificity, sensitivity and accuracy is plotted with its corresponding threshold/cut-off

Changes:
@3-1-2023, Changes to getting the unrelated and related pairwise comparisons
@3-2-2023, Removed define_epi, added clustering
@7-2-2023, Changed clustering to dendrogram
@14-2-2023, Added cluster analysis
"""

# Import packages here
import argparse
import timeit

import numpy as np
import pandas as pd
from scipy.spatial.distance import squareform, pdist
from scipy.stats import spearmanr
from scipy.stats import mannwhitneyu
import seaborn as sns
import matplotlib.pyplot as plt
from itertools import combinations_with_replacement
import os
from scipy.cluster.hierarchy import linkage, fcluster, dendrogram, cophenet
import sklearn.metrics as skm
import sklearn.cluster as skc

# Initialize parser, for all arguments given in command-line
parser = argparse.ArgumentParser()
# Main arguments
parser.add_argument("-i", "--Input",
                    help="Provide input file name of the result database in '.xlsx' format, see help for other input "
                         "options", required=True)
parser.add_argument("-o", "--Output", help="Provide path to the output folder", required=True)

# Optional arguments
parser.add_argument("-v", "--Verbose", help="When this argument is given, it will turn on debug, this will print A LOT",
                    action="store_true")
parser.add_argument("-u", "--unknown_epi", help="Will assume there is no epi link and not compare medians",
                    action="store_true")

# Variables
args = parser.parse_args()
main_input = args.Input
main_output = args.Output
verbose = args.Verbose
unknown_epi = args.unknown_epi

start = timeit.default_timer()
if not os.path.isdir(main_output):
    os.makedirs(main_output)

distances = ["pmatch", "snp", "wgMLST_seqsphere", "wgMLST_stable_seqsphere", "wgMLST_chewBACCA"] #, "cgMLST_stable_seqsphere"]
# distances = ["Core", "Accessory", "pmatch", "snp", "wgMLST_seqsphere", "wgMLST_stable_seqsphere", "wgMLST_chewBACCA"]
heatmap_out_sp = f"{main_output}sp_heatmap.png"
# heatmap_out_ps = f"{main_output}ps_heatmap.png"
df_in = pd.read_excel(main_input)
pd.options.mode.chained_assignment = None

######################################## Functions #####################################################################


def spearman_calc(df_sc):
    """
    The spearman correlation coefficient is calculated between the different methods and a heatmap is plotted for
    visualisation.
    :param df_sc: database with results of pairwise comparisons
    :return:
    """
    df_sc["pmatch"] = 1 - df_sc["pmatch"]
    # Get all pairwise comparisons excluding repetitions
    pairwise_comparisons = list(combinations_with_replacement(distances, 2))
    sp_dict = {}
    label_map = {"pmatch": "PopPUNK",
                 "snp": "Snippy",
                 "wgMLST_seqsphere": "E. coli ST131 specific wgMLST scheme SeqSphere",
                 "wgMLST_stable_seqsphere": "Standard E. coli wgMLST scheme SeqSphere",
                 "wgMLST_chewBACCA": "E. coli ST131 specific wgMLST scheme chewBACCA",
                 "cgMLST_stable_seqsphere": "Standard E. coli cgMLST scheme SeqSphere"}

    for dist_combination in pairwise_comparisons:
        sp_result = []
        dist_combination = list(dist_combination)
        column1 = df_sc[dist_combination[0]]
        column2 = df_sc[dist_combination[1]]
        coefficient, p_value = spearmanr(column1, column2, nan_policy="omit")
        sp_result.append(coefficient)
        sp_result.append(p_value)
        # print(f"{dist_combination}: {p_value}")
        sp_dict[tuple(dist_combination)] = sp_result

    # plotting the heatmap here
    df_dist_only = df_sc[distances]
    sp_result = df_dist_only.corr(method="spearman")
    plt.figure(figsize=(15, 9))
    sns.set(font_scale=2)
    heatmap = sns.heatmap(sp_result, vmin=0.8, vmax=1, annot=True)
    plt.title("Spearman Correlation Coefficient")
    # Set the labels for the y-axis and x-axis using the label_map dictionary
    heatmap.set_yticklabels([label_map[label] for label in sp_result.index])
    heatmap.set_xticklabels([label_map[label] for label in sp_result.columns], rotation=90)
    heatmap.get_figure().savefig(heatmap_out_sp, bbox_inches="tight")
    heatmap.get_figure().clf()

    # now the same for pearson correlation
    # ps_result = df_dist_only.corr(method="pearson")
    # plt.figure(figsize=(15, 9))
    # heatmap = sns.heatmap(ps_result, vmin=0.8, vmax=1, annot=True)
    # plt.title("Pearson Correlation Coefficient")
    # # Set the labels for the y-axis and x-axis using the label_map dictionary
    # heatmap.set_yticklabels([label_map[label] for label in ps_result.index])
    # heatmap.set_xticklabels([label_map[label] for label in ps_result.columns], rotation=90)
    # heatmap.get_figure().savefig(heatmap_out_ps, bbox_inches="tight")
    # heatmap.get_figure().clf()
    return pd.DataFrame.from_dict(sp_dict, orient="index")


def get_classification_cutoff(df_cc, dist_cc):
    """
    This will plot a line plot showing the percentage of misclassified comparisons (related and unrelated) as a function
    of the threshold.
    The x-axis represents the threshold, and the y-axis represents the percentage of misclassified comparisons.
    :param df_cc: database with results of pairwise comparisons
    :param dist_cc: one of the distance columns
    :return: line plots with misclassified unrelated and related and a table
    """
    # Create an empty list to store the percentage of misclassified comparisons for each threshold
    misclassified_as_related_list = []
    misclassified_as_unrelated_list = []
    threshold_r_first = True
    threshold_ur_first = True
    threshold_r = -1
    threshold_ur = -1
    plot_range_dict = {"pmatch": np.arange(0, 1.01, 0.0001), "snp": range(0, 1251, 1),
                       "wgMLST_seqsphere": range(0, 1251, 1), "wgMLST_stable_seqsphere": range(0, 1251, 1),
                       "wgMLST_chewBACCA": range(0, 1251, 1)}
    label_dict = {"pmatch": "Pairwise difference in sequence similarity (1 - pmatch)",
                  "snp": "Pairwise difference in number of SNP's",
                  "wgMLST_seqsphere": "Pairwise difference in number of alleles",
                  "wgMLST_stable_seqsphere": "Pairwise difference in number of alleles",
                  "wgMLST_chewBACCA": "Pairwise difference in number of alleles"}
    title_dict = {"pmatch": "%Misclassified using PopPUNK", "snp": "%Misclassified using snippy",
                  "wgMLST_seqsphere": "%Misclassified using a specific wgMLST scheme for E. coli ST131 in SeqSphere",
                  "wgMLST_stable_seqsphere": "%Misclassified using a standard wgMLST E. coli in SeqSphere",
                  "wgMLST_chewBACCA": "%Misclassified using a specific wgMLST scheme for E. coli ST131 in chewBACCA"}
    spec_list = []
    sens_list = []
    acc_list = []
    acc_dict = {}
    threshold_list = []

    # need to adjust the range for certain distances
    if dist_cc == "pmatch":
        df_cc[dist_cc] = 1 - df_cc[dist_cc]
    # Core plot_range = np.arange(0, 0.01, 0.00005)
    # accessory plot_range = np.arange(0, 0.3, 0.0015)
    df_cc.replace("Same patient", "Related", inplace=True)
    # Iterate over a range of thresholds
    number_of_related = df_cc[df_cc["epi_link_id"] == "Related"].shape[0]
    number_of_unrelated = df_cc[df_cc["epi_link_id"] == "Unrelated"].shape[0]

    for threshold in plot_range_dict[dist_cc]:
        # Count the number of comparisons that are misclassified as either related or unrelated
        misclassified_as_related = df_cc[(df_cc["epi_link_id"] == "Unrelated") & (df_cc[dist_cc] <= threshold)].shape[0]
        misclassified_as_unrelated = df_cc[(df_cc["epi_link_id"] == "Related") & (df_cc[dist_cc] > threshold)].shape[0]

        # Calculate the percentage of misclassified comparisons
        misclassified_as_related_list.append(misclassified_as_related / number_of_unrelated * 100)
        misclassified_as_unrelated_list.append(misclassified_as_unrelated / number_of_related * 100)

        # Calculate the threshold where the %misclassfied as unrelated = 0 (100% sensitivity)
        if (misclassified_as_unrelated / number_of_related * 100) == 0 and threshold_ur_first:
            threshold_ur = threshold
            print(
                f"cutoff = {threshold_ur} and percentage = {misclassified_as_related / number_of_unrelated * 100} of {dist_cc} 0 %misclassfied as unrelated")
            perc_misc_as_related = misclassified_as_related / number_of_unrelated * 100
            threshold_ur_first = False

        # Calculate the threshold where the %misclassfied as related = 0 (100% specificity)
        if (misclassified_as_related / number_of_unrelated * 100) > 0 and threshold_r_first:
            try:
                threshold_r = threshold_list[-1]
                print(
                    f"cutoff = {threshold_r} and percentage = {misclassified_as_unrelated / number_of_related * 100} of {dist_cc} 0 %misclassfied as related")
                perc_misc_as_unrelated = misclassified_as_unrelated / number_of_related * 100
                threshold_r_first = False
            except IndexError:
                print(f"index error threshold: {threshold}")

        # plot specificity/sensitivity
        TP = number_of_related
        FP = misclassified_as_related
        TN = number_of_unrelated
        FN = misclassified_as_unrelated

        specificity = TN / (TN + FP) * 100
        sensitivity = TP / (TP + FN) * 100
        accuracy = (TP + TN) / (TP + FP + FN + TN) * 100
        spec_list.append(specificity)
        sens_list.append(sensitivity)
        acc_list.append(accuracy)
        acc_dict[accuracy] = threshold
        threshold_list.append(threshold)

    max_acc = max(acc_list)
    print(f"{dist_cc} max accuracy = {max_acc} and threshold = {acc_dict[max_acc]}")
    misclassified_as_related = df_cc[(df_cc["epi_link_id"] == "Unrelated") & (df_cc[dist_cc] <= max_acc)].shape[0]
    misclassified_as_unrelated = df_cc[(df_cc["epi_link_id"] == "Related") & (df_cc[dist_cc] > max_acc)].shape[0]
    perc_misc_as_related = misclassified_as_related / number_of_unrelated * 100
    perc_misc_as_unrelated = misclassified_as_unrelated / number_of_related * 100
    print(f'max_accuracy: {perc_misc_as_unrelated = } and {perc_misc_as_related =}')

    # Plot the results
    # sns.set(font_scale=3)
    plt.plot(plot_range_dict[dist_cc], misclassified_as_related_list, color="b", label="Misclassified as related")
    plt.plot(plot_range_dict[dist_cc], misclassified_as_unrelated_list, color="r", label="Misclassified as unrelated")

    if threshold_ur >= 0:
        plt.vlines(threshold_ur, ymin=0, ymax=100, colors="r", linestyles="dashed")
        plt.text(x=threshold_ur, y=100, s=threshold_ur, bbox={
            'facecolor': 'red', 'alpha': 0.5, 'pad': 0.3, "boxstyle": "round"})
    if threshold_r >= 0:
        plt.vlines(threshold_r, ymin=0, ymax=100, colors="b", linestyles="dashed")
        plt.text(x=threshold_r, y=100, s=threshold_r, bbox={
            'facecolor': 'blue', 'alpha': 0.5, 'pad': 0.3, "boxstyle": "round"})
    plt.vlines(acc_dict[max_acc], ymin=0, ymax=100, colors="g", linestyles="dashed")
    plt.text(x=acc_dict[max_acc], y=100, s=acc_dict[max_acc], bbox={
        'facecolor': 'green', 'alpha': 0.5, 'pad': 0.3, "boxstyle": "round"})

    # add (ax) titles and legend
    plt.title(title_dict[dist_cc])
    plt.xlabel(label_dict[dist_cc])
    plt.ylabel("Percentage of misclassified comparisons (%)")
    plt.legend(loc="center right")
    # plt.savefig(f"{main_output}{dist}.svg")
    plt.savefig(f"{main_output}{dist_cc}_misclassification.png")
    plt.clf()

    # plot specificity and sensitivity
    plt.plot(plot_range_dict[dist], spec_list, "-b", label="Specificity")
    plt.plot(plot_range_dict[dist], sens_list, "-r", label="Sensitivity")
    plt.plot(plot_range_dict[dist], acc_list, "-g", label="Accuracy")

    # plot threshold vlines here
    if threshold_ur >= 0:
        plt.vlines(threshold_ur, ymin=0, ymax=100, colors="r", linestyles="dashed")
        plt.text(x=threshold_ur, y=100, s=threshold_ur, bbox={
            'facecolor': 'red', 'alpha': 0.5, 'pad': 0.3, "boxstyle": "round"})
    if threshold_r >= 0:
        plt.vlines(threshold_r, ymin=0, ymax=100, colors="b", linestyles="dashed")
        plt.text(x=threshold_r, y=90, s=threshold_r, bbox={
            'facecolor': 'blue', 'alpha': 0.5, 'pad': 0.3, "boxstyle": "round"})
    plt.vlines(acc_dict[max_acc], ymin=0, ymax=100, colors="g", linestyles="dashed")
    plt.text(x=acc_dict[max_acc], y=80, s=acc_dict[max_acc], bbox={
        'facecolor': 'green', 'alpha': 0.5, 'pad': 0.3, "boxstyle": "round"})

    # add (ax) titles and legend
    plt.xlabel(label_dict[dist_cc])
    plt.ylabel("Sensitivity & Specificity (in %)")
    plt.legend()
    plt.title(f"Line plot {dist_cc} of sensitivity and specificity")
    plt.savefig(f"{main_output}{dist_cc}_spec_sens.png")
    plt.clf()

    # get the percentage of samples that are misclassified instead of comparisons
    # Count the number of comparisons that are misclassified as either related or unrelated
    misc_as_related = \
    df_cc[(df_cc["epi_link_id"] == "Unrelated") & (df_cc[dist_cc] <= threshold_ur)].groupby(by="query").count().shape[0]
    misc_as_unrelated = \
    df_cc[(df_cc["epi_link_id"] == "Related") & (df_cc[dist_cc] > threshold_r)].groupby(by="query").count().shape[0]

    # Calculate the percentage of misclassified comparisons
    misc_as_related = misc_as_related / 28 * 100
    misc_as_unrelated = misc_as_unrelated / 13 * 100

    print(f"unique true misclasifications as unrelated: {misc_as_unrelated}")
    print(f"unique true misclasifications as related: {misc_as_related}")

    df_cutoff = pd.DataFrame(columns=["typing method", "cutoff all related included", "cutoff all unrelated included",
                                      "% all related included", "% all unrelated included",
                                      "%misclassfied as unrelated samples", "%misclassfied as related samples"])
    df_cutoff["typing method"] = pd.Series(dist_cc)
    df_cutoff["cutoff all related included"] = threshold_ur
    df_cutoff["cutoff all unrelated included"] = threshold_r
    df_cutoff["% all related included"] = perc_misc_as_unrelated
    df_cutoff["% all unrelated included"] = perc_misc_as_related
    df_cutoff["%misclassfied as unrelated samples"] = misc_as_unrelated
    df_cutoff["%misclassified as related samples"] = misc_as_related
    return misclassified_as_related_list, misclassified_as_unrelated_list


def median_comp(df_mc, dist_mc):
    """
    The Mann-Whitney U test is a non-parametric statistical test that can be used to compare the medians of two groups
    of data. It is based on the ranks of the data rather than the actual values, and it does not assume that the data
    follows a particular distribution.
    If the p-value is less than a predetermined significance level (e.g. 0.05), it suggests that there is a significant
    difference between the medians of the two groups. If the p-value is greater than the significance level, it suggests
    that there is no significant difference between the medians of the two groups.
    :param df_mc:
    :param dist_mc:
    :return:
    """
    label_dict = {"pmatch": "Pairwise difference in sequence similarity (1 - pmatch)",
                  "snp": "Pairwise difference in number of SNP's",
                  "wgMLST_seqsphere": "Pairwise difference in number of alleles",
                  "wgMLST_stable_seqsphere": "Pairwise difference in number of alleles",
                  "wgMLST_chewBACCA": "Pairwise difference in number of alleles",
                  "cgMLST_stable_seqsphere": "Pairwise difference in number of alleles"}
    title_dict = {"pmatch": "PopPUNK", "snp": "Snippy",
                  "wgMLST_seqsphere": "Specific wgMLST scheme for E. coli ST131 in SeqSphere",
                  "wgMLST_stable_seqsphere": "Standard wgMLST E. coli in SeqSphere",
                  "wgMLST_chewBACCA": "Specific wgMLST scheme for E. coli ST131 in chewBACCA",
                  "cgMLST_stable_seqsphere": "Standard cgMLST E. coli in SeqSphere"}
    colors = ["green", "blue"]

    if not unknown_epi:
        df_mc.replace("Same patient", "Related", inplace=True)

        # Select the two groups of data for which you want to compare the medians
        group1 = df_mc[df_mc["epi_link_id"] == "Unrelated"][dist_mc]
        group2 = df_mc[df_mc["epi_link_id"] == "Related"][dist_mc]

        # Perform the Mann-Whitney U test
        statistic, p_value = mannwhitneyu(group1, group2)
        if p_value < 0.001:
            p_value = "< 0.001"
        # Add new row to dataframe
        df_mc[f"p_value_{dist_mc}"] = p_value
        df_mc_r = df_mc[df_mc["epi_link_id"] == "Related"]
        df_mc_ur = df_mc[df_mc["epi_link_id"] == "Unrelated"]
        # print(f"Median ur: {df_mc_ur[dist_mc].median()}\nMedian r: {df_mc_r[dist_mc].median()} of {dist_mc}")
    sns.set_style("whitegrid")
    df_mc['epi_link_id'].fillna('BRMO surveillance', inplace=True)
    boxplot_mc = sns.boxplot(df_mc, x=df_mc["epi_link_id"], y=df_mc[dist_mc],
                             medianprops={"linewidth": 5, "color": "red"}, width=0.2)

    # add (ax) titles and legend
    plt.title(title_dict[dist_mc])

    if not unknown_epi:
        plt.xlabel(f"Mann-Whitney U test p-value: {p_value}")

    plt.ylabel(label_dict[dist_mc])
    # plt.savefig(f"{main_output}{dist}.svg")
    # plt.savefig(f"{main_output}boxplot_{dist}.png")
    plt.savefig(f"{main_output}boxplot_plt_{dist_mc}.png", bbox_inches="tight")
    plt.clf()

    # df_comp = df_mc.groupby(by="epi_link_id")["pmatch", "snp", "wgMLST_seqsphere", "wgMLST_stable_seqsphere", "wgMLST_chewBACCA", "p_value_pmatch", "p_value_snp", "p_value_wgMLST_seqsphere", "p_value_wgMLST_stable_seqsphere", "p_value_wgMLST_chewBACCA"].mean()
    # df_comp.to_excel("median_grouped.xlsx")
    # print(f"Mann-Whitney U test {dist_mc}: {statistic}, {p_value}")

    df_mc_table = pd.DataFrame()
    if not unknown_epi:
        df_mc_table["typing_tool"] = dist_mc
        df_mc_table[""] = df_mc_r[dist_mc].median()
        df_mc_table[""] = df_mc_ur[dist_mc].median()
    return df_mc_table


def cluster_dendrogram(df_cd, dist_cd):
    """
    Creates clusters based on a distance cutoff
    :param df_cd: input dataframe
    :return: table with number of clusters at a cutoff
    """
    df_cd1 = df_cd.loc[:, ["query", "reference", dist_cd]]
    if dist_cd == 'pmatch':
        df_cd[dist_cd] = 1 - df_cd[dist_cd]
    cluster_table = pd.DataFrame()
    cutoff_dist_dict = {"snp": {"lowest_cutoff_r_r": 120, "highest_cutoff_ur_ur": 44, "max_acc": 44}, "pmatch":
        {"lowest_cutoff_r_r": 0.081, "highest_cutoff_ur_ur": 0.0256, "max_acc": 0.0339}, "wgMLST_chewBACCA":
        {"lowest_cutoff_r_r": 561, "highest_cutoff_ur_ur": 208, "max_acc": 237}, "wgMLST_seqsphere":
        {"lowest_cutoff_r_r": 128, "highest_cutoff_ur_ur": 42, "max_acc": 42}, "wgMLST_stable_seqsphere":
        {"lowest_cutoff_r_r": 97, "highest_cutoff_ur_ur": 34, "max_acc": 34}, "cgMLST_stable_seqsphere":
        {"current_cutoff": 29}}
    cutoff_dict = cutoff_dist_dict[dist_cd]

    # Create a mapping of ids to integers
    id_map = {id: i for i, id in enumerate(set(df_cd1['query'].tolist() + df_cd1['reference'].tolist()))}

    # Replace the ids with integers in the query and reference columns
    df_cd1['query'] = df_cd1['query'].map(id_map)
    df_cd1['reference'] = df_cd1['reference'].map(id_map)

    # Pivot the dataframe to create a matrix of distances
    distance_matrix = df_cd1.pivot(index='query', columns='reference', values=dist_cd).fillna(0).values
    distance_matrix = np.maximum(distance_matrix, distance_matrix.T)
    np.fill_diagonal(distance_matrix, 0)

    # Convert the pairwise distance matrix to a squareform distance matrix
    dist_mtrx = squareform(distance_matrix)

    # Calculate the linkage matrix
    linkage_matrix = linkage(df_cd[dist_cd], method='single') # , optimal_ordering=True #dist_mtrx
    # c, coph_dists = cophenet(linkage_matrix, pdist(distance_matrix))

    for k, v in cutoff_dict.items():
        # Extract the cluster labels at the cutoff distance
        cluster_labels = fcluster(linkage_matrix, v, criterion='distance')

        # Calculate the silhouette score
        # silhouette_score = skm.silhouette_score(distance_matrix, cluster_labels, metric='precomputed')

        # Get the number of unique clusters
        n_clusters = len(set(cluster_labels))

        # AgglomerativeClustering
        agg_cluster = skc.AgglomerativeClustering(n_clusters=None, distance_threshold=v, metric='precomputed', linkage='single').fit(distance_matrix)

        # Note the number of clusters
        cluster_table['typing_tool'] = pd.Series(dist_cd)
        # cluster_table['cophenet coefficient'] = c
        cluster_table[k] = v
        cluster_table[f'{k}_n_clusters'] = n_clusters
        cluster_table[f'{k}_n_clusters'] = agg_cluster.n_clusters_
        # cluster_table[f'{k}_silhouette_score'] = silhouette_score

        # calculate full dendrogram
        plt.figure(figsize=(25, 10))
        plt.title(f'Hierarchical Clustering Dendrogram {dist_cd} for {k} at {v}')
        plt.xlabel('sample index')
        plt.ylabel('distance')

        trunk = False
        if trunk:
            fancy_dendrogram(
                linkage_matrix,
                leaf_rotation=90,
                leaf_font_size=12,
                max_d=v,
                show_contracted=True,
                truncate_mode='lastp',
                p=12
            )
        else:
            fancy_dendrogram(
                linkage_matrix,
                leaf_rotation=90,
                leaf_font_size=12,
                max_d=v,
                annotate_above=100000
            )
        plt.savefig(f'{main_output}{dist_cd}_{k}dendro.png')
        plt.clf()
    return cluster_table


def fancy_dendrogram(*args, **kwargs):
    max_d = kwargs.pop('max_d', None)
    if max_d and 'color_threshold' not in kwargs:
        kwargs['color_threshold'] = max_d
    annotate_above = kwargs.pop('annotate_above', 0)

    ddata = dendrogram(*args, **kwargs)

    if not kwargs.get('no_plot', False):
        # plt.title('Hierarchical Clustering Dendrogram (truncated)')
        # plt.xlabel('sample index or (cluster size)')
        # plt.ylabel('distance')
        for i, d, c in zip(ddata['icoord'], ddata['dcoord'], ddata['color_list']):
            x = 0.5 * sum(i[1:3])
            y = d[1]
            if y > annotate_above:
                plt.plot(x, y, 'o', c=c)
                plt.annotate("%.3g" % y, (x, y), xytext=(0, -5),
                             textcoords='offset points',
                             va='top', ha='center')
        if max_d:
            plt.axhline(y=max_d, c='k', linestyle='--')
            # xmin, xmax, ymin, ymax = plt.axis()
            # plt.annotate(f"cutoff at: {max_d}", xy=(0, ymax), bbox=dict(boxstyle='round, pad=0.3'))
    return ddata


######################################## Main ##########################################################################
# Call functions here
sp_calc = spearman_calc(pd.read_excel(main_input))
sp_calc.to_excel(f"{main_output}stats_test.xlsx")
df_cutoff_out = pd.DataFrame()
df_mc_out = pd.DataFrame()

for dist in distances:
    if not unknown_epi:
        cutoff_table = get_classification_cutoff(df_in, dist)
        # df_cutoff_out = pd.concat([df_cutoff_out, cutoff_table])

    # clustering and plotting dendrograms
    # clustering_df = cluster_dendrogram(df_in, dist)
    # median comparison and plotting boxplots
    median_comp_table = median_comp(df_in, dist)
    df_mc_out = pd.concat([df_mc_out, median_comp_table])

print("random seed =", np.random.get_state()[1][0])

# End of script
stop = timeit.default_timer()  # ends timer
print("Finished succesfully, time taken: ", stop - start)  # prints script run duration
