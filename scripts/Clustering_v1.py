"""
Credits:
- Freek de Kreek, 2022

Description:
Lorum ipsum dolor sit amet consectetuer

Changes:
@date-of-change, description of added feature/edit

"""

# Import packages here
import argparse
import timeit
from collections import defaultdict
from itertools import combinations_with_replacement

import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from sklearn.metrics import adjusted_rand_score
import seaborn as sns

# Initialize parser, for all arguments given in command-line
parser = argparse.ArgumentParser()
# Main arguments
parser.add_argument("-i", "--Input",
                    help="Provide input file name of the distance matrix in '.xlsx' format, see help for other input "
                         "options", required=False)
parser.add_argument("-o", "--Output", help="Provide path to the output folder", required=False)

# Optional arguments
parser.add_argument("-v", "--Verbose", help="When this argument is given, it will turn on debug, this will print A LOT",
                    action="store_true")

# Variables
args = parser.parse_args()
main_input = args.Input
main_output = args.Output
verbose = args.Verbose
start = timeit.default_timer()


######################################## Functions #####################################################################

def cluster_analysis(df_ca, dist_ca, df_clst_label):
    """
    Put a function description here
    :param df_ca, is the dataframe with all results
    :param dist_ca, is the name or abbreviation of the genetic distance
    :param df_clst_label, dataframe to which the labels of the different clusters are assigned to
    :return:
    """

    df_ca1 = df_ca[["query", "reference", dist_ca]]
    cutoff_dist_dict = {"snp": {"lowest_cutoff_r_r": 120, "highest_cutoff_ur_ur": 44, "max_acc": 44}, "pmatch":
        {"lowest_cutoff_r_r": 0.081, "highest_cutoff_ur_ur": 0.0256, "max_acc": 0.0339}, "wgMLST_chewBACCA":
                            {"lowest_cutoff_r_r": 561, "highest_cutoff_ur_ur": 208, "max_acc": 237}, "wgMLST_seqsphere":
                            {"lowest_cutoff_r_r": 128, "highest_cutoff_ur_ur": 42, "max_acc": 42},
                        "wgMLST_stable_seqsphere":
                            {"lowest_cutoff_r_r": 97, "highest_cutoff_ur_ur": 34, "max_acc": 34},
                        "cgMLST_stable_seqsphere":
                            {"current_cutoff": 29}}
    cutoff_dict = cutoff_dist_dict[dist_ca]
    sl_cluster_dict = {}
    df_cluster_analysis = pd.DataFrame()
    # list_unique = df_ca1['query'].to_list()
    # list_unique.extend(df_ca1['reference'].to_list())

    for cutoff, cutoff_dist in cutoff_dict.items():
        # Get the samples that are within the cutoff
        df_cluster = df_ca1[df_ca1[dist_ca] <= cutoff_dist]
        # list_unique = df_cluster['query'].to_list()
        # list_unique.extend(df_cluster['reference'].to_list())
        # list_unique = set(list_unique)

        # Create a list to store all pairs of samples within the current clustercutoff
        set_list = []
        for row in df_cluster.itertuples():
            set_list.append([row.query, row.reference])

        def merge_clusters(pairs):
            merged = False
            for i in range(len(pairs)):
                for j in range(i + 1, len(pairs)):
                    if any(sample in pairs[i] for sample in pairs[j]):
                        merged_pair = list(set(pairs[i] + pairs[j]))
                        merged_pairs = merge_clusters([merged_pair] + pairs[:i] + pairs[i + 1:j] + pairs[j + 1:])
                        return merged_pairs
            return pairs

        sl_clusters = merge_clusters(set_list)
        sl_cluster_dict[cutoff] = sl_clusters

        # Calculate all cluster parameters
        n_pws_comp = df_cluster.shape[0]
        n_unique = [s for pair in sl_clusters for s in pair]
        n_clusters = len(sl_clusters)
        len_list = [len(clst) for clst in sl_clusters]
        # print(dist_ca, sum(len_list))
        median_in_cluster = np.median(len_list)
        median_in_cluster = str(median_in_cluster).replace('.', ',')
        range_n_in_cluster = [min(len_list), max(len_list)]
        perc_in_cluster = len(n_unique) / 70 * 100
        perc_in_cluster = str(perc_in_cluster).replace('.', ',')
        n_sim = [len(n) * (len(n) - 1) for n in sl_clusters]  # sum of number of samples * (n samples -1)
        N_sim = 70 * 69  # total number of samples * (n samples - 1)
        sim_index = 1 - sum(n_sim) / N_sim
        sim_index = str(sim_index).replace('.', ',')
        sil_score = silhouette_score_calculator(sl_clusters, df_ca1, dist_ca, cutoff)
        sil_score = str(sil_score).replace('.', ',')

        # Note the number of clusters, n samples, etc
        df_clst = pd.DataFrame()
        df_clst['typing_tool'] = pd.Series(dist_ca)
        df_clst[cutoff] = cutoff_dist
        df_clst[f'{cutoff}_n_cluster'] = n_clusters
        df_clst[f'{cutoff}_n_in_cluster'] = len(n_unique)
        df_clst[f'{cutoff}_median_in_cluster'] = median_in_cluster
        df_clst[f'{cutoff}_range_n_in_cluster'] = str(range_n_in_cluster)
        df_clst[f'{cutoff}_perc_in_cluster'] = perc_in_cluster
        df_clst[f'{cutoff}_sim_index'] = sim_index
        df_clst[f'{cutoff}_sil_score'] = sil_score
        df_cluster_analysis = pd.concat([df_cluster_analysis, df_clst])

        print(f'{dist_ca}:{cutoff}:n_cluster:{n_clusters}:n_in_cluster:{len(n_unique)}:median_in_cluster:'
              f'{median_in_cluster}:range_n_in_cluster:{str(range_n_in_cluster)}:perc_in_cluster:{perc_in_cluster}:'
              f'sim_index:{sim_index}:sil_score:{sil_score}')

        def cluster_labelling():
            # Create a new column to note clusters
            new_col_name = f'{dist_ca}_{cutoff}'
            df_ca.assign(new_col_name='')

            # Create a cluster label for each sample
            cluster_dict = {}
            for i, cluster in enumerate(sl_clusters):
                for sample in cluster:
                    cluster_dict[sample] = i

            reversed_cluster_dict = defaultdict(list)
            for key, value in cluster_dict.items():
                reversed_cluster_dict[value].append(key)

            # Give each pairwise comparison in a cluster a label
            for pair in sl_clusters:
                for sample in pair:
                    df_ca.loc[(df_ca['query'].isin(reversed_cluster_dict[cluster_dict[sample]])) & (df_ca['reference'].isin(reversed_cluster_dict[cluster_dict[sample]])), new_col_name] = cluster_dict[sample]

            # Create overview of metadata of samples in a cluster
            new_column_name = f'{dist_ca}_{cutoff}'
            df_clst_label.assign(new_column_name='')

            for sample, cluster_id in cluster_dict.items():
                df_clst_label.loc[df_clst_label['Sample ID'] == sample, new_column_name] = cluster_id
            return
        cluster_labelling()
    return df_cluster_analysis, sl_cluster_dict


def cluster_comp(clst_list1, clst_list2, df_clst, dist1, dist2):
    """
    Compare the single linkage clusters of one typing method with another by calculating the Rand index.
    Rand index = (a + b)/n(n-1)/2
    a, the number of pairs of elements in S that are in the same subset in X and in the same subset in Y
    b, the number of pairs of elements in S that are in different subsets in X and in different subsets in Y
    """
    # Get the similarities and differences
    jaccard_dict = {}
    cutoff_list = ['lowest_cutoff_r_r', 'highest_cutoff_ur_ur', 'max_acc']
    ari_dict = {}
    all_samples = set(df_clst['query'].tolist() + df_clst['reference'].tolist())

    for cutoff1 in clst_list1.keys():
        for cutoff2 in cutoff_list:
            if (dist1 != 'cgMLST_stable_seqsphere') & (dist2 != 'cgMLST_stable_seqsphere'):
                clusters1 = clst_list1[cutoff1]
                clusters2 = clst_list2[cutoff2]
            elif (dist1 == 'cgMLST_stable_seqsphere') & (dist2 != 'cgMLST_stable_seqsphere'):
                clusters1 = clst_list1['current_cutoff']
                clusters2 = clst_list2[cutoff2]
            elif (dist1 != 'cgMLST_stable_seqsphere') & (dist2 == 'cgMLST_stable_seqsphere'):
                clusters1 = clst_list1[cutoff1]
                clusters2 = clst_list2['current_cutoff']
            else:
                clusters1 = clst_list1['current_cutoff']
                clusters2 = clst_list2['current_cutoff']

            # Add any missing samples to the cluster list
            missing_clusters1 = list(all_samples.difference(set([sample for cluster in clusters1 for sample in cluster])))
            missing_clusters2 = list(all_samples.difference(set([sample for cluster in clusters2 for sample in cluster])))

            for sample in missing_clusters1:
                if sample not in [s for cluster in clusters1 for s in cluster]:
                    clusters1.append([sample])
            for sample in missing_clusters2:
                if sample not in [s for cluster in clusters2 for s in cluster]:
                    clusters2.append([sample])

            # Now we can calculate the Rand index
            def adjusted_rand_index(cl1, cl2):
                """
                Calculates the Rand index between two clusters cl1 and cl2.
                the function checks for each sample in cl1 whether it is in cl2 and vice versa. It then calculates the
                number of pairs that are in the same cluster in both cl1 and cl2 (a), the number of pairs that are in
                different clusters in both cl1 and cl2 (b).the number of pairs that are in the same cluster in cl1 but
                in different clusters in cl2 (c), the number of pairs that are in different clusters in cl1 but in the
                same cluster in cl2 (d)
                """
                n = len(cl1)
                m = len(cl2)
                a = b = c = d = 0
                for i in range(n):
                    for j in range(i + 1, n):
                        try:
                            if (cl1[i] == cl1[j]) and (cl2[i] == cl2[j]):
                                a += 1
                            elif (cl1[i] != cl1[j]) and (cl2[i] != cl2[j]):
                                b += 1
                            elif (cl1[i] == cl1[j]) and (cl2[i] != cl2[j]):
                                c += 1
                            else:
                                d += 1
                        except IndexError:  # If a cluster in cl1 only consists of 1 sample
                            try:
                                if (cl1[0] == cl1[0]) and (cl2[i] == cl2[j]):
                                    a += 1
                                elif (cl1[0] != cl1[0]) and (cl2[i] != cl2[j]):
                                    b += 1
                                elif (cl1[0] == cl1[0]) and (cl2[i] != cl2[j]):
                                    c += 1
                                else:
                                    d += 1
                            except IndexError:  # If a cluster in cl2 only consists of 1 sample
                                try:
                                    if (cl1[i] == cl1[j]) and (cl2[0] == cl2[0]):
                                        a += 1
                                    elif (cl1[i] != cl1[j]) and (cl2[0] != cl2[0]):
                                        b += 1
                                    elif (cl1[i] == cl1[j]) and (cl2[0] != cl2[0]):
                                        c += 1
                                    else:
                                        d += 1
                                except IndexError:  # If both clusters in cl1 & cl2 consist of 1 sample
                                    if (cl1[0] == cl1[0]) and (cl2[0] == cl2[0]):
                                        a += 1
                                    elif (cl1[0] != cl1[0]) and (cl2[0] != cl2[0]):
                                        b += 1
                                    elif (cl1[0] == cl1[0]) and (cl2[0] != cl2[0]):
                                        c += 1
                                    else:
                                        d += 1
                return (a + b) / (a + c + d + b)

            def adjusted_rand_index4(cl1, cl2):
                """
                To use the adjusted_rand_score function, two lists of cluster labels are needed, where each
                label corresponds to the cluster that the corresponding data point belongs to. A unique label is
                assigned to each cluster in each clustering, and then the lists of clusters is flattened into a
                single list of labels for each clustering. The expression [i+1 for i in range(len(cl1)) for j in cl1[i]]
                generates a list of labels for the first clustering, where the first cluster gets label 1, the second
                cluster gets label 2, and so on.
                :param cl1:
                :param cl2:
                :return: adjusted rand index
                """
                ari = adjusted_rand_score([i + 1 for i in range(len(cl1)) for j in cl1[i]],
                                          [i + 1 for i in range(len(cl2)) for j in cl2[i]])
                print(f'{ari = } for {cutoff1} & {cutoff2}')
                return ari

            ari = adjusted_rand_index4(clusters1, clusters2)
            if (dist1 != 'cgMLST_stable_seqsphere') & (dist2 != 'cgMLST_stable_seqsphere'):
                ari_dict[f'{dist1}_{cutoff1}-{dist2}_{cutoff2}'] = ari
            elif (dist1 == 'cgMLST_stable_seqsphere') & (dist2 != 'cgMLST_stable_seqsphere'):
                ari_dict[f'{dist1}_current_cutoff-{dist2}_{cutoff2}'] = ari
            elif (dist1 != 'cgMLST_stable_seqsphere') & (dist2 == 'cgMLST_stable_seqsphere'):
                ari_dict[f'{dist1}_{cutoff1}-{dist2}_current_cutoff'] = ari
            else:
                ari_dict[f'{dist1}_current_cutoff-{dist2}_current_cutoff'] = ari
    return ari_dict


def silhouette_score_calculator(clusters, df_ss, dist_ss, cutoff):
    """
    :param clusters, list of clusters (list of lists contains the samples within a cluster)
    :param df_ss, dataframe of all samples, only contains the columns query, reference and a distance
    :param dist_ss, corresponding genetic distances of the clusters
    :param cutoff, corresponding cutoff which the clusters were based on
    :returns: overall silhouette score of all clusters within the cluster list

    First, the list of clusters is converted to a dictionary that maps each sample to its cluster number.
    This is done by using a loop that iterates over the list of clusters.
        Next, the pairwise distance is calculated between each sample and all other samples in the same cluster, as well
    as the distance between the sample and all samples in the nearest neighboring cluster. The pairwise distance matrix
    is used for this, and two dictionaries are created to store the results.
        The intra-cluster and nearest-cluster distances are used to calculate the silhouette coefficient for each sample
    and then the average of each cluster is taken.
        Finally, the overall silhouette score is calculated for the clustering solution by taking the average of the
    silhouette scores for all clusters:
    """
    unique_samples = [x for sublist in clusters for x in sublist]
    unique_samples = list(set(unique_samples))

    df_ss = df_ss[['query', 'reference', dist_ss]]
    df_ss = df_ss[df_ss['query'].isin(unique_samples) & df_ss['reference'].isin(unique_samples)]

    # Create a cluster label for each sample
    cluster_dict = {}
    for i, cluster in enumerate(clusters):
        for sample in cluster:
            cluster_dict[sample] = i

    # Calculate the intra and nearest cluster distances for each sample
    intra_cluster_distances = {}
    nearest_cluster_distances = {}

    for i, row in df_ss.iterrows():
        sample1 = row['query']
        sample2 = row['reference']
        distance = row[dist_ss]

        if cluster_dict[sample1] == cluster_dict[sample2]:
            # Samples are in the same cluster
            if sample1 not in intra_cluster_distances:
                intra_cluster_distances[sample1] = []
            if sample2 not in intra_cluster_distances:
                intra_cluster_distances[sample2] = []
            intra_cluster_distances[sample1].append(distance)
            intra_cluster_distances[sample2].append(distance)
        else:
            # Samples are in different clusters
            if sample1 not in nearest_cluster_distances:
                nearest_cluster_distances[sample1] = []
            if sample2 not in nearest_cluster_distances:
                nearest_cluster_distances[sample2] = []
            nearest_cluster_distances[sample1].append(distance)
            nearest_cluster_distances[sample2].append(distance)

    # Calculate the silhouette score for each sample: (b - a) / max(a, b)
    # a= average intra-cluster distance i.e the average distance between each point within a cluster.
    # b= average inter-cluster distance i.e the average distance between all clusters.
    silhouette_scores = []
    # for sample
    sil_scatter_score = []
    sil_scatter_cluster = []
    # for cluster
    sil_scatter_score_c = []
    sil_scatter_cluster_c = []
    for cluster in clusters:
        cluster_score = 0
        for sample in cluster:
            a = sum(intra_cluster_distances[sample]) / len(intra_cluster_distances[sample])
            b = min([sum(nearest_cluster_distances[sample]) / len(nearest_cluster_distances[sample]) for k, v in
                     intra_cluster_distances.items() if k != sample])
            s = (b - a) / max(a, b)
            cluster_score += s
            # for plotting
            sil_scatter_score.append(s)
            sil_scatter_cluster.append(cluster_dict[sample])
        # Calculate the average score of each cluster
        cluster_score /= len(cluster)
        silhouette_scores.append(cluster_score)

    # Calculate the overall silhouette_score
    overall_silhouette_score = sum(silhouette_scores) / len(silhouette_scores)

    plot_wot_mot = False
    if plot_wot_mot:
        sns.boxplot(x=sil_scatter_cluster, y=sil_scatter_score, medianprops={"linewidth": 1, "color": "red"}, width=0.2)
        plt.title(f'{dist_ss} {cutoff}: Silhouette score for individual samples and clusters')
        plt.xlabel('Cluster id')
        plt.ylabel('Silhouette score')
        plt.savefig(f'../results/BRMO_20230202_v2/cluster_analysis/{dist_ss}_{cutoff}_sil_score.png', bbox_inches='tight')
        plt.clf()
    return overall_silhouette_score


def dist_mtrx_converter(ari_dict):
    df_adjusted_rand_index = pd.DataFrame.from_dict(ari_dict, orient='index')
    df_adjusted_rand_index.reset_index(inplace=True)
    df_adjusted_rand_index.rename(columns={'index': 'comb1', 0: 'adjusted_rand_index'}, inplace=True)
    df_adjusted_rand_index[['comb1', 'comb2']] = df_adjusted_rand_index['comb1'].str.split('-', expand=True)

    # This creates a minimalistic heatmap, easier to read
    # df_adjusted_rand_index.loc[(df_adjusted_rand_index['comb1'].isin(['pmatch_lowest_cutoff_r_r', 'snp_lowest_cutoff_r_r'
    # , 'wgMLST_seqsphere_lowest_cutoff_r_r', 'wgMLST_stable_seqsphere_lowest_cutoff_r_r',
    # 'wgMLST_chewBACCA_lowest_cutoff_r_r', 'snp_highest_cutoff_ur_ur', 'wgMLST_seqsphere_highest_cutoff_ur_ur',
    # 'wgMLST_stable_seqsphere_highest_cutoff_ur_ur'])) | (df_adjusted_rand_index['comb2'].isin(['pmatch_lowest_cutoff_r_r', 'snp_lowest_cutoff_r_r'
    # , 'wgMLST_seqsphere_lowest_cutoff_r_r', 'wgMLST_stable_seqsphere_lowest_cutoff_r_r',
    # 'wgMLST_chewBACCA_lowest_cutoff_r_r', 'snp_highest_cutoff_ur_ur', 'wgMLST_seqsphere_highest_cutoff_ur_ur',
    # 'wgMLST_stable_seqsphere_highest_cutoff_ur_ur'])), ['comb1']] = np.nan
    # df_adjusted_rand_index.dropna(inplace=True)

    # Create distance matrix
    idx = sorted(list(set(df_adjusted_rand_index['comb1'].tolist() + df_adjusted_rand_index['comb2'].tolist())))
    df_dist_mtrx = pd.DataFrame(index=idx, columns=idx)

    for i, row in df_adjusted_rand_index.iterrows():
        # Get the values from the columns
        c1 = row['comb1']
        ri = row['adjusted_rand_index']
        c2 = row['comb2']

        # Add the value to the distance matrix
        df_dist_mtrx.loc[c1, c2] = ri
        df_dist_mtrx.loc[c2, c1] = ri

    # Plot heatmap
    plt.figure(figsize=(20, 20))
    sns.set(font_scale=2)
    ari_heatmap = sns.heatmap(df_dist_mtrx.astype('float'), annot=True, vmin=-0.01, vmax=1, cmap='coolwarm', square=True, annot_kws={'fontsize': 14})
    plt.title('Adjusted Rand Index')
    ari_heatmap.get_figure().savefig('../results/BRMO_20230202_v2/cluster_analysis/ari_heatmap.png', bbox_inches='tight')

    # Plot a pairplot
    plt.clf()
    plt.figure(figsize=(20, 20))
    ari_pairplot = sns.pairplot(df_dist_mtrx.astype('float'))
    plt.title('Adjusted Rand Index')
    ari_pairplot.fig.savefig('../results/BRMO_20230202_v2/cluster_analysis/ari_pairplot.png')
    return df_dist_mtrx


######################################## Main ##########################################################################
# Call functions here
distances = ["pmatch", "snp", "wgMLST_seqsphere", "wgMLST_stable_seqsphere", "wgMLST_chewBACCA",
             "cgMLST_stable_seqsphere"]
df_in = pd.read_excel('../results/BRMO_20230202_v2/db_results_20230202_v2.xlsx')
df_in['pmatch'] = 1 - df_in['pmatch']
df_pws_clusters = pd.DataFrame()
cluster_list_dict = {}
df_clst_label = pd.read_excel('../results/BRMO_20230202_v2/20230202_BRMO_ST131_metadata.xlsx')

# Getting the clusters and assigning labels to the clusters
for dist in distances:
    pws_clusters, cluster_list_dict[dist] = cluster_analysis(df_in, dist, df_clst_label)
    df_pws_clusters = pd.concat([df_pws_clusters, pws_clusters])

df_clst_label.to_excel('../results/BRMO_20230202_v2/cluster_analysis/cluster_overview.xlsx')
df_pws_clusters.to_excel('../results/BRMO_20230202_v2/cluster_analysis/cluster_analysis.xlsx')
print(df_pws_clusters)
df_in.to_excel('../results/BRMO_20230202_v2/cluster_analysis/df_cluster_labels.xlsx')

# Calculating the Rand index
pairwise_comparisons = list(combinations_with_replacement(distances, 2))
adjusted_rand_index_dict = {}

for comb1, comb2 in pairwise_comparisons:
    cluster_comparison = cluster_comp(clst_list1=cluster_list_dict[comb1], clst_list2=cluster_list_dict[comb2], df_clst=df_in, dist1=comb1, dist2=comb2)
    adjusted_rand_index_dict.update(cluster_comparison)

# Converting the dictionary (pairwise comparison) to a distance matrix and plotting a heatmap
print(adjusted_rand_index_dict)
df_ari = dist_mtrx_converter(adjusted_rand_index_dict)
df_ari.to_excel('../results/BRMO_20230202_v2/cluster_analysis/adjusted_rand_index_dist_mtrx.xlsx')

# End of script
stop = timeit.default_timer()  # ends timer
print("Finished succesfully, time taken: ", stop - start)  # prints script run duration
