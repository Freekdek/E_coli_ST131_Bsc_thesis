"""
Made by:
Freek de Kreek
def distancematrix by Joep Stohr

Description:
This script will run the entire epi-link pipeline
Processing the distance matrix as input, creating a pariwise comparison, adding metadata using a metadata database,
finding epi-links and categorize them, and plot these based on their frequency (count).
Currently, the epi_links are considered to be the institution, department, date (currently not taken into consideration)
and patient ID. Allowing up to 6 variables to be categorized.

Last edit: @3-2-2023: cleaning (removed define_epi) and other code, adjustments so a Excel file with multiple distances
can be given as input and produces histograms. Tried to make a different distance matrix to pairwise distance matrix...
@9-2-2023: added option to directly plot figures from merged result input
"""
import os.path
import subprocess
import numpy as np
import pandas as pd
import argparse
import seaborn as sns
import timeit
import logging
import matplotlib.pyplot as plt
from pandas.tseries.offsets import DateOffset
from scipy.spatial.distance import squareform, pdist

pd.options.mode.chained_assignment = None
start = timeit.default_timer()  # starts timer

# Initialize parser, for all arguments given in command-line
parser = argparse.ArgumentParser()
# Adding optional argument
parser.add_argument("-i", "--Input",
                    help="Provide input file name of the distance matrix in '.xlsx' format, see help for other input options",
                    required=True)
parser.add_argument("-t", "--Type", help="Provide the method", required=False)
parser.add_argument("-d", "--Database", help="Provide path to the database with metadata")
parser.add_argument("-r", "--Hist_range", help="Fill in the range of the histogram, may vary", default=5000, type=int)
parser.add_argument("-o", "--Output", help="Provide path to the output folder", required=True)
parser.add_argument("-v", "--Verbose", help="When this argument is given, it will turn on debug, this will print A LOT",
                    action="store_true")
parser.add_argument("--pws_comp",
                    help="When this argument is given, it will use a pairwise comparison in '.xlsx' format as input",
                    action="store_true")
parser.add_argument("--metadata", help="When this argument is given, it will use a metadata in '.xlsx' format as input",
                    action="store_true")
parser.add_argument("--pws_only", help="When this argument is given, it will only produce a pairwise distance matrix",
                    action="store_true")
parser.add_argument("--merged_results", help="Will process a Excel file with all results and plot histograms",
                    action="store_true")
# Read arguments from command line
args = parser.parse_args()

input = args.Input
output = args.Output
database = args.Database
data_type = args.Type
hist_range = args.Hist_range
Verbose = args.Verbose
if not os.path.isdir(output):
    os.makedirs(output)
# variables
bin_size = 250
bin_width = 10
if data_type == "pmatch":
    bin_width = 0.01
img_ext = "png"
debug = "False"
label_dict = {"pmatch": "Pairwise difference in sequence similarity (1 - pmatch)",
              "snp": "Pairwise difference in number of SNP's",
              "wgMLST_seqsphere": "Pairwise difference in number of alleles",
              "wgMLST_stable_seqsphere": "Pairwise difference in number of alleles",
              "wgMLST_chewBACCA": "Pairwise difference in number of alleles", "cgMLST_stable_seqsphere":
                  "Pairwise difference in number of alleles"}
title_dict = {"pmatch": "PopPUNK", "snp": "Snippy",
              "wgMLST_seqsphere": "E. coli ST131 specific wgMLST scheme in SeqSphere",
              "wgMLST_stable_seqsphere": "Standard E. coli wgMLST scheme in SeqSphere",
              "wgMLST_chewBACCA": "E. coli ST131 specific wgMLST scheme in chewBACCA", "cgMLST_stable_seqsphere":
              "Standard E. coli cgMLST scheme in SeqSphere"}
title_dict1 = {"pmatch": "Distribution of genetic distances using PopPUNK", "snp": "Distribution of genetic distances using snippy",
              "wgMLST_seqsphere": "Distribution of genetic distances using a specific wgMLST scheme for E. coli ST131 in SeqSphere",
              "wgMLST_stable_seqsphere": "Distribution of genetic distances using a standard wgMLST E. coli in SeqSphere",
              "wgMLST_chewBACCA": "Distribution of genetic distances using wgMLST scheme in chewBACCA", "cgMLST_stable_seqsphere":
               "Distribution of genetic distances using a standard cgMLST E. coli in SeqSphere"}

# this is for logging
logging.basicConfig(filename=f"{output}epi_{data_type}.log", format='%(asctime)s %(message)s', filemode='w')
logger = logging.getLogger()
logger.setLevel(logging.INFO)

if Verbose:
    debug = True


######################################## Functions #####################################################################


def distancematrix(file, type):
    """
    This function transforms a distance matrix to a pairwise comparison.
    :param file: distance matrix
    :param type: unit of comparisons
    :return: pairwise comparison
    """
    # read input and add query, reference and distances
    df = pd.read_excel(file)
    list_dist = []
    list_query = []
    list_ref = []
    col = df.columns

    for x in range(1, (len(col))):  # iterate every collumn, collumn index = x
        d = df[df.columns[x]][x:]
        com = df[df.columns[0]][x:]

        for item in com:  # appends all query sample IDs to a list
            if ".ffn" in item:
                item = item.replace(".ffn", "")
            list_query.append(item)
            if debug == "True":
                print(item)

        for it in d:  # appends all distances to a list
            list_dist.append(it)

        for f in range(x, (len(col) - 1)):  # appends all reference sample IDs to a list
            reference = df.columns[x]
            if ".ffn" in reference:
                reference = reference.replace(".ffn", "")
            list_ref.append(reference)
            if debug == "True":
                print(reference)

    if debug == "True":
        print(list_ref)
        print(len(list_query))
        print(len(list_ref))

    df_pws_comp = pd.DataFrame()
    df_pws_comp["query"] = list_query  # list with query IDs
    df_pws_comp["reference"] = list_ref  # list with reference IDs
    df_pws_comp[type] = list_dist  # notes the data type (SNP/MLST/etc) as name and the distances(values)
    print("Pairwise comparison OK")
    logger.info("Pairwise comparison OK")
    return df_pws_comp


def add_metadata(file, db):
    """
    In this function metadata is added by merging two Excel files, a file with sample ID's and their metadata and a file
    with pairwise comparisons.
    :param file: pairwise comparison
    :param db: metadata
    :return: pairwise comparison with metadata
    """
    # read database and add epi data of query and reference
    df_db = pd.read_excel(db, converters={"OIVcodering": str})
    df_db.rename(columns={"Collection Date": "Collection_date", "patient_ID": "patient_id"}, inplace=True)
    metadata_columns = ["Sample ID", "patient_id", "Collection_date"]
    df_db_new = df_db[metadata_columns]

    # fill empty rows in OIVcodering and split it in 3 sections (inst, loc, afd)
    df_db["OIVcodering"] = df_db.OIVcodering.fillna("")
    df_db_new["instelling"] = df_db["OIVcodering"].str[:2]
    df_db_new["locatie"] = df_db["OIVcodering"].str[2:4]
    df_db_new["afdeling"] = df_db["OIVcodering"].str[4:]

    # rename column and also correct sample id
    df_db_new.rename(columns={"Sample ID": "query"}, inplace=True)

    # read input and rename query, reference and distances
    df = file
    df.rename(columns={"Query": "query"}, inplace=True)
    df.rename(columns={"Reference": "reference"}, inplace=True)
    if debug:
        print(df_db_new.columns, df.columns)
        print(df_db_new, df)
    df_out = pd.merge(df, df_db_new, how="left", on="query")
    df_db_new.rename(columns={"query": "reference"}, inplace=True)
    df_out = pd.merge(df_out, df_db_new, how="left", on="reference", suffixes=("_query", "_reference"))
    print("Metadata OK")
    logger.info("Metadata OK")
    return df_out


def grouping(file, distance_type):
    """
    In this function comparisons are made based on metadata, to define the different epi links.
    After specifying the epi links different graphs are plotted (histograms and violin plot).
    :param file: pairwise comparison with metadata
    :param distance_type: type of unit to be compared between the pairwise comparisons
    :return: histogram, violin plots, statistics of input file
    """
    # read input and add query, reference and distances
    df = file
    df["Cdate_query"] = df["Collection_date_query"]
    df["Cdate_reference"] = df["Collection_date_reference"]

    # start comparing query and reference and assign epi link id
    df["epi_link_id"] = None
    df_pt1 = df.dropna(subset=["patient_id_query", "patient_id_reference"])
    df_pt1 = df_pt1[df_pt1["patient_id_query"] == df_pt1["patient_id_reference"]]
    df_pt1["epi_link_id"] = "0"
    df_pt = df.dropna(subset=["patient_id_query", "patient_id_reference"])
    df_pt = df_pt[df["patient_id_query"] == df["patient_id_reference"]]
    df_pt["epi_link_id"] = "1"
    df_afd1 = df.dropna(subset=["afdeling_query", "afdeling_reference"])
    df_afd1 = df_afd1[
        (df_afd1["afdeling_query"] == df_afd1["afdeling_reference"]) & (df_afd1["instelling_query"] == df_afd1["instelling_reference"])]
    df_afd1["epi_link_id"] = "0"
    df_afd = df.dropna(subset=["afdeling_query", "afdeling_reference"])
    df_afd = df_afd[
        (df_afd["afdeling_query"] == df_afd["afdeling_reference"]) & (df_afd["instelling_query"] == df_afd["instelling_reference"])]
    df_afd["epi_link_id"] = "1"
    df_inst1 = df.dropna(subset=["instelling_query", "instelling_reference"])
    df_inst1 = df_inst1[df_inst1["instelling_query"] == df_inst1["instelling_reference"]]
    df_inst1["epi_link_id"] = "0"
    df_inst = df.dropna(subset=["instelling_query", "instelling_reference"])
    df_inst = df_inst[df_inst["instelling_query"] == df_inst["instelling_reference"]]
    df_inst["epi_link_id"] = "1"

    # dropping duplicate columns, so we can join later
    drop_this = ["query", "reference", distance_type, "patient_id_query", "Collection_date_query",
                 "instelling_query", "locatie_query", "afdeling_query", "patient_id_reference",
                 "Collection_date_reference", "instelling_reference", "locatie_reference", "afdeling_reference"]
    df_pt1.drop(columns=drop_this, inplace=True)
    df_pt.drop(columns=drop_this, inplace=True)
    df_afd1.drop(columns=drop_this, inplace=True)
    df_afd.drop(columns=drop_this, inplace=True)
    df_inst1.drop(columns=drop_this, inplace=True)
    df_inst.drop(columns=drop_this, inplace=True)

    # adding all dfs in a single df
    df_final_group = df.join(df_pt1, how="left", rsuffix="_pt1")
    df_final_group = df_final_group.join(df_pt, how="left", rsuffix="_pt")
    df_final_group = df_final_group.join(df_afd1, how="left", rsuffix="_afd1")
    df_final_group = df_final_group.join(df_afd, how="left", rsuffix="_afd")
    df_final_group = df_final_group.join(df_inst1, how="left", rsuffix="_inst1")
    df_final_group = df_final_group.join(df_inst, how="left", rsuffix="_inst")

    # filling the epi-link columns with no epi link
    df_final_group["epi_link_id_pt1"] = df_final_group["epi_link_id_pt1"].fillna("0")
    df_final_group["epi_link_id_pt"] = df_final_group["epi_link_id_pt"].fillna("0")
    df_final_group["epi_link_id_afd1"] = df_final_group["epi_link_id_afd1"].fillna("0")
    df_final_group["epi_link_id_afd"] = df_final_group["epi_link_id_afd"].fillna("0")
    df_final_group["epi_link_id_inst1"] = df_final_group["epi_link_id_inst1"].fillna("0")
    df_final_group["epi_link_id_inst"] = df_final_group["epi_link_id_inst"].fillna("0")

    # join epi link ids (=sum of epi links)
    df_final_group["epi_link_id"] = df_final_group["epi_link_id_pt1"] + df_final_group["epi_link_id_pt"] + \
                                    df_final_group["epi_link_id_afd1"] + df_final_group["epi_link_id_afd"] + \
                                    df_final_group["epi_link_id_inst1"] + df_final_group["epi_link_id_inst"]

    # this switches the x-axis for a more intuitive analysis (makes it more similar to other methods other than PopPUNK)
    if distance_type == "pmatch":
        df_final_group[distance_type] = 1 - df_final_group[distance_type]

    # assign new epi link id to 1010101 (sum of epi links)
    combinations = ["000000", "000001", "000010", "000011", "000100", "000101", "000110", "000111", "001000", "001001",
                    "001010", "001011", "001100", "001101", "001110", "001111", "010000", "010001", "010010", "010011",
                    "010100", "010101", "010110", "010111", "011000", "011001", "011010", "011011", "011100", "011101",
                    "011110", "011111", "100000", "100001", "100010", "100011", "100100", "100101", "100110", "100111",
                    "101000", "101001", "101010", "101011", "101100", "101101", "101110", "101111", "110000", "110001",
                    "110010", "110011", "110100", "110101", "110110", "110111", "111000", "111001", "111010", "111011",
                    "111100", "111101", "111110", "111111"]
    combinations_new_id = {}
    annotated_combinations = {}
    new_id = 0

    for code in combinations:
        combinations_new_id[code] = new_id  # assign new id to epi links combination
        # this starts describing the new ids (just for reference)
        annotation_list = []
        if code[0] != "0":
            annotation_list.append("pt1")
        else:
            annotation_list.append("NaN")
        if code[1] != "0":
            annotation_list.append("pt")
        else:
            annotation_list.append("NaN")
        if code[2] != "0":
            annotation_list.append("afd1")
        else:
            annotation_list.append("NaN")
        if code[3] != "0":
            annotation_list.append("afd")
        else:
            annotation_list.append("NaN")
        if code[4] != "0":
            annotation_list.append("inst1")
        else:
            annotation_list.append("NaN")
        if code[5] != "0":
            annotation_list.append("inst")
        else:
            annotation_list.append("NaN")

        annotation_desc = ", ".join(annotation_list)
        annotated_combinations[new_id] = annotation_desc
        new_id += 1

    if debug == "True":
        print(combinations_new_id, "\n" + "\n" + "\n", annotated_combinations)
        for key, value in annotated_combinations.items():
            print(key, value)

    for epi_link in combinations_new_id.keys():  # here we replace the sum of epi links (1010101) to the final epi link id (0-63)
        df_final_group["epi_link_id"].replace(epi_link, combinations_new_id[epi_link], inplace=True)

    # here we replace the epi link id for a readable epi link id and group epi link ids for patient
    df_final_group["epi_link_id"].replace(annotated_combinations, regex=True, inplace=True)
    readable_epi_link = {"NaN, NaN, NaN, NaN, NaN, NaN": "BRMO surveillance", #Unrelated
                         "NaN, NaN, NaN, NaN, NaN, inst": "BRMO surveillance", #Same institution
                         "NaN, NaN, NaN, afd, NaN, inst": "BRMO surveillance", #Same department
                         "NaN, pt, NaN, NaN, NaN, NaN": "BRMO surveillance"} #Same patient

    for readable_link in readable_epi_link.keys():
        df_final_group.loc[df_final_group["epi_link_id"] == readable_link, "epi_link_id"] = readable_epi_link[
            readable_link]

    df_final_group.loc[(df_final_group["epi_link_id"] == "NaN, pt, NaN, afd, NaN, inst") | (
                df_final_group["epi_link_id"] == "NaN, pt, NaN, NaN, NaN, inst"), "epi_link_id"] = "Same patient"

    # let's drop some columns again
    drop_this1 = ["epi_link_id_pt1", "epi_link_id_pt", "epi_link_id_afd1", "epi_link_id_afd", "epi_link_id_inst1",
                  "epi_link_id_inst"]
    df_final_group.drop(columns=drop_this1, inplace=True)
    # df_final_group.to_excel(f"{output}epi_{data_type}_metadata_epi.xlsx")

    # now group for the final output
    df_final_group["epi_link_id"].replace("Same patient", "BRMO surveillance", inplace=True) #Related
    df_group = df_final_group.groupby("epi_link_id").agg(["count", "median", "mean", "min", "max", "std"])[distance_type]

    ######## Start plotting #######
    # make a nice histogram
    sns.set(font_scale=5)


    if distance_type == "pmatch":
        hist_range = 1
    elif distance_type == 'snp':
        hist_range = 10000
    else:
        hist_range = 2500

    fig_output = f"{output}epi_{distance_type}_hist.{img_ext}"
    sb_hist_out = f"{output}epi_{distance_type}_layered.{img_ext}"
    sb_out = f"{output}epi_{distance_type}_violin.{img_ext}"
    ax = df_final_group.hist(column=distance_type, by="epi_link_id", figsize=(30, 30), range=(0, hist_range),
                             bins=bin_size, color="b")  # change range to 'zoom' in or out
    try:
        if distance_type == "pmatch":
            ax[0][0].set_xlabel("1-pmatch", fontsize=30)
        else:
            ax[0][0].set_xlabel(distance_type, fontsize=30)
        ax[0][0].set_ylabel("frequency", fontsize=30)
        ax[0][0].get_figure().savefig(fig_output)
        ax[0][0].get_figure().clf()
    except TypeError:
        try:
            if distance_type == "pmatch":
                ax[0].set_xlabel("1-pmatch", fontsize=30)
            else:
                ax[0].set_xlabel(distance_type, fontsize=30)
            ax[0].set_ylabel("frequency")
            ax[0].get_figure().savefig(fig_output)
            ax[0].get_figure().clf()
        except TypeError:
            if distance_type == "pmatch":
                ax.set_xlabel("1-pmatch", fontsize=30)
            else:
                ax.set_xlabel(distance_type, fontsize=30)
            ax.set_ylabel("frequency")
            ax.get_figure().savefig(fig_output)
            ax.get_figure().clf()
        # plot layered histogram
        sb_hist = sns.histplot(data=df_final_group, x=distance_type, hue="epi_link_id",
                               multiple="layer", binwidth=bin_width, palette="CMRmap", binrange=(0, 1000))
        plt.title(title_dict[distance_type])
        plt.ylabel("frequency")
        plt.xlabel(label_dict[distance_type])
        sb_hist.get_figure().savefig(sb_hist_out)
        sb_hist.get_figure().clf()
        # plot violin plot
        sb_fig = sns.violinplot(data=df_final_group, x=distance_type, y="epi_link_id", scale="area", cut=0, inner="point")
        sb_fig.get_figure().savefig(sb_out, bbox_inches="tight")
        sb_fig.get_figure().clf()

        # df_final_group.to_excel(f"{output}{distance_type}_epilink.xlsx")
        print("Plotting OK")
        logger.info("Plotting OK")
    return df_group


######################################## Main ##########################################################################

# Define the output files names here
output_pws_comp = f"{output}epi_{data_type}_pws_comp.xlsx"
output_metadata = f"{output}epi_{data_type}_metadata.xlsx"
output_plots = f"{output}epi_{data_type}_stats.xlsx"

# this will create all the results, start position is dependent on the provided arguments (--pws_comp or --metadata)
print("Starting epi_pipeline...")
logger.info("Starting epi_pipeline...")
if args.pws_comp:
    df_pws_out = pd.read_excel(input)
    df_metadata_out = add_metadata(df_pws_out, database)
    df_group_out = grouping(df_metadata_out, data_type)
    df_metadata_out.to_excel(output_metadata)
    df_group_out.to_excel(output_plots)
elif args.metadata:  # uses a metadata file as input
    df_metadata_out = pd.read_excel(input,
                                    converters={"instelling_query": str, "locatie_query": str, "afdeling_query": str,
                                                "instelling_reference": str, "locatie_reference": str,
                                                "afdeling_reference": str, "Collection_date_query": str,
                                                "Collection_date_reference": str})
    df_group_out = grouping(df_metadata_out, data_type)
    df_group_out.to_excel(output_plots)
elif args.pws_only:
    df_pws_comp_out = distancematrix(input, data_type)
    df_pws_comp_out.to_excel(output_pws_comp)
elif args.merged_results:
    distances = ["pmatch", "snp", "wgMLST_seqsphere", "wgMLST_stable_seqsphere", "wgMLST_chewBACCA", "cgMLST_stable_seqsphere"]
    for dist in distances:
        df_group_out = grouping(pd.read_excel(input), dist)
else:  # will start from the distance matrix (the beginning) as default
    df_pws_comp_out = distancematrix(input, data_type)  # input=distance_matrix.xlsx, data_type=method + _Difference
    df_metadata_out = add_metadata(df_pws_comp_out, database)  # input=pws_comp_out, database=file with meta
    df_group_out = grouping(df_metadata_out, data_type)  # input=df_metadata_out

    df_pws_comp_out.to_excel(output_pws_comp)
    df_metadata_out.to_excel(output_metadata)
    df_group_out.to_excel(output_plots)

# end of pipeline |--_--__-|
stop = timeit.default_timer()  # ends timer
print("Finished succesfully, time taken: ", stop - start)  # prints script run duration
logger.info("Finished succesfully, time taken: ", stop - start)
