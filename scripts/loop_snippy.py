# Created by Freek de Kreek, 3-11-2022
# run snippy command: snippy --ctgs snippy_input.fa --reference reference_assemblies/st131_assemblies/assembly01.fasta --cpus 10 --ram 24 --report --outdir results/snippy_results/

import os
import argparse
import subprocess
import timeit

start = timeit.default_timer()  # starts timer
# For each file:

# initialize parser
parser = argparse.ArgumentParser()
# adding optional argument
parser.add_argument("-i", "--Input", help="Show Input sample reports directory")
parser.add_argument("-o", "--Output", help="Show Output directory")
parser.add_argument("-v", "--Verbose", help="Show more information")
# read arguments from command line
args = parser.parse_args()
input = args.Input
output = args.Output
verbose = args.Verbose

exclusion_list = ["MWGS220152", "MWGS220167", "MWGS220170", "MWGS220172", "MWGS220173", "MWGS220191", "MWGS220192", "MWGS220267", "MWGS220268", "MWGS220314", "MWGS220332", "MWGS220333", "MWGS220338", "MWGS220349", "MWGS220351", "MWGS220353", "MWGS220357", "MWGS220371", "MWGS220386", "MWGS220400", "MWGS220404", "MWGS220407", "MWGS220408", "MWGS220414", "MWGS220425", "MWGS220428", "MWGS220436", "MWGS220439", "MWGS220450", "MWGS220473", "MWGS220474", "MWGS220483", "MWGS220491", "MWGS220505", "MWGS220520", "MWGS220527", "MWGS220541", "MWGS220545", "MWGS220549", "MWGS220558", "MWGS220570"]


file_list = os.listdir(input)  # Creates a list of all files in directory (specified at input)
# define variable to load the wookbook
count_excl_files = 0

for file in file_list:  # every sample report is iterated
    if file.endswith(".fasta"):
        split1 = file.split(".fasta")
        new_name = split1[0]
        if new_name in exclusion_list:
            count_excl_files += 1
            print(f"{new_name} is excluded {count_excl_files}/41")
        file_path = os.path.join(input + file)  # path of the file (one of the sample reports)
        destination = output + new_name  # adds the new name to the output folder
        if os.path.exists(destination):
            print(f"{new_name} already exists")
        else:
            subprocess.run(f"snippy --ctgs {file_path} --reference fastas/reference_assemblies/st131_assemblies/assembly01.fasta --cpus 10 --ram 8 --outdir results/{destination}", shell=True)

stop = timeit.default_timer()  # ends timer

print('Time: ', stop - start)  # prints script run duration
