# Created by Freek de Kreek, 3-11-2022
# For each file:
# run bakta command

import os
import argparse
import subprocess
import timeit

start = timeit.default_timer()  # starts timer

# initialize parser
parser = argparse.ArgumentParser()
# adding optional argument
parser.add_argument("-i", "--Input", help="Show Input sample reports directory")
parser.add_argument("-o", "--Output", help="Show Output directory")
parser.add_argument("-v", "--Verbose", help="Show more information")
# read arguments from command line
args = parser.parse_args()
input = args.Input
database = args.Database
output = args.Output
verbose = args.Verbose

if verbose == "True":
    debug = True
else:
    debug = False

if input and database and output:
    file_list = os.listdir(input)  # Creates a list of all files in directory (specified at input)
    # define variable to load the wookbook

    for file in file_list:  # every sample report is iterated
        if file.endswith(".fa"):
            file_path = os.path.join(input + file)  # path of the file (one of the sample reports)
            if debug:
                print("Will now run bakta for " + file)
            split1 = file.split(".fa")
            new_name = split1[0]
            destination = output + new_name  # adds the new name to the output folder
            if debug:
                print(new_name)
            if os.path.exists(destination):
                print(f"{new_name} already exists")
            else:
                subprocess.run("bakta reference_assemblies/st131_assemblies/" + file + "--db bakta_db/db/ --genus "
                                                                                       "Escherichia --species coli "
                                                                                       "--complete --threads 20 "
                                                                                       "--output "
                                                                                       "results/bakta_results/" + new_name, shell=True)

stop = timeit.default_timer()  # ends timer

print('Time: ', stop - start)  # prints script run duration
