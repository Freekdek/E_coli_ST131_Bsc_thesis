"""
Credits:
- Freek de Kreek, 2022

Description:
Create a sammple list in either PopPUNK format or chewBBACA format.

Changes:
@date-of-change, description of added feature/edit

"""

# Import packages here
import argparse
import os
import timeit

# Initialize parser, for all arguments given in command-line
parser = argparse.ArgumentParser()
# Main arguments
parser.add_argument("-i", "--Input",
                    help="Provide input file name of the distance matrix in '.xlsx' format, see help for other input "
                         "options", required=True)
parser.add_argument("-o", "--Output", help="Provide sample list file name", required=True)

# Optional arguments
parser.add_argument("-v", "--Verbose", help="When this argument is given, it will turn on debug, this will print A LOT",
                    action="store_true")
parser.add_argument("-c", "--chewBBACA", help="Will create a sample list for chewBBACA, default is in poppunk format",
                    action="store_true")
parser.add_argument("-s", "--snippy", help="Will create a sample list for snippy, default is in poppunk format",
                    action="store_true")

# Variables
args = parser.parse_args()
main_input = args.Input
main_output = args.Output
verbose = args.Verbose
start = timeit.default_timer()


######################################## Functions #####################################################################

def write_sample_list(input_dir, outfile):
    """
    Put a function description here
    :param input_dir:
    :param outfile:
    :return:
    """
    file_list = os.listdir(input_dir)
    with open(outfile, "w") as output_file:

        for file_name in file_list:
            if file_name.endswith(".fa"):
                sample_name = file_name.split(".fa")[0]
            if file_name.endswith(".fasta"):
                sample_name = file_name.split(".fasta")[0]
            print(f"{sample_name} is being written to sample list")
            if args.snippy:
                output_file.write(f"{sample_name} ")
            elif args.chewBBACA:
                output_file.write(f"results/bakta_results/{sample_name}/{sample_name}.ffn\n")
            else:
                output_file.write(f"{sample_name}\t{input_dir}{file_name}\n")
    return print("Finished creating a sample list...")


######################################## Main ##########################################################################
# Call functions here

write_sample_list(main_input, main_output)


# End of script
stop = timeit.default_timer()  # ends timer
print("Finished succesfully, time taken: ", stop - start)  # prints script run duration
