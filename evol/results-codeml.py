######### Extract results from codeml runs #########
####################################################

# Author: Letícia Magpali

# Loop through run folder for that model
# For each folder, go through files 
# if file has "out" in name, read file
# locate lnL, np, omega

import sys
import os
import glob
import csv
from scipy.stats import chi2
import re

# Define the path to the master folder
RUN_FOLDER_ALT = sys.argv[1]
RUN_FOLDER_NULL = sys.argv[2]
OUTPUT_FILE = sys.argv[3]

def get_codeml_results(run_folder):
    # Creating a dictionary to store all results in a master folder
    folder_results = {}
    # Loop through each subfolder in the master folder
    for subfolder in os.listdir(run_folder): 
        subfolder_path = os.path.join(run_folder, subfolder)
        
        # Check if the path is a directory (to skip any files in the master folder itself)
        if os.path.isdir(subfolder_path):
            # Getting gene name from run folder
            gene_name = subfolder.split("_")[1]
            model_name = subfolder.split("_")[0]
            # Use glob to find files with "out" in their names in the current subfolder
            outfiles = glob.glob(os.path.join(subfolder_path, '*out*'))

            # Loop through the found files and read their contents
            if outfiles:  # Proceed only if a matching file is found
                outfile_path = outfiles[0]  # Get the first (and only) match
                with open(outfile_path, "r") as outfile:
                    for line in outfile:
                        # Check if "lnL" is in the line 
                        # splitting the line and looking for the lnL and np values
                        if "lnL" in line:
                            if "lnL" in line:
                                # Regular expression to extract lnL value and np value
                                match = re.search(r"lnL\(ntime:\s*\d+\s+np:\s*(\d+)\):\s*(-?\d+\.\d+)", line)
                                if match:
                                    np_value = match.group(1)  # Captures the value of np
                                    lnL_value = match.group(2)  # Captures the lnL value
                                    print(f"np: {np_value}, lnL: {lnL_value}")
                                else:
                                    print("Could not parse the lnL line.")
                
                # Store results in dictionary by (gene_name, model_name) key
                folder_results[(gene_name)] = {
                    "Gene": gene_name,
                    "Model": model_name,
                    "lnL": float(lnL_value),
                    "np": int(np_value),
                }
    return folder_results


# Gather results from both folders
results_folder1 = get_codeml_results(RUN_FOLDER_ALT)
results_folder2 = get_codeml_results(RUN_FOLDER_NULL)
                

# Merge results for calculations and write to CSV
with open(OUTPUT_FILE, mode='w', newline='') as csvfile:
    fieldnames = ["Gene", "Model", "lnL", "LRT", "np", "df", "p-value"]
    writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
    writer.writeheader()
    
    for key in results_folder1.keys():
        if key in results_folder2:
            # Get corresponding results from both folders
            data1 = results_folder1[key]
            data2 = results_folder2[key]

            # Calculate LRT and df
            lrt = 2 * (data1["lnL"] - data2["lnL"])
            df = data1["np"] - data2["np"]

            # Calculate the p-value using the Chi-squared distribution
            p_value = chi2.sf(lrt, df)

            # Write data for RUN_FOLDER_1 with calculated LRT and df
            writer.writerow({
                "Gene": data1["Gene"],
                "Model": data1["Model"],
                "lnL": data1["lnL"],
                "LRT": lrt,
                "np": data1["np"],
                "df": df,
                "p-value": p_value
            })

            # Write data for RUN_FOLDER_2 with calculated LRT and df
            writer.writerow({
                "Gene": data2["Gene"],
                "Model": data2["Model"],
                "lnL": data2["lnL"],
                "LRT": "N/A", # Indicating no LRT calculation for this row
                "np": data2["np"],
                "df": "N/A", # Indicating no df calculation for this row
                "p-value": "N/A" # Indicating no p_value calculation for this row
            })