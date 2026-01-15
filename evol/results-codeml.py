######### Extract results from codeml runs #########
####################################################

# Author: Letícia Magpali

# Loop through run folder for that model
# For each folder, go through files 
# if file has "out" in name, read file
# locate lnL, np, omega

import argparse
import os
import glob
import csv
from scipy.stats import chi2
import re

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Extract codeml results from two run folders and calculate LRT.")
    parser.add_argument("--alt", help="Path to the first run folder (alternative model).")
    parser.add_argument("--null", help="Path to the second run folder (null model).")
    parser.add_argument("--out", help="Path to the output CSV file.")
    parser.add_argument("--m", help="Type of model being analyzed (branch, site or branchsite).")
    args = parser.parse_args()

    # Define the path to the master folder
    RUN_FOLDER_ALT = args.alt
    RUN_FOLDER_NULL = args.null
    OUTPUT_FILE = args.out
    MODEL_TYPE = args.m

    def extract_lnl(line):
        # Regular expression to extract lnL value and np value
        match = re.search(r"lnL\(ntime:\s*\d+\s+np:\s*(\d+)\):\s*(-?\d+\.\d+)", line)
        if match:
            np_value = match.group(1)  # Captures the value of np
            lnL_value = match.group(2)  # Captures the lnL value
            print(f"np: {np_value}, lnL: {lnL_value}")
        else:
            print("Could not parse the lnL line.")
        return lnL_value, np_value

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
                    if MODEL_TYPE == "branch":
                        with open(outfile_path, "r") as outfile:
                            for line in outfile:
                                # Check if "lnL" is in the line and extract values
                                if "lnL" in line:
                                    lnL_value, np_value = extract_lnl(line)

                        # Store results in dictionary by (gene_name, model_name) key
                        folder_results[(gene_name)] = {
                            "Gene": gene_name,
                            "Model": model_name,
                            "lnL": float(lnL_value),
                            "np": int(np_value)
                        }

                    if MODEL_TYPE == "branchsite" or MODEL_TYPE == "site":
                        with open(outfile_path, "r") as outfile:
                            positive_sites = []
                            for line in outfile:
                                # Check if "lnL" is in the line and extract values
                                if "lnL" in line:
                                    lnL_value, np_value = extract_lnl(line)
                                # Check for positive sites section
                                if "Bayes Empirical Bayes (BEB)" in line:
                                    positive_sites = []
                                    # Extract positive sites information
                                    for next_line in outfile:
                                        ps_match = re.match(r"^\s*(\d+)\s+([A-Za-z])\s+([01]\.\d+)\s*$", next_line)
                                        if ps_match:
                                            ps_string = f"{ps_match.group(1)}({ps_match.group(2)},{ps_match.group(3)})"
                                            positive_sites.append(ps_string)
                                        else:
                                            break  # Exit loop if the line does not match the expected format

                    
                        # Store results in dictionary by (gene_name, model_name) key
                        folder_results[(gene_name)] = {
                            "Gene": gene_name,
                            "Model": model_name,
                            "lnL": float(lnL_value),
                            "np": int(np_value),
                            "Positive_Sites": positive_sites
                        }

        return folder_results


    # Gather results from both folders
    results_folder1 = get_codeml_results(RUN_FOLDER_ALT)
    results_folder2 = get_codeml_results(RUN_FOLDER_NULL)
                    

    # Merge results for calculations and write to CSV
    with open(OUTPUT_FILE, mode='w', newline='') as csvfile:
        fieldnames = ["Gene", "Model", "lnL", "LRT", "np", "df", "p-value"]
        if MODEL_TYPE in ["branchsite", "site"]:
            fieldnames.append("PSS")
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

                # Build row for RUN_FOLDER_1, conditionally including PSS
                row1 = {
                    "Gene": data1["Gene"],
                    "Model": data1["Model"],
                    "lnL": data1["lnL"],
                    "LRT": lrt,
                    "np": data1["np"],
                    "df": df,
                    "p-value": p_value,
                }
                if "Positive_Sites" in data1:
                    row1["PSS"] = " ".join(data1["Positive_Sites"])
                writer.writerow(row1)

                # Build row for RUN_FOLDER_2, with PSS set to "NA" (since null model typically has no positive sites)
                row2 = {
                    "Gene": data2["Gene"],
                    "Model": data2["Model"],
                    "lnL": data2["lnL"],
                    "LRT": "NA",  # Indicating no LRT calculation for this row
                    "np": data2["np"],
                    "df": "NA",  # Indicating no df calculation for this row
                    "p-value": "NA",  # Indicating no p_value calculation for this row
                }
                writer.writerow(row2)