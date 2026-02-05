######### Extract results from codeml runs #########
####################################################

# Author: Letícia Magpali

""" This script extract results from codeml output files.

INPUTS:
It requires as input two folders: one containing the results for the alternative model, and another for the null model.
The folders should be organized into subfolders, one for each gene, with the codeml output files inside.
Before running the script, the subfolders should be named in the format: MODEL_GENE (e.g., branch_gene1, site_gene2). 
The model name and gene name will be extracted from the folder name, and used in the output.
Additionally, the user must specify the type of model being analyzed: branch, site, or branchsite, 
and a name/path for the output CSV file.

RESULTS:
The script generates a CSV file containing the extracted results, including lnL, np, LRT, 
degrees of freedom (df), p-value, and positive selection sites (PSS) if applicable."""

import os 
import re
import json 
import csv
import argparse
from scipy.stats import chi2
from itertools import dropwhile

# Precompile regex for speed

LNL = re.compile(r"lnL\(ntime:\s*\d+\s+np:\s*(\d+)\):\s*(-?\d+\.\d+)", re.IGNORECASE)
PSS = re.compile(r"^\s*(\d+)\s+([A-Za-z])\s+([01]\.\d+\**)\s*$", re.IGNORECASE | re.VERBOSE)
SITE_CLASS = re.compile(r"[\d.]+")
OMEGA = re.compile(r"(\d+\.\d+)")
SITE_CLASS_LINE = re.compile(r"site class\s+", re.IGNORECASE)



# ---------------------- Functions ---------------------- #


# Functions to extract specific information from codeml output files 
# ------------------------------------------------------------------ # 

# The next functions extract specific pieces of information 
# (like lnL, np, omeg, positive sites) from the codeml output files
# by searching for regular expressions and capturing groups.

def extract_lnl(line):
    """ Extracts the lnl value and np value from a given line of text. """
    match = LNL.search(line)
    if match:
        # Returns lnl_value (match group 2), np_value (match group 1)
        return match.group(2), match.group(1)  
    else:
        return None, None

def extract_positive_sites(beb_line):
    """ Extracts positive sites from a given line of text in BEB section. """
    
    matches = PSS.findall(beb_line)
    if matches:
        # Returns a list of formatted strings: site(aa,prob)
        return [f'{site}({aa},{prob})' for site, aa, prob in matches]
    else:
        return []

def parse_siteclass_table(lines):
    """ Parses the site class table from codeml branch-site output file. """
    
    table_data = {}
    proportion = background_w = foreground_w = None
    for line in lines:
        if "proportion" in line:
            proportion = SITE_CLASS.findall(line)
        if "background w" in line:
            background_w = SITE_CLASS.findall(line)
        if "foreground w" in line:
            foreground_w = SITE_CLASS.findall(line)

    if proportion and background_w and foreground_w:
        table_data["proportion"] = proportion
        table_data["background_w"] = background_w
        table_data["foreground_w"] = foreground_w   
    
    return table_data



# Main function to extract results from codeml output files 
#----------------------------------------------------------- #

def get_codeml_results(run_folder, model_type):
    """ Extracts results from codeml output files located in subfolders of a specified run folder. """

    # 1) Getting paths from output files ########
    #############################################

    # Create a dictionary to store all results from the run folder
    folder_results = {}
    
    # Check MODEL_TYPE once outside the loop
    # assign boolean flags for model types
    is_branch_model = model_type == "branch"
    is_branchsite = model_type == "branchsite"

    # Loops through each gene folder in the run folder, and creates a list of their paths
    genefolder_paths = [os.path.join(run_folder, d) for d in os.listdir(run_folder) if os.path.isdir(os.path.join(run_folder, d))]


    # Extraction loop 
    #----------------
    
    # --> The following loop goes inside each gene folder, locates the output file,
    # and extracts the relevant information based on the model type.
    
    for path in genefolder_paths:

        # Extract gene name and model name from folder name
        gene_name = os.path.basename(path).split("_")[1]
        model_name = os.path.basename(path).split("_")[0]

        # Locate the codeml output file in the gene folder
        outfile_path = None
        with os.scandir(path) as genefolder:
            for entry in genefolder:
                if entry.is_file() and 'out' in entry.name:
                    outfile_path = entry.path
                    break
        if outfile_path is None:
            raise FileNotFoundError(f'No codeml output file found in folder: {path}')
        
        
    # 2) Extracting parameters, branch models ######
    ################################################

        if is_branch_model:
            # Initialize variables to store extracted values
            lnl_value = np_value = None
            omega_alt = omega_null = None

            # Reading the output file with dropwhile to skip lines until "lnl" is found
            with open(outfile_path, "r") as out:
                lines_to_read = dropwhile(lambda l: "lnl" not in l.lower(), out)
                for line in lines_to_read:
                    lower_case = line.lower()

                    # Check if "lnl" is in the line and extract values
                    if "lnl" in lower_case:
                        lnl_value, np_value = extract_lnl(line)
                    
                    # Extract omega values from alternative model output
                    elif "w (dn/ds) for branches:" in lower_case:
                        omega_alt = OMEGA.findall(line)

                    # Extract omega values from null model output
                    elif "omega (dn/ds) =" in lower_case:
                        omega_null = OMEGA.findall(line)
                    
                    # Break the loop if all values have been found
                    if lnl_value and np_value and (omega_alt or omega_null):
                        break

            # If values were found, store them in the results dictionary
            # This is a dictionary of dictionaries, with gene_name as key, 
            # and another dictionary (made of paramater name and value) as value
            if lnl_value and np_value and (omega_alt or omega_null):
                    folder_results[gene_name] = {
                        "gene": gene_name,
                        "model": model_name,
                        "lnl": float(lnl_value),
                        "np": int(np_value),
                        "omega_alt": omega_alt,
                        "omega_null": omega_null
                    }
            
            # Raise errors if any value is missing
            else:
                if lnl_value is None:
                    raise ValueError(f'Could not find lnL value for gene {gene_name}')
                if np_value is None:
                    raise ValueError(f'Could not find np value for gene {gene_name}') 
                if not (omega_alt or omega_null):
                    raise ValueError(f'No omega values found for gene {gene_name}')  

    
    # 3) Extracting parameters, branch-site models #####
    ####################################################
        
        if is_branchsite:
            lnl_value = None
            np_value = None
            positive_sites = []
            table_data = {}
            with open(outfile_path, "r") as out:
                lines_to_read = dropwhile(lambda l: "lnl" not in l.lower(), out)
                for line in lines_to_read:
                    lower_case = line.lower()
                    
                    # Check if "lnl" is in the line and extract values
                    if "lnl" in lower_case:
                        lnl_value, np_value = extract_lnl(line)
                    
                    # Checks for the first line of the site class table
                    # and loops through the next 3 lines to extract the table data
                    elif SITE_CLASS_LINE.search(lower_case):
                        lines_to_parse = [line]
                        for _ in range(3): 
                            lines_to_parse.append(next(out))
                        table_data = parse_siteclass_table(lines_to_parse)
                    
                    # Check for BEB section to extract positive sites
                    elif "bayes empirical bayes (beb) analysis" in lower_case:
                        for beb_line in out:
                            positive_sites.extend(extract_positive_sites(beb_line))              

            if lnl_value is None or np_value is None:
                raise ValueError(f'Could not find lnL or np values for gene {gene_name}')
                
            
            # Extract omega values and proportion from site class table
            if lnl_value and np_value and table_data:
                prop0, prop1, prop2a, prop2b = map(float, table_data["proportion"])
                bg_w0, bg_w1, bg_w2a, bg_w2b = map(float, table_data["background_w"])
                fg_w0, fg_w1, fg_w2a, fg_w2b = map(float, table_data["foreground_w"])
                
                # Store results in dictionary by (gene_name) key
                folder_results[gene_name] = {
                    "gene": gene_name,
                    "model": model_name,
                    "lnl": float(lnl_value),
                    "np": int(np_value),
                    "pss": positive_sites,
                    "foreground_w2a": fg_w2a,
                    "foreground_w2b": fg_w2b,
                    "proportion_2a": prop2a,
                    "proportion_2b": prop2b,
                    "branchsite_omegas": {"0": {"proportion": prop0, "background_w": bg_w0, "foreground_w": fg_w0},
                                        "1": {"proportion": prop1, "background_w": bg_w1, "foreground_w": fg_w1},
                                        "2a": {"proportion": prop2a, "background_w": bg_w2a, "foreground_w": fg_w2a},
                                        "2b": {"proportion": prop2b, "background_w": bg_w2b, "foreground_w": fg_w2b}
                                        }
                }
            else:
                if lnl_value and np_value:
                    folder_results[gene_name] = {
                    "gene": gene_name,
                    "model": model_name,
                    "lnl": float(lnl_value),
                    "np": int(np_value)
                }
    print(f'Extracted results from {run_folder}:')

    # Print all lnL and np values
    for gene_name, data in sorted(folder_results.items()):
        pss = " ".join(data["pss"]) if "pss" in data else ""
        pss_str = f', PSS: {pss}' if pss else ""
        print(f'{gene_name}, lnl: {data["lnl"]}, np: {data["np"]}{pss_str}')

    return folder_results

                
def write_codeml_results(results_alt, results_null, output_file, model_type):
    """ Merges results from two codeml result dictionaries and writes them to a CSV file."""
    
    # Merge results for calculations and write to CSV
    with open(output_file, mode='w', newline='') as csvfile:
    
        if  model_type == "branch":
            fieldnames = ["gene", "model", "lnl_alt", "lnl_null", "lrt", "np_alt", "np_null", "df", "pvalue", "omega_alt", "omega_null"]
        if model_type in ["branchsite", "site"]:
            fieldnames = ["gene", "model", "lnl_alt", "lnl_null", "lrt", "np_alt", "np_null", "df", "pvalue", "pss",
                            "fg_omega_2a (prop)", "fg_omega_2b (prop)", "fg_omega_max", "has_pos_sel", "branch_site_omegas"]
        
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        writer.writeheader()
        
        for key in sorted(results_alt.keys(), key=str.lower):
            if key in results_null:
                # Get corresponding results from both folders
                data_alt = results_alt[key]
                data_null = results_null[key]
                
                # Calculate LRT and df
                lrt = 2 * (data_alt["lnl"] - data_null["lnl"])
                df = data_alt["np"] - data_null["np"]

                # Calculate the p-value using the Chi-squared distribution
                p_value = chi2.sf(lrt, df)

                # Build row to write to CSV
                row = {
                    "gene": data_alt["gene"],
                    "model": data_alt["model"],
                    "lnl_alt": data_alt["lnl"],
                    "lnl_null": data_null["lnl"],
                    "lrt": lrt,
                    "np_alt": data_alt["np"],
                    "np_null": data_null["np"],
                    "df": df,
                    "pvalue": p_value,
                }
                if "omega_alt" in data_alt:
                    row["omega_alt"] = " ".join(data_alt["omega_alt"])
                if "omega_null" in data_null:
                    row["omega_null"] = " ".join(data_null["omega_null"])
                if "pss" in data_alt:
                    row["pss"] = " ".join(data_alt["pss"])
                if "branchsite_omegas" in data_alt:
                    row["fg_omega_2a (prop)"] = f'{data_alt["foreground_w2a"]} ({data_alt["proportion_2a"]})' 
                    row["fg_omega_2b (prop)"] = f'{data_alt["foreground_w2b"]} ({data_alt["proportion_2b"]})'
                    row["fg_omega_max"] = max(data_alt["foreground_w2a"], data_alt["foreground_w2b"])
                    row["has_pos_sel"] = "Yes" if (data_alt["foreground_w2a"] > 1 or data_alt["foreground_w2b"] > 1) and (data_alt["proportion_2a"] > 0 or data_alt["proportion_2b"] > 0) else "No"
                    row["branch_site_omegas"] = json.dumps(data_alt["branchsite_omegas"])
                
                writer.writerow(row)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Extract codeml results (lnL, np, omega, PSS) from two models and calculate LRT.")
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

    results_alt = get_codeml_results(RUN_FOLDER_ALT, MODEL_TYPE)
    results_null = get_codeml_results(RUN_FOLDER_NULL, MODEL_TYPE)

    write_codeml_results(results_alt, results_null, OUTPUT_FILE, MODEL_TYPE)