######### Extract results from codeml runs #########
####################################################

# Author: Letícia Magpali

""" This script extract results from codeml output files.

INPUTS:
It requires as input two folders: one containing the results for the alternative model, and another for the null model.
The folders should be organized into subfolders, one for each gene, with the codeml output files inside.
Before running the script, the subfolders should be named in the format: GENE_MODEL (e.g., gene1_branch, gene2_site). 
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
import warnings
from scipy.stats import chi2
from itertools import dropwhile, islice

# Precompile regex for speed

LNL = re.compile(r"lnL\(ntime:\s*\d+\s+np:\s*(\d+)\):\s*(-?\d+\.\d+)", re.IGNORECASE)
PSS = re.compile(r"^\s*(\d+)\s+([A-Za-z-])\s+([01]\.\d+\*{0,2})", re.IGNORECASE | re.VERBOSE)
SITE_CLASS = re.compile(r"[\d.]+")
OMEGA = re.compile(r"(\d+\.\d+)")
SITE_CLASS_LINE = re.compile(r"site class\s+", re.IGNORECASE)
MODEL_TYPES = {"branch", "site", "branchsite"}



# ---------------------- Functions ---------------------- #


# Functions to extract specific information from codeml output files 
# ------------------------------------------------------------------ # 

# The next functions extract specific pieces of information 
# (like lnL, np, omeg, positive sites) from the codeml output files
# by searching for regular expressions and capturing groups.

def extract_lnl(line: str) -> tuple[str | None, str | None]:
    """ Extracts the lnl value and np value from a given line of text. """
    match = LNL.search(line)
    if match:
        return match.group(2), match.group(1)  
    return None, None

def extract_positive_sites(beb_line):
    """ Extracts positive sites from a given line of text in BEB section. """
    
    match = PSS.search(beb_line)
    if match:
        site, aa, prob = match.groups()
        return [f'{site}({aa},{prob})']
    else:
        return []

def calculate_lrt_pvalue(lnl_alt, lnl_null, np_alt, np_null, gene_name):
    """Calculate an LRT p-value after validating the model comparison."""
    lrt = 2 * (lnl_alt - lnl_null)
    df = np_alt - np_null

    if df <= 0:
        raise ValueError(f"Degrees of freedom must be positive for gene {gene_name}; got {df}")
    if lrt < 0:
        warnings.warn(
            f"LRT is negative for gene {gene_name}; got {lrt}",
            RuntimeWarning,
            stacklevel=2,
        )

    return lrt, df, chi2.sf(lrt, df)

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
    if model_type not in MODEL_TYPES:
        raise ValueError(f"Unsupported model type: {model_type}")

    # 1) Getting paths from output files ########
    #############################################

    # Create a dictionary to store all results from the run folder
    folder_results = {}
    
    # Check MODEL_TYPE once outside the loop
    # assign boolean flags for model types
    is_branch = model_type == "branch"
    is_branchsite = model_type == "branchsite"
    is_site = model_type == "site"

    # Loops through each gene folder in the run folder, and creates a list of their paths
    with os.scandir(run_folder) as run_entries:
        genefolder_paths = [entry.path for entry in run_entries if entry.is_dir()]


    # Extraction loop 
    #----------------
    
    # --> The following loop goes inside each gene folder, locates the output file,
    # and extracts the relevant information based on the model type.
    
    for path in genefolder_paths:

        # Extract gene name and model name from folder name
        gene_name = os.path.basename(path).split("_")[0]
        model_name = os.path.basename(path).split("_")[1] 

        # Locate the codeml output file in the gene folder
        with os.scandir(path) as genefolder:
            output_files = [
                entry.path
                for entry in genefolder
                if entry.is_file()
                and entry.name.startswith("out_")
                and entry.name.endswith(".txt")
            ]
        if not output_files:
            raise FileNotFoundError(f'No codeml output file found in folder: {path}')
        if len(output_files) > 1:
            raise ValueError(f"Multiple codeml output files found in folder: {path}")
        outfile_path = output_files[0]
        
        
    # 2) Extracting parameters, branch models ######
    ################################################

        if is_branch:

            '''This section extracts the parameters for branch models. It reads the output file and looks for the 
            ln, omega and np values.'''
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
                    if lnl_value is not None and np_value is not None and (omega_alt or omega_null):
                        break

            # If values were found, store them in the results dictionary
            # This is a dictionary of dictionaries, with gene_name as key, 
            # and another dictionary (made of paramater name and value) as value
            if lnl_value is not None and np_value is not None and (omega_alt or omega_null):
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

    
    # 3) Extracting parameters, site models #####
    #############################################
        
        if is_site:
        
            '''This section extracts results from site model output files.
            It looks for the model name, lnL value, np value, and positive selection sites (if applicable) in the output file. 
            The results are stored in a nested dictionary with the following structure: 
            {gene_name: [{model_name: ..., lnl: ..., np: ..., pss: [...]}]}.'''
            
            with open(outfile_path, "r") as out:
                for line in out:
                    current_model: dict[str, object] = {}
                    model_name = None
                    lnl_value = None
                    np_value = None
                    positive_sites = []

                    header = dropwhile(lambda l: "nssites model" not in l.lower(), out)
                    for line in header:
                        # Check for model name in the line
                        # If found, begin a new current_model dictionary to store results for that model
                        if "nssites model" in line.lower():
                            model_name = line.strip().split(":")[0]
                            model_name = model_name.lower()
                            current_model = {"model": model_name}
                            break

                    if model_name is None:
                        continue

                    # Fast-forward to lnl section
                    # Loops through the output file until it finds the line containing "lnl" 
                    # and extracts the values
                    for line in out:
                        if "lnl" in line.lower():
                            lnl_value, np_value = extract_lnl(line)
                            if lnl_value is not None and np_value is not None:
                                current_model["lnl"] = float(lnl_value)
                                current_model["np"] = int(np_value)
                            break

                    # This step is only for models that contain BEB sites (i.e., model 2 and model 8)
                    if "model 2" in model_name or "model 8" in model_name:
                        beb_lines = dropwhile(lambda l: "bayes empirical bayes (beb) analysis" not in l.lower(), out)
                        for line in beb_lines:
                            if PSS.search(line):
                                positive_sites.extend(extract_positive_sites(line))
                            if "the grid" in line.lower():
                                break

                    # Add the list of positive sites to the current_model dictionary
                        current_model["pss"] = positive_sites

                    if lnl_value is not None and np_value is not None:
                        folder_results.setdefault(gene_name,[]).append(current_model)   


    # 4) Extracting parameters, branch-site models #####
    ####################################################
        
        if is_branchsite:
            '''This section extracts parameters for branch-site models.'''
            lnl_value = None
            np_value = None
            # Intilize an empty list to store positively selected sites (PSS) if found
            positive_sites = []
            # Initialize an empty dictionary to store the site class table data, 
            # this table is found only in branch-site model output files
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
                        lines_to_parse = [line, *islice(out, 3)]
                        if len(lines_to_parse) != 4:
                            raise ValueError(
                                f"Incomplete site class table for gene {gene_name} "
                                f"in {outfile_path}"
                            )
                        table_data = parse_siteclass_table(lines_to_parse)
                    
                    # Check for BEB section to extract positive sites
                    elif "bayes empirical bayes (beb) analysis" in lower_case:
                        for beb_line in out:
                            if PSS.search(beb_line):
                                positive_sites.extend(extract_positive_sites(beb_line))

            if lnl_value is None or np_value is None:
                raise ValueError(f'Could not find lnL or np values for gene {gene_name}')
                
            
            # Extract omega values and proportion from site class table
            if lnl_value is not None and np_value is not None and table_data:
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
                if lnl_value is not None and np_value is not None:
                    folder_results[gene_name] = {
                    "gene": gene_name,
                    "model": model_name,
                    "lnl": float(lnl_value),
                    "np": int(np_value)
                }

    # 5) Callout of results for each run folder #####
    #################################################

    print(f'Extracted results from {run_folder}:')

    # Print all lnL and np values
    if is_branch or is_branchsite:
        for gene_name, data in sorted(folder_results.items()):
            pss = " ".join(data["pss"]) if "pss" in data else ""
            pss_str = f', PSS: {pss}' if pss else ""
            print(f'{gene_name}, lnl: {data["lnl"]}, np: {data["np"]}{pss_str}')

    return folder_results


# Function to create a table from extracted codeml results
#----------------------------------------------------------- #
                
def write_codeml_results(output_file, model_type, *,
                         results_alt = None, results_null = None, results_site = None):
    """ Merges results from two codeml result dictionaries and writes them to a CSV file."""
    if model_type not in MODEL_TYPES:
        raise ValueError(f"Unsupported model type: {model_type}")

    # Check MODEL_TYPE once outside the loop
    # assign boolean flags for model types
    is_branch = model_type == "branch"
    is_branchsite = model_type == "branchsite"
    is_site = model_type == "site"

    if (is_branch or is_branchsite) and (results_alt is None or results_null is None):
        raise ValueError(f"results for {model_type} not generated")
    if is_site and results_site is None:
        raise ValueError("Results for site model not generated")
    
    # Merge results for calculations and write to CSV
    with open(output_file, mode='w', newline='') as csvfile:
    

        # Defining the fieldnames for the CSV file based on the model type
        if  is_branch:
            fieldnames = ["gene", "model", "lnl_alt", "lnl_null", "lrt", "np_alt", "np_null", "df", "pvalue", "omega_alt", "omega_null"]
        if is_branchsite:
            fieldnames = ["gene", "model", "lnl_alt", "lnl_null", "lrt", "np_alt", "np_null", "df", "pvalue", "pss",
                            "fg_omega_2a (prop)", "fg_omega_2b (prop)", "fg_omega_max", "has_pos_sel", "branch_site_omegas"]
        if is_site:
            fieldnames = ["gene", "model", "lnl", "lrt", "np", "df", "pvalue", "pss"]
        
        # Create a CSV DictWriter object and write the header row
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        writer.writeheader()

        if is_branch or is_branchsite:
            if results_alt is None or results_null is None:
                raise ValueError(f"results for {model_type} not generated")
            else:
                for key in sorted(results_alt.keys(), key=str.lower):
                    if key not in results_null:
                        continue

                    # Get corresponding results from both folders
                    data_alt = results_alt[key]
                    data_null = results_null[key]

                    # Calculate LRT and df
                    lrt, df, p_value = calculate_lrt_pvalue(
                        data_alt["lnl"],
                        data_null["lnl"],
                        data_alt["np"],
                        data_null["np"],
                        key)

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
                        # Make change here to add a third condition for positive selection: p-value < 0.05, even if omega is not >1 or proportion is not >0
                        row["has_pos_sel"] = "Yes" if (data_alt["foreground_w2a"] > 1 or data_alt["foreground_w2b"] > 1) and (data_alt["proportion_2a"] > 0 or data_alt["proportion_2b"] > 0) else "No"
                        row["branch_site_omegas"] = json.dumps(data_alt["branchsite_omegas"])

                    writer.writerow(row)

        if is_site:
            if results_site is None:
                raise ValueError("Results for site model not generated")
            else:
                for gene_key in sorted(results_site, key=str.lower):
                    models = results_site[gene_key]
                    if len(models) != 4:
                        raise ValueError(
                            f"Expected four site-model results for gene {gene_key}; "
                            f"found {len(models)}"
                        )
                    results_model1 = models[0]
                    results_model2 = models[1]
                    results_model7 = models[2]
                    results_model8 = models[3]
                    lrt1, df1, p_value1 = calculate_lrt_pvalue(
                        results_model2["lnl"],
                        results_model1["lnl"],
                        results_model2["np"],
                        results_model1["np"],
                        gene_key)
                    
                    lrt2, df2, p_value2 = calculate_lrt_pvalue(
                        results_model8["lnl"],
                        results_model7["lnl"],
                        results_model8["np"],
                        results_model7["np"],
                        gene_key)
                    
                    row1 = {
                        "gene": gene_key,
                        "model": results_model2["model"],
                        "lnl": results_model2["lnl"],
                        "lrt": lrt1,
                        "np": results_model2["np"],
                        "df": df1,
                        "pvalue": p_value1,
                        "pss": " ".join(results_model2["pss"])
                    }

                    row2 = {
                        "gene": gene_key,
                        "model": results_model8["model"],
                        "lnl": results_model8["lnl"],
                        "lrt": lrt2,
                        "np": results_model8["np"],
                        "df": df2,
                        "pvalue": p_value2,
                        "pss": " ".join(results_model8["pss"])}

                    writer.writerow(row1)
                    writer.writerow(row2)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Extract codeml results (lnL, np, omega, PSS) from two models and calculate LRT.")
    parser.add_argument("--alt", help="Path to the first run folder (alternative model).")
    parser.add_argument("--null", help="Path to the second run folder (null model).")
    parser.add_argument("--site", help="Path to the site model run folder (if applicable).")
    parser.add_argument("--out", required=True, help="Path to the output CSV file.")
    parser.add_argument(
        "--m",
        required=True,
        choices=("branch", "site", "branchsite"),
        help="Type of model being analyzed.",
    )
    args = parser.parse_args()

    if args.m == "site":
        if not args.site:
            parser.error("--site is required when --m is 'site'.")
        results_site = get_codeml_results(args.site, args.m)
        write_codeml_results(args.out, args.m, results_site=results_site)
    else:
        if not args.alt or not args.null:
            parser.error("--alt and --null are required for branch and branchsite models.")
        results_alt = get_codeml_results(args.alt, args.m)
        results_null = get_codeml_results(args.null, args.m)
        write_codeml_results(
            args.out,
            args.m,
            results_alt=results_alt,
            results_null=results_null,
        )