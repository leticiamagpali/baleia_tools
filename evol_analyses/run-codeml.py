####### Run codeml #########
############################

# Author: Letícia Magpali

'''This script runs the program codeml for multiple control files,
each one inside its own folder'''

import subprocess
import sys
import os

# This is the master folder where all folders with control files for each run are located
RUN_FOLDER = sys.argv[1]

SUBPROCESSES = []

for folder in os.listdir(RUN_FOLDER):
    folder_path = os.path.join(RUN_FOLDER, folder)
    if os.path.isdir(folder_path):
        for control_file in os.listdir(folder_path):
            if control_file.endswith(".ctl"):
                control_file_path = os.path.join(folder_path, control_file)
                SUBPROCESSES.append(subprocess.Popen(
                    ["nohup", "codeml", control_file_path], cwd=folder_path))

for process in SUBPROCESSES:
    process.wait()
