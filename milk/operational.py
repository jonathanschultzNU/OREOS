"""
operational.py

This module contains utility functions for file handling, logging, and 
saving simulation results for the OREOS optical response simulator.

Functions:
- initialize_logfile: Creates a new log file with a welcome message.
- print_to_log_file: Appends messages to the log file.
- read_input_file: Parses the simulation input file into a dictionary.
- check_dir_exists: Ensures that necessary directories exist or creates them if missing.
- save_simulation_results: Saves the simulation results and metadata to files.

"""

import os
import pickle

# . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . 

def initialize_logfile(logfilename):
    """
    Initializes a new log file with a welcome message.

    Parameters:
        logfilename (str): Path to the log file.
    """
    
    with open(logfilename, "w") as logfile: 
        logfile.write(r"""
####### ######  ####### #######  #####  
#     # #     # #       #     # #     # 
#     # #     # #       #     # #       
#     # ######  #####   #     #  #####  
#     # #   #   #       #     #       # 
#     # #    #  #       #     # #     # 
####### #     # ####### #######  ##### 

Thank you for using OREOS!

""")

    pass

# . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . 

def print_to_log_file(filename, string):
    """
    Appends a message to the log file.

    Parameters:
        filename (str): Path to the log file.
        string (str): Message to be written to the log.
    """
    
    with open(filename, "a") as logfile: 
        logfile.write(f'{string}\n')
        
# . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . 
        
def read_input_file(inputfilename):
    """
    Reads the simulation input file and parses it into a dictionary.

    Parameters:
        inputfilename (str): Path to the input file.

    Returns:
        dict: Parsed input parameters.
    """
    
    p = {}
    f = open(inputfilename)
    data = f.readlines()
    for line in data:
        key, value = line.split("=")
        if key.strip() == "vib_freq" or key.strip() == "displacement" or key.strip() == "vmax":
           entries = value.count(",") + 1 
           key_n = {}
           value_n = value.split(",")
           for i in range(entries):
               key_n[str(i+1)] = value_n[i].strip()   
           p[key.strip()] = key_n
        else:
           p[key.strip()] = value.strip()
    f.close()
    
    return p

# . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . 

def check_dir_exists(path):
    """
    Checks if a directory exists; creates it if missing.

    Parameters:
        path (str): Directory path to check or create.
    """
    
    if os.path.exists(path):
        print(f"'{path}' detected.")
    else:
        try:
            os.makedirs(path)
            print(f"Directory '{path}' created.")
        except OSError as e:
            print(f"Error creating directory '{path}': {e}")

# . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . 

def save_simulation_results(results, time_elapsed, filenames):
    """
    Saves simulation results and metadata to specified files.

    Parameters:
        results (dict): Dictionary containing simulation results.
        time_elapsed (dict): Dictionary with 'CPU' and 'wall' time elapsed.
        filenames (dict): Dictionary containing paths for saving results.

    Outputs:
        - Saves results as a pickle file.
        - Logs timing and metadata information.
    """

    print_to_log_file(filenames['log file'],'Now saving data.')
    with open(filenames['results'], 'wb') as f:
        pickle.dump(results, f)
    
    print_to_log_file(filenames['log file'], f"Data saved to {filenames['results']}")
    
    log_pkl = {'parameters': results['p'],
               'cpu_time_elapsed': time_elapsed['CPU'],
               'wall_time_elapsed': time_elapsed['wall']}
    
    with open(filenames['results log'], 'wb') as g:
        pickle.dump(log_pkl, g)
    
    print_to_log_file(filenames['log file'], f'''Data log saved to {filenames['results log']}
All outputs successfully saved!

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
~~~ OREOS operations complete! Take care! ~~~
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~''')  

# . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . 
