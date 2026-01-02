"""
main.py

This script serves as the main execution entry point for the OREOS optical response simulator. 
It initializes directories, reads input parameters, sets defaults, generates the system Hamiltonian, 
performs response function simulations, and processes the results.

Workflow:
1. Reads the simulation input file and sets default values where needed.
2. Generates the system Hamiltonian.
3. Simulates third-order nonlinear optical response functions.
4. Applies post-processing (lineshapes, normalization, and windowing).
5. Saves the processed results.

Usage:
    python main.py <runname>

Arguments:
    runname (str): Name of the simulation run, used for input and output file naming.
"""

import sys
import os
import time
from milk.operational import initialize_logfile, read_input_file, check_dir_exists, save_simulation_results
from milk.simulations import process_time_domain_response, response_simulation
from milk.system_model import get_Hamiltonian
from milk.defaults import check_inputs

# Retrieve run name from command-line arguments
runname = sys.argv[1]

paths = {'working': os.getcwd(),
        'outputs': os.path.join(os.getcwd(), 'outputs'),
        'inputs': os.path.join(os.getcwd(), 'inputs'),
        'log file': f"{os.getcwd()}/outputs/{runname}.log",
        'input file': f"{os.getcwd()}/inputs/{runname}.txt",
        'results': f"{os.getcwd()}/outputs/{runname}.pkl",
        'results log': f"{os.getcwd()}/outputs/{runname}LOG.pkl"
        }

for directory in [paths['outputs'], paths['inputs']]:
    check_dir_exists(directory)

initialize_logfile(paths['log file'])
print(f"See {paths['log file']} for simulation details.", flush=True)

# Read and validate input parameters
p = read_input_file(paths['input file'])
p = check_inputs(p, paths['log file'])

start_cpu_time = time.process_time()
start_wall_time = time.time()

# Create system Hamiltonian
H = get_Hamiltonian(p, paths['log file'])

# Perform response function simulations
R, Rt1t3 = response_simulation(p, H, paths['log file'])

# Post-processing (normalization, lineshapes, windowing)
results = process_time_domain_response(p, R, Rt1t3, H, paths['log file'])

# Calculate the elapsed times
time_elapsed = {'CPU': time.process_time() - start_cpu_time,
                'wall': time.time() - start_wall_time}

# Save results
save_simulation_results(results, time_elapsed, paths)
