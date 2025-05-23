"""
defaults.py

This module defines default parameters and input validation functions for the OREOS optical 
response simulator. It ensures that required parameters are set and assigns default values 
for missing ones.

Functions:
- check_inputs: Validates user input parameters and assigns default values where necessary.
- get_defaults: Returns a dictionary of default system parameters for different models.
- convert_input_types: Converts input parameter values to appropriate data types.
"""

from milk.operational import print_to_log_file

# . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . 

def check_inputs(p, logfilename):
    """
    Validates input parameters and ensures required keys are present.
    Modifies `p` in-place by setting default values for missing parameters.

    Parameters:
        p (dict): Dictionary containing user-defined parameters.
        logfilename (str): Path to the log file for tracking missing parameters.

    Returns:
        dict: Updated dictionary with missing parameters assigned default values.
    """
        
    print_to_log_file(logfilename, 
                      "Now checking for missing input parameters.\n"
                      "(If a parameter is missing, the default value will be used and you will be notified)")
    
    # Ensure 'model' key exists and has a valid value
    if 'model' not in p.keys():
        p['model'] = 'monomer'
        print_to_log_file(logfilename, "'model' set to default value of 'monomer'")
    elif p['model'] not in ['monomer', 'Frenkel exciton dimer']:
            raise Exception('ERROR: Invalid model requested. (Options: monomer, Frenkel exciton dimer)')
    
    default_params = get_defaults()
    
    # Assign model-specific default parameters
    for key, value in default_params['available models'][p['model']].items():
        if key not in p:
            p[key] = value
            print_to_log_file(logfilename, f"'{key}' set to default value of '{value}'")
    
    # Assign general default parameters
    for key, value in default_params['other'].items():
        if key not in p:
            p[key] = value
            print_to_log_file(logfilename, f"'{key}' set to default value of '{value}'")   
            
    print_to_log_file(logfilename, 'Input successfully read and parameters set.')
 
    # Convert input types
    p = convert_input_types(p)
           
    return p

    # . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . 

def convert_input_types(p):
    """
    Converts input parameter values to their appropriate data types.
    
    Ensures that numeric parameters are properly cast as `float` or `int`, 
    and maintains dictionary structures where necessary.

    Parameters:
        p (dict): Dictionary containing user-defined parameters.

    Returns:
        dict: Updated dictionary with standardized data types.
    """

    int_keys = ['t1_points', 't2_points', 't3_points', 'padding value']
    float_keys = ['t1_step', 't2_step', 't3_step', 'bath_t1t3_disorder', 'bath_t1t3_time', 
                  'bath_t2_disorder', 'bath_t2_time', 'e1', 'e2', 'J', 'temperature']
    
    # Convert simple int and float parameters
    for key in int_keys:
        if key in p:
            p[key] = int(p[key])
    
    for key in float_keys:
        if key in p:
            p[key] = float(p[key])
    
    # Convert dictionary-style parameters (e.g., vibrational parameters)
    vib_param_keys = ['vib_freq', 'displacement', 'vmax']
    
    for key in vib_param_keys:
        if key in p and isinstance(p[key], dict):
            for subkey in p[key]:
                p[key][subkey] = float(p[key][subkey])

    return p

# . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . 

def get_defaults():
    """
    Returns a dictionary containing default system parameters for different models.

    Returns:
        dict: A nested dictionary containing:
            - 'monomer': Default parameters for a monomer system.
            - 'Frenkel exciton dimer': Default parameters for a dimer system.
            - 'other': General default parameters used in all models.
    """
    
    default_dict = {'available models': {}}
       
    # General simulation parameters
    default_dict['other'] = {
        'bath_t1t3_disorder': 1300,  # [cm^-1]
        'bath_t1t3_time': 40,        # [fs]
        'bath_t2_disorder': 125,     # [cm^-1]
        'bath_t2_time': 300,         # [fs]
        't2_step': 5,                # [fs]
        't2_points': 10,
        't1_step': 3,                # [fs]
        't1_points': 62,
        't3_step': 3,                # [fs]
        't3_points': 62,
        'padding value': 256         # Number of points for zero-padding
    }
    
    # Default parameters for a monomer system
    default_dict['available models']['monomer'] = {
        'e1': 14500,                 # Energy of the excited state [cm^-1]
        'vib_freq': {'1': 1300},     # Vibrational frequency [cm^-1]
        'displacement': {'1': 0.3},  # Huang-Rhys factor (dimensionless)
        'vmax': {'1': 5}             # Maximum vibrational quantum number
    }
        
    # Default parameters for a Frenkel exciton dimer
    default_dict['available models']['Frenkel exciton dimer'] = {
        'e1': 14500,                 # Energy of monomer 1 excited state [cm^-1]
        'e2': 14500,                 # Energy of monomer 2 excited state [cm^-1]
        'J': 250,                    # Electronic coupling between monomers [cm^-1]
        'vib_freq': {'1': 1300},     # Vibrational frequency [cm^-1]
        'displacement': {'1': 0.3},  # Huang-Rhys factor (dimensionless)
        'vmax': {'1': 5}             # Maximum vibrational quantum number
    }
    
    return default_dict

# . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . 
