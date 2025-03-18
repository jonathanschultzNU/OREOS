# OREOS (Optical REspOnse Simulator)

OREOS is a Python-based simulation package designed for computing third-order nonlinear optical response functions. It enables the study of ultrafast spectroscopy by generating response functions in the time domain, transforming them into the frequency domain, and applying post-processing techniques such as lineshape modeling and windowing.

## Key Features

- **System Hamiltonian Construction**: Supports monomer and Frenkel exciton dimer models with electronic and vibronic coupling.
- **Third-Order Nonlinear Response Functions**: Simulates four key response functions relevant to 2D spectroscopy.
- **Fourier Transform Processing**: Converts time-domain response functions into frequency-domain spectra.
- **Lineshape and Windowing Corrections**: Implements dephasing and windowing techniques to improve spectral features.
- **Logging and Input Validation**: Ensures consistency in input parameters and provides detailed logs.
- **Modular Design**: Organized into distinct modules for Hamiltonian generation, simulations, signal processing, and file handling.

## Modules Overview

- `main.py` - Main execution script that initializes directories, processes input files, and runs simulations.
- `system_model.py` - Constructs the system Hamiltonian for monomers and dimers.
- `simulations.py` - Computes response functions and processes time-domain signals.
- `signal_processing.py` - Applies Fourier transforms, normalization, and other signal analysis techniques.
- `defaults.py` - Provides default parameters and input validation.
- `operational.py` - Handles logging, input file parsing, and saving results.

---

## Setup and Installation

To use the OREOS Optical Response Simulator, follow these steps to set up the environment.

### 1. Clone or Download the Repository
If you haven't already, either download the repository or clone it to your local machine:
```bash
git clone https://github.com/your-repo/oreos-simulator.git
cd oreos-simulator
```
### 2. Create and Activate the Conda Environment
Use the provided environment.yml file to set up a Conda environment with all required dependencies:
```bash
conda env create -f environment.yml
conda activate oreos-simulator
```

### 3. Preparing for Simulation
Generate a .txt input file according to the input parameters described below. Place the input file in a directory called 'inputs'

### 4. Running the Simulation
To run the simulation, execute:

```bash
python main.py <runname>
```
where <runname> is the identifier for the input and output files. To run the included sample_input job, execute:

```bash
python main.py sample_input
```

---

## Input Parameters

The simulation requires a set of parameters provided in a dictionary format. Below is a summary of key parameters:

| Parameter | Type | Description |
|-----------|------|-------------|
| `model` | `str` | Specifies the system model (`monomer` or `Frenkel exciton dimer`). |
| `e1` | `float` | Excitation energy of monomer 1 (cm⁻¹). |
| `e2` | `float` | Excitation energy of monomer 2 (cm⁻¹) (dimer model only). |
| `J` | `float` | Electronic coupling between monomers (cm⁻¹) (dimer model only). |
| `vib_freq` | `dict` | Vibrational mode frequencies (cm⁻¹). |
| `displacement` | `dict` | Huang-Rhys factors defining vibrational-electronic coupling. |
| `vmax` | `dict` | Maximum vibrational quantum number for each mode. |
| `t1_step`, `t2_step`, `t3_step` | `float` | Time steps for each delay time (fs). |
| `t1_points`, `t2_points`, `t3_points` | `int` | Number of time points for each delay. |
| `bath_t1t3_disorder`, `bath_t2_disorder` | `float` | Spectral diffusion parameters (cm⁻¹). |
| `bath_t1t3_time`, `bath_t2_time` | `float` | Bath correlation times (fs). |
| `padding value` | `int` | Number of points for zero-padding in Fourier transforms. |


## Sample Input File:

model = Frenkel exciton dimer
e1 = 14500
vib_freq = 1300.0
vmax = 3
e2 = 14500
J = 250
displacement = 0.6
t2_step = 5
t2_points = 10
t1_step = 3
t1_points = 62
t3_step = 3
t3_points = 62
padding value = 256
bath_t1t3_disorder = 1300
bath_t1t3_time = 40
bath_t2_disorder = 125
bath_t2_time = 300

---

## Output
Results are saved in the outputs/ directory and include:

A log file (<runname>.log) with execution details.
A results file (<runname>.pkl) containing simulation data.

---

## **Citing This Work**
If you use this code in your research, please cite:
[manuscript citation info pending]

### **License**
This project is licensed under the MIT License. See LICENSE for details.

### **Contact**
For questions or contributions, please reach out to the developers.