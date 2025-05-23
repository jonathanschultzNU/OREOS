"""
signal_processing.py

This module provides functions for processing and analyzing 2D optical spectra 
in the OREOS optical response simulator. It includes normalization utilities and  
Fourier transformations.

Functions:
- multcompreal: Multiplies a complex-valued array with a real-valued array.
- normalize2d: Normalizes a 2D array to its maximum absolute value.
- normalize_array_0_to_1: Scales a 1D array to the range [0,1].
- response_FFT: Computes the 2D Fourier transform of response functions.
"""

import numpy as np

clight = 2.9979*10**(-5) # speed of light[cm/fs]

# . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . 

def multcompreal(mat1,mat2):
    """
    Multiplies a complex-valued array with a real-valued array.

    Parameters:
        mat1 (np.array): Complex-valued numpy array (dtype = complex128).
        mat2 (np.array): Real-valued numpy array (dtype = float64).

    Returns:
        np.array: Complex-valued array resulting from element-wise multiplication.
    """
    realpart = np.multiply(np.real(mat1),mat2)
    imagpart = np.multiply(np.imag(mat1),mat2)
    outmat = realpart+1j*imagpart
    return outmat

# . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . 

def normalize2d(matrix):
    """
    Normalizes a 2D matrix by its maximum absolute value.
    Parameters:   matrix (np.array): Input 2D numpy array.
    Returns:      np.array: Normalized 2D array.
    """
    absmat = np.abs(matrix)
    outmat = matrix/absmat.max()
    return outmat

# . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . 

def normalize_array_0_to_1(arr):
    """
    Scales a 1D array to the range [0,1].
    Parameters:   arr (np.array): Input 1D numpy array.
    Returns:      np.array: Normalized array in the range [0,1].
    """
    
    normalized_arr = (arr - arr.min()) / (arr.max() - arr.min())
    return normalized_arr

# . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . 

def response_FFT(p, R, Rt1t3d):
    """
    Computes the 2D Fourier transform of the response functions.

    Applies a trapezoidal rule correction (Hamm & Zanni, Ch. 9) before computing 
    the Fourier transform and organizing the response functions into different 
    contributions.

    Parameters:
        p (dict): Simulation parameters.
        R (dict): Dictionary containing response function time grids.
        Rt1t3d (dict): Dictionary containing pre-processed response functions.

    Returns:
        tuple:
            - R (dict): Updated dictionary with frequency-domain response functions.
            - Rt1t3c (dict): Corrected time-domain response functions.
            - Rw1w3 (dict): Frequency-domain response functions.
    """
    
    # Apply the trapezoidal rule correction (see Hamm and Zanni ch. 9)
    Rt1t3c = {}
    for j in range(R['rnum']):
        Rt1t3c[str(j+1)] = Rt1t3d[str(j+1)].copy()
        Rt1t3c[str(j+1)][:,0,:] = Rt1t3c[str(j+1)][:,0,:]/2 
        Rt1t3c[str(j+1)][0,1:,:] = Rt1t3c[str(j+1)][0,1:,:]/2
    
    Rw1w3 = {}
    for j in range(R['rnum']):
        Rw1w3[str(j+1)] = np.zeros([R['nw1'],R['nw3'],R['nt2']],dtype=complex)

    for k in range(R['nt2']):
        for j in range(R['rnum']):
            temp = np.fft.fft2(np.pad(Rt1t3c[str(j+1)][:,:,k].copy(), (0, R['nw3']-R['nt3']), 'constant'))
            if j+1 == 1 or j+1 == 4:
                Rw1w3[str(j+1)][:,:,k] = np.fft.fftshift(np.rot90(temp,2))
                
            elif j+1 == 2 or j+1 == 3:
                temp2 = np.fliplr(np.roll(temp.copy(),-1,axis=1)) 
                Rw1w3[str(j+1)][:,:,k] = np.fft.fftshift(np.rot90(temp2,2))
    
    R.update({'Rw1w3Rephasing': np.sum([Rw1w3['2'], 
                                        Rw1w3['3']], axis=0),
              'Rw1w3NonRephasing': np.sum([Rw1w3['1'], 
                                           Rw1w3['4']], axis=0),
              'Rw1w3Absorptive': np.real(np.sum([Rw1w3['1'], 
                                                 Rw1w3['2'], 
                                                 Rw1w3['3'], 
                                                 Rw1w3['4']], axis=0))
              })
    
    maxval = np.max(np.abs(R['Rw1w3Absorptive']))
    
    R['Rw1w3 normalized'] = {'Rw1w3Absorptive': R['Rw1w3Absorptive']/maxval,
                             'Rw1w3Rephasing': multcompreal(R['Rw1w3Rephasing'], 1/maxval),
                             'Rw1w3NonRephasing': multcompreal(R['Rw1w3NonRephasing'], 1/maxval)
                             }
    
    return R, Rt1t3c, Rw1w3

# . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . 