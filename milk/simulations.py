"""
simulations.py

This module contains functions for simulating time-domain optical response functions,
processing signal data, and calculating linear response spectra for the OREOS optical 
response simulator.

Functions:
- response_simulation: Computes third-order nonlinear response functions.
- process_time_domain_response: Applies lineshapes, normalization, and windowing to response functions.
- linear_response: Computes the linear absorption spectrum in the time domain.
- direct_linear_response: Computes the absorption spectrum using oscillator strengths and Gaussian broadening.

"""

import numpy as np
import numpy.matlib
import scipy.optimize
from milk.signal_processing import multcompreal, response_FFT
from milk.operational import print_to_log_file

clight = 2.9979*10**(-5) # speed of light[cm/fs]

# . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . 

def response_simulation(p, H, logfilename):
    
    """
    Simulates the third-order nonlinear optical response functions in the time domain.

    Response functions:
    - R1: Non-rephasing ground state bleach
    - R2: Rephasing stimulated emission
    - R3: Rephasing ground state bleach
    - R4: Non-rephasing stimulated emission

    Parameters:
        p (dict): Dictionary containing simulation parameters.
        H (dict): Hamiltonian representation of the system.
        logfilename (str): Path to the log file for tracking simulation progress.

    Returns:
        tuple:
            - R (dict): Dictionary of time-domain response function parameters and grids.
            - Rt1t3 (dict): Time-domain third-order response functions indexed as '1' to '4'.
    """
    
    # Initialize response function dictionary
    R = {
        'rnum': 4,  # Number of response functions
        'dt1': p['t1_step'],
        'dt2': p['t2_step'],
        'dt3': p['t3_step'],
        't1_max': p['t1_step'] * p['t1_points'],
        't2_max': p['t2_step'] * p['t2_points'],
        't3_max': p['t3_step'] * p['t3_points'],
        }
    
    R['t1'] = np.arange(0, R['t1_max'], R['dt1'])
    R['t2'] = np.arange(0, R['t2_max'], R['dt2'])
    R['t3'] = np.arange(0, R['t3_max'], R['dt3'])
    
    R['nt1'], R['nt2'], R['nt3'] = len(R['t1']), len(R['t2']), len(R['t3'])

    # Define frequency w3 axis for Fourier Transform
    padnum = p['padding value']
    R['nw1'] = R['nw2'] = R['nw3'] = padnum
    R['dw1'] = 1/(clight*R['dt1']*R['nw1'])
    R['dw2'] = 1/(clight*R['dt2']*R['nw2'])
    R['dw3'] = 1/(clight*R['dt3']*R['nw3'])

    # frequency axes [cm-1] (corrected for rotating frame)
    R['w1'] = np.arange((-R['nw1']*R['dw1'])/2, (R['nw1']*R['dw1'])/2, R['dw1']) + p['e1']
    R['w2'] = np.arange((-R['nw2']*R['dw2'])/2, (R['nw2']*R['dw2'])/2, R['dw2'])
    R['w3'] = np.arange((-R['nw3']*R['dw3'])/2, (R['nw3']*R['dw3'])/2, R['dw3']) + p['e1']

        
    # define initial state
    ket = H['fockgg'][:,0].copy()  # assuming no thermal distribution
    bra = ket.transpose()

    Rt1t3 = {}
    for j in range(R['rnum']):
        Rt1t3[str(j+1)] = np.zeros([R['nt1'],R['nt3'],R['nt2']], dtype=np.complex128)

    # add dimension to Hgg, Hee to make Hggni* and Heeni*
    Hggni3 = np.repeat(H['Hgg'][np.newaxis,...], R['nt3'], axis=0) 
    Heeni3 = np.repeat(H['Hee'][np.newaxis,...], R['nt3'], axis=0)

    Hggni2 = np.repeat(H['Hgg'][np.newaxis,...], R['nt2'], axis=0)
    Heeni2 = np.repeat(H['Hee'][np.newaxis,...], R['nt2'], axis=0)

    Hggni1 = np.repeat(H['Hgg'][np.newaxis,...], R['nt1'], axis=0)
    Heeni1 = np.repeat(H['Hee'][np.newaxis,...], R['nt1'], axis=0)

    # prepare Hggnt* and Heent* to be same size as Hggni* and Heeni* and complex
    lg1 = H['Hgg'].shape[0]
    le1 = H['Hee'].shape[0]

    Hggnt3 = np.empty((R['nt3'],lg1,lg1)).astype(complex)
    Heent3 = np.empty((R['nt3'],le1,le1)).astype(complex)

    Hggnt2 = np.empty((R['nt2'],lg1,lg1)).astype(complex)
    Heent2 = np.empty((R['nt2'],le1,le1)).astype(complex)

    Hggnt1 = np.empty((R['nt1'],lg1,lg1)).astype(complex)
    Heent1 = np.empty((R['nt1'],le1,le1)).astype(complex)

    # populate Hggnt* and Heent* with Unitary eqn.
    for k in range(R['nt3']):
        Hggnt3[k,:,:] = scipy.linalg.expm(-1j*Hggni3[k,:,:]*R['t3'][k])
        Heent3[k,:,:] = scipy.linalg.expm(-1j*Heeni3[k,:,:]*R['t3'][k])

    for k in range(R['nt2']):
        Hggnt2[k,:,:] = scipy.linalg.expm(-1j*Hggni2[k,:,:]*R['t2'][k])
        Heent2[k,:,:] = scipy.linalg.expm(-1j*Heeni2[k,:,:]*R['t2'][k])

    for k in range(R['nt1']):
        Hggnt1[k,:,:] = scipy.linalg.expm(-1j*Hggni1[k,:,:]*R['t1'][k])
        Heent1[k,:,:] = scipy.linalg.expm(-1j*Heeni1[k,:,:]*R['t1'][k])

    # add two dimensions to Hggnt* and Heent* so have shape: [nt3, nt1, nt2, :, :]
    Hggnt3 = np.repeat(Hggnt3[:,np.newaxis,...], R['nt1'], axis=1)
    Heent3 = np.repeat(Heent3[:,np.newaxis,...], R['nt1'], axis=1)

    Hggnt1 = np.repeat(Hggnt1[np.newaxis,...], R['nt3'], axis=0)
    Heent1 = np.repeat(Heent1[np.newaxis,...], R['nt3'], axis=0)

    # generate beginning (lhs) and end (rhs) parts of several response functions -- the dimensions need to work with the Hggnt* and Heent* stacks of matrices

    lhs1 = bra
    lhs1 = np.repeat(lhs1[np.newaxis,...], 1, axis=0)
    lhs1 = np.repeat(lhs1[np.newaxis,...], R['nt1'], axis=0)
    lhs1 = np.repeat(lhs1[np.newaxis,...], R['nt3'], axis=0)

    rhs1 = ket
    rhs1 = np.repeat(rhs1[:,np.newaxis], 1, axis=1)
    rhs1 = np.repeat(rhs1[np.newaxis,...], R['nt1'], axis=0)
    rhs1 = np.repeat(rhs1[np.newaxis,...], R['nt3'], axis=0)
    
    print_to_log_file(logfilename,'Starting third-order response simulations...')

    for k in range(R['nt2']):
        '''
        Steps to generate response functions R#:
            -calculate 4D tensor iR#t1t3 (dim [nt3, nt1, 1, 1])
            -flatten to R#t1t3 (dim [nt3, nt1, nt2])
            -delete 4D tensor to save memory
            
        Equation forms:
            R1(t1,t2,t3) = <G| UGG†(t1) UGG†(t2) UGG†(t3) μEG† UEE(t3) μGE UGG(t2) μEG† UEE(t1) μGE |G>
            R2(t1,t2,t3) = <G| μGE† UEE†(t1) UEE†(t2) μEG UGG†(t3) μEG† UEE(t3) UEE(t2) μGE UGG(t1) |G>
            R3(t1,t2,t3) = <G| μGE† UEE†(t1) μEG UGG†(t2) UGG†(t3) μEG† UEE(t3) μGE UGG(t2) UGG(t1) |G>
            R4(t1,t2,t3) = <G| UGG†(t1) μGE† UEE†(t2) μEG UGG†(t3) μEG† UEE(t3) UEE(t2) UEE(t1) μGE |G>

        '''
        
        iR1t1t3 = np.matmul(lhs1, np.matmul(np.conj(Hggnt1), np.matmul(np.conj(Hggnt2[k]), np.matmul(np.conj(Hggnt3), np.matmul(H['MUeg'].transpose(), np.matmul(Heent3, np.matmul(H['MUge'].transpose(), np.matmul(Hggnt2[k], np.matmul(H['MUeg'].transpose(), np.matmul(Heent1, np.matmul(H['MUge'].transpose(), rhs1)))))))))))
        Rt1t3['1'][:,:,k] = iR1t1t3[:,:,0,0]
        del iR1t1t3
        
        iR2t1t3 = np.matmul(lhs1, np.matmul(H['MUge'], np.matmul(np.conj(Heent1), np.matmul(np.conj(Heent2[k]), np.matmul(H['MUeg'], np.matmul(np.conj(Hggnt3), np.matmul(H['MUeg'].transpose(), np.matmul(Heent3, np.matmul(Heent2[k], np.matmul(H['MUge'].transpose(), np.matmul(Hggnt1, rhs1)))))))))))
        Rt1t3['2'][:,:,k] = iR2t1t3[:,:,0,0]
        del iR2t1t3
        
        iR3t1t3 = np.matmul(lhs1, np.matmul(H['MUge'], np.matmul(np.conj(Heent1), np.matmul(H['MUeg'], np.matmul(np.conj(Hggnt2[k]), np.matmul(np.conj(Hggnt3), np.matmul(H['MUeg'].transpose(), np.matmul(Heent3, np.matmul(H['MUge'].transpose(), np.matmul(Hggnt2[k], np.matmul(Hggnt1, rhs1)))))))))))
        Rt1t3['3'][:,:,k] = iR3t1t3[:,:,0,0]
        del iR3t1t3
        
        iR4t1t3 = np.matmul(lhs1, np.matmul(np.conj(Hggnt1),np.matmul(H['MUge'], np.matmul(np.conj(Heent2[k]), np.matmul(H['MUeg'], np.matmul(np.conj(Hggnt3), np.matmul(H['MUeg'].transpose(), np.matmul(Heent3, np.matmul(Heent2[k], np.matmul(Heent1, np.matmul(H['MUge'].transpose(), rhs1)))))))))))
        Rt1t3['4'][:,:,k] = iR4t1t3[:,:,0,0]
        del iR4t1t3
    
    print_to_log_file(logfilename, "Third-order response simulations complete!")
    
    return R, Rt1t3

# . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . 

def process_time_domain_response(p, R, Rt1t3, H, logfilename):
    """
    Applies post-processing steps to time-domain response functions:
    - Applies lineshape function (dephasing effects)
    - Normalizes response functions
    - Applies a window function

    Parameters:
        p (dict): Dictionary containing simulation parameters.
        R (dict): Time-domain response function parameters.
        Rt1t3 (dict): Raw response functions before post-processing.
        H (dict): Hamiltonian representation of the system.
        logfilename (str): Path to the log file for tracking progress.

    Returns:
        dict: Processed results containing response functions and system parameters.
    """
    
    print_to_log_file(logfilename,'Beginning signal processing...')
    
    b = {'bfluct': p['bath_t1t3_disorder'] * clight,     # [1/fs]
         'btime': p['bath_t1t3_time'],                   # [fs]
         'bt2fluct': p['bath_t2_disorder'] * clight,     # [1/fs]
         'bt2time': p['bath_t2_time'],                   # [fs]
         'T1mesh': np.meshgrid(R['t1'],R['t1'], indexing='ij'), # Generate arrays for lineshape function (w1 and w3 domains)
         'T3mesh': np.meshgrid(R['t3'],R['t3'], indexing='ij')
         }                 

    b['t1t3ls'] = np.exp(-(b['bfluct'] ** 2) * (b['btime'] ** 2) * (
                  np.exp(-((b['T1mesh'][1] + b['T3mesh'][0]) / b['btime'])) + (
                  (b['T1mesh'][1] + b['T3mesh'][0]) / b['btime'] - 1)))
    
    b['t1t3ls3D'] = np.dstack([b['t1t3ls']] * R['nt2'])

           
    b['t2ls'] = np.exp(-(b['bt2fluct']**2)*(b['bt2time']**2)*(np.exp(-(R['t2'])/b['bt2time'])+((R['t2'])/b['bt2time'])-1))
    b['t2ls2'] = np.ones([R['nt1'], R['nt3'], R['nt2']], dtype=float)
    b['t2ls3D'] = np.multiply(b['t2ls2'], b['t2ls'])
    b['relaxmat'] = np.multiply(b['t1t3ls3D'].copy(), b['t2ls3D'].copy())

    # Normalizing to absolute maximum among all response functions (assuming this occurs at t2 = 0)
    stacked_arrays = np.dstack((np.abs(Rt1t3['1'][:,:,0]), 
                                np.abs(Rt1t3['2'][:,:,0]), 
                                np.abs(Rt1t3['3'][:,:,0]), 
                                np.abs(Rt1t3['4'][:,:,0])))
    max_of_stack = np.max(stacked_arrays.max(2))

    # Generate window function
    window_sym = np.hanning(2*R['nt1'])**0.5                        # create the symmetric Hann window
    window_half = np.flipud(window_sym[:len(window_sym)//2])        # cut the Hann Window in half to match FID form
    window_matrix = np.outer(window_half,window_half)
    b['window3D'] = np.dstack([window_matrix] * R['nt2'])

    # Normalization (to absolute max across all four response functions)
    Rt1t3n = {}                 
    for j in range(R['rnum']):
        Rt1t3n[str(j+1)] = multcompreal(Rt1t3[str(j+1)].copy(),1/max_of_stack)

    # Lineshapes
    Rt1t3d_temp = {}            
    for j in range(R['rnum']):
        Rt1t3d_temp[str(j+1)] = multcompreal(-Rt1t3n[str(j+1)].copy(),b['relaxmat'].copy())

    # Windowing
    Rt1t3d = {}                 
    for j in range(R['rnum']):
        Rt1t3d[str(j+1)] = multcompreal(Rt1t3d_temp[str(j+1)].copy(),b['window3D'].copy())    # Apply window function along t1 and t3
        
    # Warn if windowed data are significantly different than un-windowed data
    win_resid = np.sum(np.abs(np.real(Rt1t3d['1'][:,:,0])) - np.abs(np.real(Rt1t3d_temp['1'][:,:,0])))
    win_resid_norm = np.abs(win_resid)/np.sum(np.abs(np.real(Rt1t3d_temp['1'][:,:,0])))
    if win_resid_norm > 0.14:
        print_to_log_file(logfilename, '''WARNING:
    Bath parameters resulted in weak dephasing along t1 and t3. The lineshape may therefore be  
    dominated by the window function. Consider increasing the strength of the lineshape function.''')  
        
    # Convert time-domain signals into frequency domain
    R, Rt1t3c, Rw1w3 = response_FFT(p, R, Rt1t3d)
    
    R['Rt1t3'] = Rt1t3
    R['Rw1w3'] = Rw1w3
    
    keys_to_remove = ['t1_max', 't2_max', 't3_max', 'w2']
    for key in keys_to_remove:
        del R[key]
    
    results = {'p': p,
               'H': H,
               'R': R}
    
    print_to_log_file(logfilename,'Signal processing complete!')
    
    return results

# . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . 