"""
system_model.py

This module defines the system Hamiltonians and basis representations for the 
OREOS optical response simulator. It includes functions for constructing the 
electronic and vibronic Hamiltonians of monomers and Frenkel exciton dimers, 
as well as the creation and annihilation operators needed for quantum mechanical 
calculations. 

Functions:
- get_Hamiltonian: Selects the appropriate Hamiltonian based on input parameters.
- M_Ham: Constructs the vibronic monomer Hamiltonian.
- FE_Ham: Constructs the Frenkel exciton dimer Hamiltonian.
- FEground, FEexcited: Generate basis states for ground and excited manifolds.
- Melecre, FEelecre, Mvibcre, FEvibcre: Define creation operators for electronic 
  and vibrational states.

"""


import numpy as np
import numpy.matlib
from math import pi
from copy import deepcopy
from milk.operational import print_to_log_file

clight = 2.9979*10**(-5) # speed of light[cm/fs]

# . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . 

def get_Hamiltonian(p, logfilename):
    """
    Determines the appropriate Hamiltonian model based on user input.
    
    Parameters:
        p (dict): Dictionary containing system parameters.
        logfilename (str): Path to the log file for error reporting.
    
    Returns:
        dict: Hamiltonian representation of the system.
    """
    
    if p['model'] == 'Frenkel exciton dimer':
        H = FE_dimer_Hamiltonian(p)
        
    elif p['model'] == 'monomer':
        H = monomer_Hamiltonian(p)
    else:
        print_to_log_file(logfilename, 'ERROR: Unsupported molecular model specified')
        raise Exception('ERROR: Unsupported molecular model specified')
        
    return H

# %% Monomer

# . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . 

def monomer_basis(E,nvib,vmax):
    """
    Constructs electronic manifolds for a monomer system.
    
    Parameters:
        E (float): Energy level.
        nvib (int): Number of vibrational modes.
        vmax (dict): Dictionary specifying the maximum vibrational quanta for each mode.
    
    Returns:
        list: Basis states representing the electronic manifolds.
    """
    
    k = 0
    temp = np.zeros([1,1+nvib])
    basis = []
    for i in range(nvib):
        for v in range(int(vmax[str(i+1)])+1):  
            basis.append(temp.copy())
            basis[k][0,0] = E
            basis[k][0,i+1] = v
            k = k+1
            if i>0:
                if v == 0:
                   basis.pop()
                   k = k-1
    return basis

# . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . 

def monomer_creation_electronic_total(basis,states):
    """
    Constructs the electronic creation operator for a monomer system.
    
    Parameters:
        basis (list): Basis states.
        states (int): Number of states in the system.
    
    Returns:
        np.array: Electronic creation operator matrix.
    """
    
    c = np.zeros([states,states],dtype=float)
    for i in range(states):
        for j in range(states):
            bra = basis[i].copy()
            ket = basis[j].copy()
            if (ket[0,0] == 1):
                ket[0,0] = 0 
                if (ket == bra).all():
                    c[i,j] = 1      
    return c

# . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . 

def monomer_creation_electronic(basis,states,num,quanta):
    """
    Constructs the monomer electronic creation operator (vibration-specific).
    
    Parameters:
        basis (list): Basis states.
        states (int): Number of states in the system.
        num (int): Index of the vibrational mode.
        quanta (int): Number of electronic quanta to be created.
    
    Returns:
        np.array: Electronic creation operator matrix.
    """
    
    c = np.zeros([states,states],dtype=float)
    for i in range(states):
        for j in range(states):
            bra = basis[i].copy()
            ket = basis[j].copy()
            if (ket[0,0] == quanta):
                ket[0,0] = 0 
                if np.count_nonzero(ket[0,1:]) == 0:
                    if (ket == bra).all():
                        c[i,j] = 1    
                elif ket[0,num+1] > 0:
                    if (ket == bra).all():
                        c[i,j] = 1   
    return c


# . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . 

def monomer_creation_vibrational_total(basis,states,num):
    """
    Constructs the vibrational creation operator for a monomer system.
    
    Parameters:
        basis (list): Basis states.
        states (int): Number of states in the system.
        num (int): Index of the vibrational mode.
    
    Returns:
        np.array: Vibrational creation operator matrix.
    """
    
    b = np.zeros([states,states],dtype=float)
    for i in range(states):
        for j in range(states):
            bra = basis[i].copy()
            ket = basis[j].copy()
            ket[0,num+1] = ket[0,num+1]-1
            if (ket == bra).all():
               b[i,j] = np.sqrt(ket[0,num+1]+1)    
    return b

# . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . 

def monomer_Hamiltonian(p):
    """
    Constructs the Hamiltonian for a vibronic monomer system.
    
    Parameters:
        p (dict): Dictionary containing system parameters.
    
    Returns:
        dict: Hamiltonian representation of the monomer system with keys:
            - 'H': Full Hamiltonian matrix
            - 'e': Eigenvalues of the Hamiltonian
            - 'v': Eigenvectors of the Hamiltonian
            - 'MU': Transition dipole moment matrix
            - 'Hrf': Rotating frame Hamiltonian
    """
    
    H = {}
    nvib = len(p['vib_freq'])
      
    # Create basis kets (occupancy number basis)
    H['gbasis'] = monomer_basis(0,nvib,p['vmax'])   
    H['ebasis'] = monomer_basis(1,nvib,p['vmax'])
    H['basis'] = H['gbasis'] + H['ebasis']
    H['states'] = len(H['basis'])
    H['fock'] = np.identity(H['states']);
      
    # construct electronic creation operator 
    cm_total = monomer_creation_electronic_total(H['basis'], H['states'])
      
    # construct vibrational creation operator
    bm_total = {}
    cm = {}
    lambda_sqr = {}
    
    for i in range(nvib):
        bm_total[str(i+1)] = monomer_creation_vibrational_total(H['basis'],H['states'],i)
        cm[str(i+1)] = monomer_creation_electronic(H['basis'],H['states'],i,1)
        
        lambda_sqr[str(i+1)] = np.zeros([H['states'],H['states']],dtype=float)
        for j in range(H['states']):
            if np.matmul(cm[str(i+1)].transpose(),cm[str(i+1)])[j,j] == 1:
                lambda_sqr[str(i+1)][j,j] = float(p['displacement'][str(i+1)])**2
  
    # Build Hamiltonian 
    Hetot = float(p['e1'])*np.matmul(cm_total.transpose(), cm_total)    # generate electronic Hamiltonian for m1
    H['H'] = Hetot.copy()
    Hvm = {}
    Hevm = {}
    for i in range(nvib):
        Hvm[str(i+1)] = float(p['vib_freq'][str(i+1)])*np.matmul(bm_total[str(i+1)].transpose(),bm_total[str(i+1)])                                                             
        Hevm[str(i+1)] = float(p['vib_freq'][str(i+1)])*np.matmul(np.matmul(cm[str(i+1)].transpose(),cm[str(i+1)]),float(p['displacement'][str(i+1)])*(bm_total[str(i+1)].transpose()+bm_total[str(i+1)])+lambda_sqr[str(i+1)])   # generate vibronic coupling Hamiltonian for m1
        H['H'] = H['H']+np.sum([Hvm[str(i+1)], Hevm[str(i+1)]], axis=0)
  
    H['e'], H['v'] = np.linalg.eig(H['H'].copy())
    H['MU'] = cm_total.transpose()+cm_total   
    H['Hrf'] = H['H'].copy()
      
    # switch into rotating frame
    for i in range(H['states']):
      if  np.matmul(cm_total.transpose(), cm_total)[i,i] > 0:
          H['Hrf'][i,i] = H['Hrf'][i,i] - float(p['e1'])
          
    # change units      
    H['Hrf'] = clight*H['Hrf']
    H['MU'] = clight*H['MU']
        
    # partition Hamiltonian, Fock, and MU matrices
    n = int(H['states']/2)
    H['Hgg'] = H['Hrf'][0:n,0:n]
    H['Hee'] = H['Hrf'][n:,n:]
    H['fockgg'] = H['fock'][0:n,0:n]
    H['fockee'] = H['fock'][n:,n:]
    H['MUge'] = H['MU'][0:n,n:]
    H['MUeg'] = H['MU'][n:,0:n]
    
    return H

# %% Frenkel exciton dimer

# . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . 

def FE_dimer_ground_basis(nvib,vmax):
    """
    Constructs the electronic ground-state manifold for a Frenkel exciton dimer.
    
    Parameters:
        nvib (int): Number of vibrational modes.
        vmax (dict): Dictionary specifying the maximum vibrational quanta for each mode.
    
    Returns:
        list: Basis states representing the ground-state manifold.
    """
    
    k = 0
    temp = np.zeros([1,2+2*nvib])
    basis = []
    for i in range(nvib):
        for v1 in range(int(vmax[str(i+1)])+1):  
            for v2 in range(int(vmax[str(i+1)])+1-v1):
                basis.append(temp.copy())
                basis[k][0,i+1] = v1
                basis[k][0,i+2+nvib] = v2
                k = k+1
                if i>0:
                    if v1==0 and v2 == 0:
                        basis.pop()
                        k = k-1
    return basis

# . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . 

def FE_dimer_excited_basis(nvib, vmax):
    """
    Constructs the singly excited electronic manifold for a Frenkel exciton dimer.
    
    Parameters:
        nvib (int): Number of vibrational modes.
        vmax (dict): Dictionary specifying the maximum vibrational quanta for each mode.
    
    Returns:
        list: Basis states representing the excited-state manifold.
    """
    
    k = 0
    temp = np.zeros([1,2+2*nvib])
    basis = []
    for i in range(nvib):
        for n in range(2): 
            m = not(n)
            for v1 in range(int(vmax[str(i+1)])+1):  
                for v2 in range(int(vmax[str(i+1)])+1-v1):
                    basis.append(temp.copy())
                    basis[k][0,0] = n
                    basis[k][0,1+nvib] = m
                    basis[k][0,i+1] = v1
                    basis[k][0,i+2+nvib] = v2
                    k = k+1
                    if i>0:
                        if v1==0 and v2 == 0:
                            basis.pop()
                            k = k-1
    return basis

# . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . 

def FE_dimer_electronic_creation_total(basis, states, molnum, nvib, quanta):
    """
    Constructs the electronic creation operator for a Frenkel exciton dimer.
    
    Parameters:
        basis (list): Basis states.
        states (int): Number of states in the system.
        molnum (int): Index of the molecule (1 or 2).
        nvib (int): Number of vibrational modes.
        quanta (int): Number of quanta in the electronic state.
    
    Returns:
        np.array: Electronic creation operator matrix.
    """
    
    c = np.zeros([states,states],dtype=float)
    ind = (nvib+1)*(molnum-1)
    for i in range(states):
        for j in range(states):
            bra = basis[i].copy()
            ket = basis[j].copy()
            if (ket[0,ind] == quanta):
               ket[0,ind] = 0 
               if (ket == bra).all():
                  c[i,j] = 1      
    return c


# . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . 

def FE_dimer_electronic_creation(basis, states, molnum, nvib, num, quanta):
    """
    Constructs the Frenkel exciton dimer electronic creation operator (vibration-specific).
    
    Parameters:
        basis (list): Basis states.
        states (int): Number of states in the system.
        molnum (int): Index of the molecule (1 or 2).
        nvib (int): Number of vibrational modes.
        num (int): Index of the vibrational mode.
        quanta (int): Number of electronic quanta to be created.
    
    Returns:
        np.array: Electronic creation operator matrix.
    """
    
    c = np.zeros([states,states],dtype=float)
    ind = (nvib+1)*(molnum-1)
    
    if molnum == 1:
        conjindex = nvib+1
    elif molnum == 2:
        conjindex = 0
        
    if quanta == 1:
        lownum = 1
        
    for i in range(states):
        for j in range(states):
            bra = basis[i].copy()
            ket = basis[j].copy()
            
            if (ket[0,ind] == quanta):
               ket[0,ind] = ket[0,ind]-lownum
               # ket[0,conjindex] = 0 
               
               vibinfo = ket.copy()
               vibinfo = np.delete(vibinfo,[ind, conjindex])
                   
               if molnum == 1:
                   #first scenario: all vibrations have zero quanta
                   if np.count_nonzero(vibinfo) == 0:
                       if (ket == bra).all():
                           c[i,j] = 1  
                   
                   #second scenario: vibration of interest has non-zero quanta on either monomer
                   elif ket[0,ind+num] > 0 or ket[0,ind+num+nvib+1] > 0:
                       if (ket == bra).all():
                           c[i,j] = 1  
                           
               if molnum == 2:
                   #first scenario: all vibrations have zero quanta
                   if np.count_nonzero(vibinfo) == 0:
                       if (ket == bra).all():
                           c[i,j] = 1  
                   
                   #second scenario: vibration of interest has non-zero quanta on either monomer
                   elif ket[0,ind+num] > 0 or ket[0,num] > 0:
                       if (ket == bra).all():
                           c[i,j] = 1  
    return c


# . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . 

def FE_dimer_vibrational_creation_total(basis,states,molnum,nvib,num):    
    """
    Constructs the vibrational creation operator for a Frenkel exciton dimer.
    
    Parameters:
        basis (list): Basis states.
        states (int): Number of states in the system.
        molnum (int): Index of the molecule (1 or 2).
        nvib (int): Number of vibrational modes.
        num (int): Index of the vibrational mode.
    
    Returns:
        np.array: Vibrational creation operator matrix.
    """
            
    b = np.zeros([states,states],dtype=float)
    ind = (nvib+1)*(molnum-1)+num
    
    for i in range(states):
        for j in range(states):
            bra = basis[i].copy()
            ket = basis[j].copy()

            ket[0,ind] = ket[0,ind]-1 
            if (ket == bra).all():
               b[i,j] = np.sqrt(ket[0,ind]+1)    
    return b

# . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . 

def FE_dimer_Hamiltonian(p):
    """
    Constructs the Hamiltonian for a vibronic dimer (Frenkel exciton model).
    
    Parameters:
        p (dict): Dictionary containing system parameters.
    
    Returns:
        dict: Hamiltonian representation of the dimer system with keys:
            - 'H': Full Hamiltonian matrix
            - 'Hcouplings': Coupling elements in the Hamiltonian
            - 'e': Eigenvalues of the Hamiltonian
            - 'v': Eigenvectors of the Hamiltonian
            - 'MU': Transition dipole moment matrix
            - 'Hrf': Rotating frame Hamiltonian
    """
    
    H = {}
    nvib = len(p['vib_freq'])
    
    # Create basis kets (occupancy number basis)
    H['gbasis'] = FE_dimer_ground_basis(nvib,p['vmax'])   
    H['ebasis'] = FE_dimer_excited_basis(nvib,p['vmax'])
    H['basis'] = H['gbasis']+H['ebasis']
    H['states'] = len(H['basis'])
    H['fock'] = np.identity(H['states']);
    
    # construct electronic operators   
    cm1TOT = FE_dimer_electronic_creation_total(H['basis'],H['states'],1,nvib,1) 
    cm2TOT = FE_dimer_electronic_creation_total(H['basis'],H['states'],2,nvib,1)
    cTOT = np.matmul(cm1TOT.transpose(),cm1TOT)+np.matmul(cm2TOT.transpose(),cm2TOT)  #need np.add?
  
    # construct creation operators for molecules 1 (m1) and 2 (m2) - specific to each vibration
    bm1TOT = {}
    bm2TOT = {}
    cm1 = {}
    cm2 = {}
    
    # displacement matrices for convenience
    lambda_sqr = {}
    for i in range(nvib):
        cm1[str(i+1)] = FE_dimer_electronic_creation(H['basis'],H['states'],1,nvib,i+1,1)
        cm2[str(i+1)] = FE_dimer_electronic_creation(H['basis'],H['states'],2,nvib,i+1,1)
        bm1TOT[str(i+1)] = FE_dimer_vibrational_creation_total(H['basis'],H['states'],1,nvib,i+1)
        bm2TOT[str(i+1)] = FE_dimer_vibrational_creation_total(H['basis'],H['states'],2,nvib,i+1)
        
        lambda_sqr[str(i+1)] = np.zeros([H['states'],H['states']],dtype=float)
        for j in range(H['states']):
          if cTOT[j,j] == 1:
            lambda_sqr[str(i+1)][j,j] = float(p['displacement'][str(i+1)])**2 
  
    # Build Hamiltonian
      
    # purely electronic Hamiltonians (Condon approximation valid)
    Hem1 = float(p['e1'])*np.matmul(cm1TOT.transpose(),cm1TOT)    
    Hem2 = float(p['e2'])*np.matmul(cm2TOT.transpose(),cm2TOT)
    Hetot = np.add(Hem1,Hem2)
    Hem1em2 = float(p['J'])*(np.matmul(cm1TOT.transpose(),cm2TOT)+np.matmul(cm2TOT.transpose(),cm1TOT)) 
  
    # vibronic Hamiltonians
    H['H'] = np.sum([Hetot.copy(),Hem1em2.copy()],axis=0)
    H['Hcouplings'] = Hem1em2.copy()
    Hvm1 = {}
    Hvm2 = {}
    Hevm1 = {}
    Hevm2 = {}
    for i in range(nvib):
        Hvm1[str(i+1)] = float(p['vib_freq'][str(i+1)])*np.matmul(bm1TOT[str(i+1)].transpose(),bm1TOT[str(i+1)])       
        Hvm2[str(i+1)] = float(p['vib_freq'][str(i+1)])*np.matmul(bm2TOT[str(i+1)].transpose(),bm2TOT[str(i+1)])                                                         
        Hevm1[str(i+1)] = float(p['vib_freq'][str(i+1)])*np.matmul(np.matmul(cm1[str(i+1)].transpose(),cm1[str(i+1)]),float(p['displacement'][str(i+1)])*(bm1TOT[str(i+1)].transpose()+bm1TOT[str(i+1)])+lambda_sqr[str(i+1)])   # generate vibronic coupling Hamiltonian for m1
        Hevm2[str(i+1)] = float(p['vib_freq'][str(i+1)])*np.matmul(np.matmul(cm2[str(i+1)].transpose(),cm2[str(i+1)]),float(p['displacement'][str(i+1)])*(bm2TOT[str(i+1)].transpose()+bm2TOT[str(i+1)])+lambda_sqr[str(i+1)])   # generate vibronic coupling Hamiltonian for m2
  
        H['H'] = H['H']+np.sum([Hvm1[str(i+1)], Hvm2[str(i+1)], Hevm1[str(i+1)], Hevm2[str(i+1)]], axis=0)
        H['Hcouplings'] = H['Hcouplings'] + np.sum([Hevm1[str(i+1)], Hevm2[str(i+1)]], axis=0)
  
    H['e'], H['v'] = np.linalg.eig(H['H'].copy())
    H['MU'] = 0.5*(cm1TOT.transpose()+cm1TOT+cm2TOT.transpose()+cm2TOT)    # transition dipole operator
    H['Hrf'] = deepcopy(H['H'])
    H['MU site basis'] = deepcopy(H['MU'])
    
    # switch into rotating frame
    for i in range(H['states']):
      if cTOT[i,i] == 1:
          H['Hrf'][i,i] = H['Hrf'][i,i]-float(p['e1'])
    
    # change units
    c_light = 2.998*10**-5    # [cm/fs]
    H['Hrf'] = (2*pi)*c_light*H['Hrf']
    H['MU'] = (2*pi)*c_light*H['MU']
    
    # partition Hamiltonian, Fock, and MU matrices
    n = int(H['states']/3)
    H['Hgg'] = H['Hrf'][0:n,0:n]
    H['Hee'] = H['Hrf'][n:,n:]
    H['fockgg'] = H['fock'][0:n,0:n]
    H['fockee'] = H['fock'][n:,n:]
    H['MUge'] = H['MU'][0:n,n:]
    H['MUeg'] = H['MU'][n:,0:n]
    
    return H

# . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . 