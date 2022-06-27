function Hmat = RSMonomer(p)

%   Description: Formulation of Frenkel-Exciton Hamiltonian for a dimer with multiple vibrations

c = 2.9979 * 10^-5;      % speed of light [cm/fs]
p.numvib = length(p.wfreq);

%Create site basis set
Hmat.basis_g = monomer_basis(0, p.numvib,p.vmax);
Hmat.basis_s = monomer_basis(1, p.numvib,p.vmax);
Hmat.basis = [Hmat.basis_g; Hmat.basis_s];

Hmat.n_states = length(Hmat.basis(:,1));
Hmat.fock = eye(Hmat.n_states);

% construct electronic operators   
cs1 = elecreTot_mon(Hmat.basis,Hmat.n_states,1);  %both vibrations share a common ground excited state - don't want to double the diagonal energy

ctots1 = (cs1')*cs1;
b = cell(1,p.numvib);
s1 = cell(1,p.numvib);
lambda_sqr = cell(1,p.numvib);

for i=1:p.numvib
    
    b{i} = vibcre_mon(Hmat.basis,Hmat.n_states,i); %construct vibrational annhiliation operators
    s1{i} = elecre_mon(Hmat.basis,Hmat.n_states,p.numvib,i,1);
    
    lambda_sqr{i} = zeros(Hmat.n_states,Hmat.n_states);
    for k = 1:Hmat.n_states
        
        if ctots1(k,k) == 1
            lambda_sqr{i}(k,k) = (p.lambda(i)^2);             %Matrix of Huang-Rhys factor squared for Hamiltonian term
        end
    end
    
end


% Create Hamiltonian piecewise

Hs1tot = p.e1*(cs1')*cs1;     %Hamiltonian term 1: electronic states

Hmat.H = zeros(Hmat.n_states);

Hvib = cell(1,p.numvib);
Hs1vib = cell(1,p.numvib);

for i=1:p.numvib

    Hvib{i} = p.wfreq(i)*(b{i}')*b{i};
    Hs1vib{i} = p.wfreq(i)*(s1{i}')*s1{i}*(p.lambda(i).*(b{i}'+b{i})+lambda_sqr{i});
    Hmat.H = Hmat.H+Hvib{i}+Hs1vib{i};

end

Hmat.H = Hmat.H+Hs1tot;

%construct transition dipole operators
mu1E = 1;   
Hmat.mus1 = 0.5*mu1E*(cs1'+cs1);
Hmat.mu = Hmat.mus1;

%Rotate all matrices and vectors to eigenbasis of Hamiltonian
[Hmat.V,Hmat.D] = eig(Hmat.H);    %V being the rotation matrix, D being the diagonalized Hamiltonian
% [Hmat.V,Hmat.D] = eigenshuffle(Hmat.H);
% Removal of fast-oscillating frequencies from diagonal elements of Hamiltonian; 
% added back to Fourier Transform axis before plotting

for i=1:Hmat.n_states
    if ctots1(i,i) > 0
        Hmat.H(i,i) = Hmat.H(i,i)-p.e1;
    end
end

Hmat.mu = Hmat.mu*c;
Hmat.Hwvn = Hmat.H;
Hmat.H = Hmat.H*c;
BC = @(A) Hmat.V'*A*Hmat.V;           %Define change of basis function
Hmat.mueig = BC(Hmat.mu);   %Rotate transition dipole matrix   

n = Hmat.n_states/(p.nelec);


Hmat.Hgg = Hmat.H(1:n,1:n);
Hmat.siteGG = diag(Hmat.Hgg);

Hmat.Hee = Hmat.H(n+1:2*n,n+1:2*n);
Hmat.siteSE = diag(Hmat.Hwvn(n+1:2*n,n+1:2*n));

Hmat.eigvals = diag(Hmat.H);

Hmat.fockgg = Hmat.fock(1:n,1:n);
Hmat.fockee = Hmat.fock(n+1:2*n,n+1:2*n);
Hmat.MUge = Hmat.mu(1:n,n+1:2*n);
Hmat.MUeg = Hmat.mu(n+1:2*n,1:n);
end



