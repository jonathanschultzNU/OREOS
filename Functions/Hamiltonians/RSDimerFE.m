function Hmat = RSDimerFE(p)

%   Description: Formulation of Frenkel-Exciton Hamiltonian for a dimer with multiple vibrations

c = 2.9979 * 10^-5;      % speed of light [cm/fs]
p.numvib = length(p.wfreq);

%Create site basis set
Hmat.basis_g = ground(p.numvib,p.vmax);
Hmat.basis_s = singlyexcited(p.numvib,p.vmax);
Hmat.basis = [Hmat.basis_g; Hmat.basis_s];

Hmat.n_states = length(Hmat.basis(:,1));
Hmat.fock = eye(Hmat.n_states);

% construct electronic operators   
ctotm1 = elecreTot_dim(Hmat.basis,Hmat.n_states,1,p.numvib,1);  %both vibrations share a common ground excited state - don't want to double the diagonal energy
ctotm2 = elecreTot_dim(Hmat.basis,Hmat.n_states,2,p.numvib,1);

ctot = (ctotm1')*ctotm1+(ctotm2')*ctotm2;
b1 = cell(1,p.numvib);
b2 = cell(1,p.numvib);
s1 = cell(1,p.numvib);
s2 = cell(1,p.numvib);
lambda_sqr = cell(1,p.numvib);

for i=1:p.numvib
    
    b1{i} = vibcre_dim(Hmat.basis,Hmat.n_states,1,i,p.numvib); %construct vibrational annhiliation operators
    b2{i} = vibcre_dim(Hmat.basis,Hmat.n_states,2,i,p.numvib);
    s1{i} = elcre_dim(Hmat.basis,Hmat.n_states,1,i,p.numvib,1);
    s2{i} = elcre_dim(Hmat.basis,Hmat.n_states,2,i,p.numvib,1);
    
    lambda_sqr{i} = zeros(Hmat.n_states,Hmat.n_states);
    for k = 1:Hmat.n_states
        
        if ctot(k,k) == 1
            lambda_sqr{i}(k,k) = (p.lambda(i)^2);             %Matrix of Huang-Rhys factor squared for Hamiltonian term
        end
    end
    
end


% Create Hamiltonian piecewise

He1tot = p.e1*(ctotm1')*ctotm1;     %Hamiltonian term 1: electronic states
He2tot = p.e2*(ctotm2')*ctotm2;     %Hamiltonian term 1: electronic states
He1e2 = p.Jcoul*((ctotm1')*ctotm2+(ctotm2')*ctotm1);     %Hamiltonian term 1: electronic states

Hmat.H = zeros(Hmat.n_states);

Hvib1 = cell(1,p.numvib);
Hvib2 = cell(1,p.numvib);
Helvib1 = cell(1,p.numvib);
Helvib2 = cell(1,p.numvib);

for i=1:p.numvib

    Hvib1{i} = p.wfreq(i)*(b1{i}')*b1{i};
    Hvib2{i} = p.wfreq(i)*(b2{i}')*b2{i};
    Helvib1{i} = p.wfreq(i)*(s1{i}')*s1{i}*(p.lambda(i).*(b1{i}'+b1{i})+lambda_sqr{i});
    Helvib2{i} = p.wfreq(i)*(s2{i}')*s2{i}*(p.lambda(i).*(b2{i}'+b2{i})+lambda_sqr{i});
    Hmat.H = Hmat.H+Hvib1{i}+Hvib2{i}+Helvib1{i}+Helvib2{i};

end

Hmat.H = Hmat.H+He1tot+He2tot+He1e2;

%construct transition dipole operators
mu1E = 1;   
Hmat.mu = 0.5*mu1E*(ctotm1'+ctotm1+ctotm2'+ctotm2);

%Rotate all matrices and vectors to eigenbasis of Hamiltonian
[Hmat.V,Hmat.D] = eig(Hmat.H);    %V being the rotation matrix, D being the diagonalized Hamiltonian

% Removal of fast-oscillating frequencies from diagonal elements of Hamiltonian; 
% added back to Fourier Transform axis before plotting

for i=1:Hmat.n_states
    if ctot(i,i) > 0
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
m = length(Hmat.siteGG);

Hmat.Hee = Hmat.H(n+1:end,n+1:end);
Hmat.siteSE = diag(Hmat.Hwvn(n+1:end,n+1:end));
Hmat.siteSE1 = Hmat.siteSE(1:m);
Hmat.siteSE2 = Hmat.siteSE(m+1:end);

Hmat.eigvals = diag(Hmat.H);

Hmat.fockgg = Hmat.fock(1:n,1:n);
Hmat.fockee = Hmat.fock(n+1:end,n+1:end);
Hmat.MUge = Hmat.mu(1:n,n+1:end);
Hmat.MUeg = Hmat.mu(n+1:end,1:n);

end


