function [p] = modelparameters(model)

p.c = 2.9979 * 10^-5;      % speed of light [cm/fs]

switch model
    
    case 'monomer'
        p.model = model;
        p.nelec = 2;
        p.e1 = 17550;            % center energy of local excited electronic state [cm-1] 
        p.wfreq = [1350];   % vibrations coupled to electronic state [cm-1]
        p.lambda = [0.5];    % square root of Huang-Rhys factor for each vibration
        p.vmax = [2];         % maximum vibrational quanta for each vibration;
        
    case 'FE dimer'
        p.nelec = 3;
        p.e1 = 10000;            % center energy of molecule 1 local excited electronic state [cm-1] 
        p.e2 = 9000;             % center energy of molecule 2 local excited electronic state [cm-1] 
        p.Jcoul = 100;           % Coulombic dipole-dipole coupling
        p.wfreq = 1000;          % vibrations coupled to electronic state [cm-1]  
        p.lambda = 0.5;          % square root of Huang-Rhys factor for each vibration
        p.vmax = 3;              % maximum vibrational quanta for each vibration;  
        
    otherwise
        error('Error: invalid model specification')
end

end