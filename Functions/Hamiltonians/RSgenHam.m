function Hmat = RSgenHam(model,p)

switch model
    case 'monomer'
        Hmat = RSMonomer(p);

    case 'FE dimer'
        Hmat = RSDimerFE(p);  
        
    otherwise
        error('Error: invalid model specification')
end

end