function output = vibcre_mon_mv_v1(basis,n_states,nvib)

b = zeros(n_states);

for i = 1:n_states
    for j = 1:n_states
        bra = basis(i,:);
        ket = basis(j,:);
        
           ket(nvib+1) = ket(nvib+1)-1; 

           tf = isequal(ket,bra);
           if tf>0
                b(i,j) = sqrt(ket(nvib+1)+1); 
           end
    end
end
     

output = b;
end