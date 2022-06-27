function output = vibcre_dim_2v(basis,n_states,molnum,nvib,numvib)

b = zeros(n_states);

for i = 1:n_states
    for j = 1:n_states
        bra = basis(i,:);
        ket = basis(j,:);
       
        ket(molnum+nvib+numvib*(molnum-1)) = ket(molnum+nvib+numvib*(molnum-1))-1; 

           tf = isequal(ket,bra);
           if tf>0
                b(i,j) = sqrt(ket(molnum+nvib+numvib*(molnum-1))+1); 
           end
    end
end
     

output = b;