function output = elcregenMON_v1(basis,n_states,numvib,nvib,quanta)

c = zeros(n_states);

for i = 1:n_states
    for j = 1:n_states
        bra = basis(i,:);
        ket = basis(j,:);
        
        if ket(1) == quanta
           ket(1) = ket(1)-1; 
           
           if ket(2:numvib+1) == 0 
                tf = isequal(ket,bra);
                if tf>0   
                    c(i,j) = 1;       
                end
                
           elseif ket(nvib+1) > 0 
                tf = isequal(ket,bra);
                if tf>0   
                    c(i,j) = 1;       
                end
           end
        end
        
    end
end 

output = c;
