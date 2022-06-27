function output = elcregen_v1(basis,n_states,molnum,nvib,numvib,quanta)

index = (numvib+1)*(molnum-1)+1;

if molnum == 1
    conjindex = (numvib+1)*(2-1)+1;
elseif molnum == 2
    conjindex = (numvib+1)*(1-1)+1;
end

c = zeros(n_states);

if quanta == 1
    lownum = 1;     %separation between S0 and S1 is 1
    conjnum = 0;
    opt = 1;
elseif quanta == 5
    lownum = 4;
    conjnum = 0;
    opt = 1;
else
    opt = 0;
end

for i = 1:n_states
    for j = 1:n_states
        bra = basis(i,:);
        ket = basis(j,:);
        
        if ket(index) == quanta
           ket(index) = ket(index)-lownum; 
           ket(conjindex) = 0;
           
           if molnum == 1
           
               test = ket; 
               test(1) = [];
               test(numvib+1) = [];
               
               if test == zeros(size(test))
                    tf = isequal(ket,bra);
                    if tf>0   
                        c(i,j) = 1;       
                    end

               elseif ket(index+nvib) > 0 || ket(index+nvib+numvib+1) > 0  
                    tf = isequal(ket,bra);
                    if tf>0   
                        c(i,j) = 1;       
                    end
                   
               end
               
           elseif molnum == 2
               
               test = ket; 
               test(1) = [];
               test(numvib+1) = [];
               
               if test == zeros(size(test))
                    tf = isequal(ket,bra);
                    if tf>0   
                        c(i,j) = 1;       
                    end

               elseif ket(index+nvib) > 0 || ket(1+nvib) > 0  
                    tf = isequal(ket,bra);
                    if tf>0   
                        c(i,j) = 1;       
                    end
               end
           end
           
        else 
            
        end
    end
end

output = c;
end
