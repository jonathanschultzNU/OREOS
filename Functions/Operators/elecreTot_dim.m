function output = elecregenTot_v1(basis,n_states,molnum,numvib,quanta)

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
elseif quanta == 2
    lownum = 1; 
    conjnum = 0;
elseif quanta == 5
    lownum = 4;
    conjnum = 0;
elseif quanta == 6
    lownum = 1;
    conjnum = quanta-1;
end

for i = 1:n_states
    for j = 1:n_states
        bra = basis(i,:);
        ket = basis(j,:);
        
        if ket(index) == quanta
           ket(index) = ket(index)-lownum; 
           ket(conjindex) = conjnum;

                    tf = isequal(ket,bra);
                    if tf>0   
                        c(i,j) = 1;       
                    end
            
        end
    end
end

output = c;
end
