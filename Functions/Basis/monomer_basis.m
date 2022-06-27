function output = vibron_mon_mv(E,numvib,vmax)

length = sum(vmax)+1;
basis = zeros(length,1+numvib);

k = 1;
for i=1:numvib
    
    if i==1
        for v = 0:vmax(i)
            basis(k,1) = E;
            basis(k,i+1) = v;
            k = k+1;
        end
    end
    
    if i>1
        for v = 1:vmax(i)
            basis(k,1) = E;
            basis(k,i+1) = v;
            k = k+1;
        end
    end
    
end

output = basis;

end


