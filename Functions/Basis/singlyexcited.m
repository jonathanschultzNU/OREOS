function output = singlyexcited(numvib,vmax)

% length = 2*(vmax(1)^2 + 3*vmax(1)+1);
% basis = zeros(length,2+2*numvib);

k = 1;

for i=1:numvib
    
for n = 0:1    %index over ground and singly excited state
    m = not(n);
    for v_1 = 0:vmax(i)
        for v_2 = 0:vmax()-v_1
        
        basis(k,1) = n;
        basis(k,2+numvib) = m;
        basis(k,i+1) = v_1;
        basis(k,i+2+numvib) = v_2;
        k = k+1;
        if i>1
            if v_1==0 && v_2 == 0
                k = k-1;
            end
        end
        end
    end
end

end

output = basis;


