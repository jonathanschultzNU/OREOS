function p = findind(vect,val)
    [~,p] = min(abs(vect-val));
end