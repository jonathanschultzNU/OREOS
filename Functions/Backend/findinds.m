function results = findinds(vect,minval,maxval)
    [~,pmin] = min(abs(vect-minval));
    [~,pmax] = min(abs(vect-maxval));
    results = [pmin,pmax];
end