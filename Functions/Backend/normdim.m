function data = normdim(data)

d = ndims(data);

if d == 1
    data = data./max(abs(real(data)));
end

if d == 2
    data = data./max(max((abs(real(data)))));
end

if d == 3
    data = data./max(max(max(abs(real(data)))));
end

end