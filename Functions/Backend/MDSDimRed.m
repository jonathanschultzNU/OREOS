function datadimred = MDSDimRed(data_w,w1_rel,w3_rel,t2)

%Converting data to 2D matrix

num_w1 = length(w1_rel);
num_w3 = length(w3_rel);
num_t2 = length(t2);
h = num_w1*num_w3;

data2D = zeros(num_t2,h);
k = 1;
for j=1:num_w1
    for i=1:num_w3
        data2D(:,k) = data_w(i,j,:);
        k=k+1;
    end
end

data2D = data2D';

%bFormat 2D matrix for TASVD

wavelength_ind = linspace(1,h,h);
wavelength_ind = wavelength_ind';

data_TASVD = zeros(h+1,length(t2)+1);
data_TASVD(2:end,1) = wavelength_ind;
data_TASVD(1,2:end) = t2;

data_TASVD(2:end,2:end) = data2D;

datadimred = data_TASVD;

end