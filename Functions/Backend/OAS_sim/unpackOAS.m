function dataw1w3fits = unpackOAS(mat,num_w3,num_w1,num_t2)

dataw1w3fits = zeros(num_w3,num_w1,num_t2);

k = 1;
for j=1:num_w1
    for i=1:num_w3
        dataw1w3fits(i,j,:) = mat(k,:);
        k=k+1;
    end
end

end