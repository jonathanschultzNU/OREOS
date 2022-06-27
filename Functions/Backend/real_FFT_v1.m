function dataw1w2w3 = real_FFT_v1(beatiso,t2,w1,w3,npad,nzero)

num_w1 = length(w1);
num_w3 = length(w3);

han2 = 0.5;
win_2 = hanning(2 .* length(t2)) .^ han2;       %create the symmetric Hann window
win_2 = win_2((length(win_2) / 2) + 1 : end);     %cut the Hann Window in half to match FID form

dataw1w2w3 = zeros(num_w3,num_w1,npad);

for i = 1:num_w3
    for j = 1:num_w1
        
        q = zeros(1,npad);                            %initialize q
        tempvect(1,:)=beatiso(i,j,:);
        q(nzero+1:nzero+length(t2)) = tempvect.*win_2';   %load beat into q r = r.*win_lpsvd';
        q(nzero+1) = q(nzero+1)./2;                                 %divide initial value by two to compensate trapezoidal rule
        data_beat_w = fftshift(fft(q));                 %Fourier transform and load into initialized/padded vector
        dataw1w2w3(i,j,:) = data_beat_w;
        
    end
end

end