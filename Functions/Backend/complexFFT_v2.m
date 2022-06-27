function output = complexFFT_v2(beat_iso_real,beat_iso_imag,t2,w1,w3,npad,n_cut_zero)

num_w1 = length(w1);
num_w3 = length(w3);

beat_iso_complex = beat_iso_real+1i.*beat_iso_imag;
han2 = 0.5;
win_2 = hanning(2 .* length(t2)) .^ han2;       %create the symmetric Hann window
win_2 = win_2((length(win_2) / 2) + 1 : end);     %cut the Hann Window in half to match FID form

for i = 1:num_w3
    for j = 1:num_w1
        
        data_beat_w = zeros(1,npad);                  %initialize Fourier transformed data matrix
        q = zeros(1,npad);                            %initialize q
        tempvect(1,:)=beat_iso_complex(i,j,:);
        q(n_cut_zero+1:n_cut_zero+length(t2)) = tempvect.*win_2';   %load beat into q r = r.*win_lpsvd';
        q(n_cut_zero+1) = q(n_cut_zero+1)./2;                                 %divide initial value by two to compensate trapezoidal rule
        data_beat_w = fftshift(fft(q));                 %Fourier transform and load into initialized/padded vector
        data_beat_freq(i,j,:) = data_beat_w;
    end
end

output.data_w1w2w3 = data_beat_freq;