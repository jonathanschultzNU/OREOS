function Rw1w3 = rescalc(Rt1t3,Res,i)

Rw1w3 = zeros(Res.nw3,Res.nw1,length(Res.t2));

for k = 1:length(Res.t2)
    
    switch i
        case 1
            temp = fft2(Rt1t3(:,:,k),Res.nw3,Res.nw1);   
            temp = rot90(temp,2);
            Rw1w3(:,:,k) = fftshift(temp);
        case 2
            temp = fft2(Rt1t3(:,:,k),Res.nw3,Res.nw1); 
            temp = fliplr(circshift(temp,-1,2)); 
            temp = rot90(temp,2);
            Rw1w3(:,:,k) = fftshift(temp);
        case 3
            temp = fft2(Rt1t3(:,:,k),Res.nw3,Res.nw1); 
            temp = fliplr(circshift(temp,-1,2));
            temp = rot90(temp,2);
            Rw1w3(:,:,k) = fftshift(temp);    
        case 4
            temp = fft2(Rt1t3(:,:,k),Res.nw3,Res.nw1);  
            temp = rot90(temp,2);
            Rw1w3(:,:,k) = fftshift(temp);
        case 5
            temp = fft2(Rt1t3(:,:,k),Res.nw3,Res.nw1); 
            temp = fliplr(circshift(temp,-1,2)); 
            temp = rot90(temp,2);
            Rw1w3(:,:,k) = -fftshift(temp);
        case 6
            temp = fft2(Rt1t3(:,:,k),Res.nw3,Res.nw1);  
            temp = rot90(temp,2);
            Rw1w3(:,:,k) = -fftshift(temp);
        otherwise
            disp('Invalid response function specification')   
    end 
    clear temp
end