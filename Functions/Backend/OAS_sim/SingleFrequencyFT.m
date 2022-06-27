function [amp, ang] = SingleFrequencyFT(x, y, period)
%SingleFrequencyFT Calculate the discrete Fourier transform at an
%individual wavelength
% 
% [amp, ang] = SingleFrequencyFT(x, y, period)
%
% INPUT
%   x             Time data (this function assumes an even sampling rate)
%   y             Amplitude data
%   period        Oscillation period of interest
%
% OUTPUT
%   amp           Magnitude of Fourier Transform
%   ang           Phase of Fourier Transform


%   Author: Matthew S. Kirschner
%   Email: kirschner.21 (at) gmail.com
%   Last revision date: August 6, 2019
%
%   Copyright: Matthew S. Kirschner, 2019

% sine integral
a = sum(sin(2 * pi * x / period) .* y);

%cosine integral
b = sum(cos(2 * pi * x / period) .* y);

% FT magnitude
amp = sqrt(a^2 + b^2);

% FT phase
ang = atan2(a, b);

end

