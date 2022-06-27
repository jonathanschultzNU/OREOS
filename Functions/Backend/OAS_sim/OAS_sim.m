function [ExponentialFits,OscFits] = OAS_sim(datain)

time = datain(1,2:end);
wavelength = datain(2:end,1);
data = datain(2:end,2:end);

first_wavelength = min(wavelength); % the first wavelength in the OAS
last_wavelength = max(wavelength); % the final wavelength you want in your OAS
first_time = min(time); % the first time point for your Fourier Transform
last_time = max(time); % the final time point for your fit
time_resolution = mean(diff(time)); % the time resolution of your measurements, your data
% will be resampled/interpolated at this rate
criterion = 10; %-1; %-2; %Information criterion for itcmp, -1 uses AIC, -2 MDL 
% see itcmp for more info
oscillation_ratio = 5; %75; % Threshold for defining which components are 
% oscillatory.  If the exponential damping rate is less than oscThresh 
% times the frequency, it is considered an oscillation. See 
% itcmpFilterOscillations for more details  (testing = wavenumber units)
c = 2.9979 * 10^-5;       %speed of light[cm/fs]
lowwavenumber = 10;
highwavenumber = 50;
initial_period  = 1./(highwavenumber*c); % the shortest period you are interested in
final_period = 1./(lowwavenumber*c); % the longest period you are interested in
% res_wavenumber = 3;
resolution = 8; %1./(res_wavenumber*c); % the resolution you want

% now we convert your input values to a matrix
periods = initial_period : resolution : final_period;

% The actual OAS calculation
progressbar2 % Init single bar
[ExponentialFits,OscFits] = OscillationAssociatedSpectrum(data, ...
    wavelength, time, first_wavelength, last_wavelength, first_time,...
    last_time, time_resolution, periods, criterion, oscillation_ratio);

progressbar2(1)            % Close

end
