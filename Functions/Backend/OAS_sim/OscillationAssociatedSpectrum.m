function [ExponentialFits,OscFits] = OscillationAssociatedSpectrum(data, ...
    wavelength, time, first_wavelength, last_wavelength, first_time,...
    last_time, time_resolution, periods, criterion, oscillation_ratio)
%OscillationAssociatedSpectrum Calculates OAS of data
% 
% [Amplitude, Phase, w] = OscillationAssociatedSpectrum(data, ...
%    wavelength, time, first_wavelength, last_wavelength, first_time,...
%    last_time, time_resolution, periods, criterion, oscillation_ratio)
%
% INPUT
%   data                2D-array (wavelength x time) of transient 
%                       absorption (TA) data
%   wavelength          1D-array corresponding to the wavelengths in 
%                       the TA data
%   time                1D-array corresponding to the times in the TA data 
%   first_wavelength    First wavelength in the OAS
%   last_wavelength     Final wavelength in the OAS
%   first_time          First time to be used in the Fourier transform
%   last_time           Final time to be used in the Fourier transform
%   time_resolution     Time dependent data will be resampled/interpolated
%                       at this rate
%   periods             Oscillation periods to be considered
%   crieterion          Information criterion for itcmp, -1 uses AIC, 
%                       -2 MDL see itcmp for more info
%   oscillation_ratio   Threshold for defining which components are 
%                       oscillatory.  If the exponential damping rate is
%                       less than oscThresh times the frequency, it is 
%                       considered an oscillation. See 
%                       itcmpFilterOscillations for more details
%
% OUTPUT
%   Amplitude           2D-amplitude map of the OAS at every
%                       wavelength/period
%   Phase               2D-phase map of the OAS at every period
%   w                   1D-array of wavelengths used in the OAS


%   Author: Matthew S. Kirschner
%   Email: kirschner.21 (at) gmail.com
%   Last revision date: August 6, 2019
%
%   Copyright: Matthew S. Kirschner, 2019


% find the wavelength/time range of interest
wavelength_indices = sort(findinds(wavelength, first_wavelength,...
    last_wavelength));
first_wavelength_index = wavelength_indices(1); 
last_wavelength_index= wavelength_indices(2);

time_indices = sort(findinds(time, first_time, last_time));
first_time_index = time_indices(1);
last_time_index = time_indices(2);


% initialize amplitude and phase matrices
Amplitude = zeros(last_wavelength_index - first_wavelength_index + 1,...
    length(periods));

Phase = zeros(last_wavelength_index - first_wavelength_index + 1,...
    length(periods));

ExponentialFits = zeros(last_wavelength_index-first_wavelength_index+1,last_time_index-first_time_index+1);
OscFits = zeros(last_wavelength_index-first_wavelength_index+1,last_time_index-first_time_index+1);
% LowOscFits = zeros(last_wavelength_index-first_wavelength_index+1,last_time_index-first_time_index+1);

% the old time values
old_time = time(first_time_index : last_time_index);

% the new time values that the data will be interpolated to
% new_time = old_time(1) : time_resolution : old_time(end);
new_time = old_time;

count = 1;
% the actual fit, performed at every wavelength independently
k = 1; 
% figure
for current_wavelength = first_wavelength_index : last_wavelength_index
    
    % resample the data
%     resampled_data = interp1(old_time, data(current_wavelength, ...
%         first_time_index : last_time_index), new_time);
    resampled_data = data(current_wavelength,first_time_index : last_time_index);
    
    % set NaN values to 0
    resampled_data(~isfinite(resampled_data)) = 0;
   
    % matrix pencil method
    temppara = itcmp(resampled_data, criterion);
    
    % figure out which decays are complex
    signalFitExp = itcmpEval(temppara, length(new_time), 'components',...
        'nonOsc','OscThresh',oscillation_ratio); %change here if want to output oscillatory fits
    
    ExponentialFits(count,:) = real(signalFitExp);
%     subplot(3,1,1)
%     plot(new_time,ExponentialFits(count,:))
%     hold on
%     plot(new_time,resampled_data)
%     hold off
    
    signalFitOsc = itcmpEval(temppara, length(new_time), 'components',...
        'osc','OscThresh',oscillation_ratio); %change here if want to output oscillatory fits

    OscFits(count,:) = real(signalFitOsc);
%     subplot(3,1,2)
%     plot(new_time,OscFits(count,:))
%     hold on
%     plot(new_time,resampled_data)
%     hold off
    
%     subplot(3,1,3)
%     plot(new_time,ExponentialFits(count,:)+OscFits(count,:))
%     hold on
%     plot(new_time,resampled_data)
%     hold off
%     signalFitLowOsc = itcmpEval(temppara, length(new_time), 'components',...
%         'lowosc','OscThresh',oscillation_ratio); %change here if want to output oscillatory fits
% %     
%     LowOscFits(count,:) = real(signalFitLowOsc);
%     subplot(1,3,3)
%     plot(time,LowOscFits(count,:))
    
    count = count+1;
    
    % perform the Fourier transformations
%     for period = 1 : length(periods)
%         [Amplitude(current_wavelength + 1 - first_wavelength_index, period),...
%             Phase(current_wavelength + 1 - first_wavelength_index, period)]...
%             = SingleFrequencyFT(new_time, resampled_data - real(signalFitExp),...
%             periods(period));
%     end

if ismember(k,(0:500:100000)) == 1    
        progressbar2(k/(length(wavelength))) % Update progress bar
end

k=k+1;
    
end

% extract the wavelengths used
w = wavelength(first_wavelength_index : last_wavelength_index);

end