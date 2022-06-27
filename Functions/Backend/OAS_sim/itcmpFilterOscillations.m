function [nonOsc, osc] = itcmpFilterOscillations( params, n, varargin)
%ITCMPFILTEROSCILLATIONS Separate components of ITCMP fit by frequency
%	Categorizes and separates itcmp fit components as "oscillatory" or
%	"non-oscillatory" based on the ratio of damping rate to frequency.
%
% [ nonOsc, osc ] = itcmpFilterOscillations(params, n, [NAME, VALUE])
%
% INPUT
%	params			(real) Output parameters from itcmp with columns for
%					damping rate, frequency, amplitude, and phase
%	n				(real) number of time series points in the data 
%					modelled by itcmp
%	*oscThreshold	(real) Threshold for defining which components are 
%					oscillatory.  If the exponential damping rate is less
%					than oscThresh times the frequency, it is considered an
%					oscillation.
%					Default: 5
%
% * indicates optional name-value ('key', value) input parameter
%
% OUTPUT
%	nonOsc	(real) Matrix of non-oscillatory itcmp components with
%			columns for damping rate, frequency, amplitude, and phase.
%	osc		(real) Matrix of oscillatory itcmp components with columns
%			for damping rate, frequency, amplitude, and phase.
%
%	See also itcmp, itcmpEval, itcmpParamConversion.
%

%	Author: Austin P. Spencer
%	Email: austin.spencer (at) northwestern.edu
%	Last revision date: April 5, 2019
%
%	Copyright: Austin P. Spencer, 2019

c = 2.9979 * 10^-5;       %speed of light[cm/fs]

	% Parse function inputs
	ip = inputParser;
	ip.addRequired('params', @isreal);
	ip.addRequired('n', @isreal);
	ip.addParameter('oscThresh', 5, @isreal);
	ip.parse(params, n, varargin{:});
	inputs = ip.Results;

	% Functional form: A*exp(-rate*t + 1i*2*pi*f*t + 1i*phi)
	A = params(:,3);		% Amplitude
	rate = params(:,1);		% Damping rate
	f = params(:,2);		% Frequency
    fwav = f./c;
	phi = params(:,4);		% Phase
	
	% Determine which components are oscillatory/non-oscillatory
	isOsc = false(1, size(params, 1));
%     islowOsc = false(1, size(params, 1));
%     ishighOsc = false(1, size(params, 1));
	for ii=1:size(params, 1)
        
        %FREQUENCY CUTOFF METHOD
%         if abs(fwav(ii)) > abs(10)
% 				isOsc(ii) = true;
%                 if abs(fwav(ii)) < abs(150)
%                     islowOsc(ii) = true;
%                 elseif abs(fwav(ii)) >= abs(150)
%                     ishighOsc(ii) = true;
%                 end
%         end
        
        %OSCTHRESH METHOD
		% Regardless of damping rate, if the oscillation period is long
		% relative to the time range of the data, treat the component as
		% non-oscillatory.  In these cases, the damping rate is poorly
		% constrained and effectively meaningless.
		if abs(f(ii)) > (1/(2*n))	% period < 2*range(time) ?
			% If oscThresh*f > rate, classify component as oscillatory.
			% In words, if the component lifetime is more than 1/oscThresh
			% (1/5 by default) of a period, it is oscillatory.
            if inputs.oscThresh*abs(f(ii)) > abs(rate(ii))
				isOsc(ii) = true;
            end
           
		end

	end
	
	% Separate components into oscillatory and non-oscillatory
	nonOsc = params(~isOsc,:);
	osc = params(isOsc,:);
%     highosc = params(ishighOsc);
%     lowosc = params(islowOsc,:);

end

