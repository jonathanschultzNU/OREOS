function [phi] = PhaseMean(phases)
%PhaseMean finds the average phase in a region
% It is important to avoid a 2 pi shift when calculating the average phase.
% Here we check if the phases in a region are close to pi or - pi, and if 
% so we shift the negative phase values by 2 pi so that all of the phases 
% are positive, avoiding a 2 pi shift.
% 
% [phi] = PhaseMean(phases)
%
% INPUT
%   phases        An array of phases
%
% OUTPUT
%   phi           Average phase

%   Author: Matthew S. Kirschner
%   Email: kirschner.21 (at) gmail.com
%   Last revision date: August 6, 2019
%
%   Copyright: Matthew S. Kirschner, 2019


if mean(abs(phases)) > (3 * pi / 4)
    phases(phases < 0) = phases(phases < 0) + 2 * pi;
end

% after that shift, we just find the mean value
phi = mean(phases);

end