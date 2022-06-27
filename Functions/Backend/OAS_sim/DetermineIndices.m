function indices = DetermineIndices(in_array, first, last)
%DetermineIndices find the indices that correspond to a range in an array
%   Identify firt and last indices of an array within the given value range
%   indices = DetermineIndices(in_array, first, last)
% 
% indices = DetermineIndices(in_array, first, last)
%
% INPUT
%   in_array      The data array you are looking for the indices of
%   first         First value of interest
%   last          Final value of interest
%
% OUTPUT
%   indices       First and last indices of the in_array that correspond to
%                 the range of interest


%   Author: Matthew S. Kirschner
%   Email: kirschner.21 (at) gmail.com
%   Last revision date: August 6, 2019
%
%   Copyright: Matthew S. Kirschner, 2019


% handle if the input array is not sorted
[sorted_array, indices] = sort(in_array);

% finds indices in an array
fi = find(sorted_array > first, 1); %find lower bound
li = find(sorted_array > last, 1) - 1; %find upper bound

% assign an output value if the last value is too high
if isempty(li)
    li = length(sorted_array);
end

% get the output values
indices = [indices(fi), indices(li)];

end

