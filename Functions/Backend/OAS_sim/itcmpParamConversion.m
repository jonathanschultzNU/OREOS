function [paramsStruct, paramsTable] = itcmpParamConversion(params, dt, varargin)
%ITCMPPARAMCONVERSION Convert itcmp parameters to physical units
%   Convert the raw "params" output by itcmp into values with physically
%   meaningful units.  Converted parameters are returned as a structure
%   and/or a table and are optionally printed to the command window.
%	
% [paramsStruct, paramsTable] = itcmpParamConversion(params, dt, [silent], 
%									[NAME, VALUE])
%
% INPUT
%	params			(real) Output parameters from itcmp with columns for
%					damping rate, frequency, amplitude, and phase
%	dt				(real) Step size of time axis (in seconds) for data fit
%					by itcmp
%	silent			(logical, optional) Flag specifying whether to print
%					output to command window or not
%					Default: false
%	*frequencyUnits	(char) Frequency units to convert to
%					Options: '1/cm', 'Hz'
%					Default: '1/cm'
%
% * indicates optional name-value ('key', value) input parameter
% [] indicates optional input
%
% OUTPUT
%	paramsStruct	(struct) Structure containing unit-converted params.
%					Fields:	tau (lifetime), f (frequency), A (amplitude), 
%							phi (phase)
%	paramsTable 	(table) Table containing unit-converted params
%
%	See also itcmp, itcmpEval, itcmpFilterOscillations.
%

%	Author: Austin P. Spencer
%	Email: austin.spencer (at) northwestern.edu
%	Last revision date: April 5, 2019
%
%	Copyright: Austin P. Spencer, 2019

	% Parse function inputs
	ip = inputParser;
	ip.addRequired('params', @isreal);
	ip.addRequired('dt', @isreal);
	ip.addOptional('silent', false, @islogical);
	ip.addParameter('frequencyUnits', '1/cm', ...
		@(x) any(validatestring(x, ...
		{'wavenumbers', '1/cm', 'cm^-1', 'Hz', '1/s'})));
	ip.parse(params, dt, varargin{:});
	params = ip.Results.params;
	dt = ip.Results.dt;
	silent = ip.Results.silent;
	frequencyUnits = ip.Results.frequencyUnits;
	
	% Select scaling factor and name for requested frequency units
	switch frequencyUnits
		case {'wavenumbers', '1/cm', 'cm^-1'}
			freqUnitFactor = 1/(c*100*dt);
			freqUnitName = '1/cm';
		case {'Hz', '1/s'}
			freqUnitFactor = 1/dt;
			freqUnitName = 'Hz';
	end
	
	% Convert time and frequency units.  Sort components by frequency.
	params(:,1) = dt./params(:,1);
	params(:,2) = params(:,2)*freqUnitFactor;
	params = sortrows(params, -2);
	
	% Generate outputs
	if nargout==1 && silent
		paramsStruct = struct('tau', params(:,1)', ...
			'f', params(:,2)', ...
			'A', params(:,3)', ...
			'phi', params(:,4)');
	else
		paramsTable = array2table(params, 'VariableNames', ...
			{'tau', 'f', 'A', 'phi'});
		paramsTable.Properties.VariableUnits = {'s', freqUnitName, '', 'radian'};
		paramsStruct = table2struct(paramsTable);
		if ~silent
			disp(paramsTable);
		end
	end
	
end

