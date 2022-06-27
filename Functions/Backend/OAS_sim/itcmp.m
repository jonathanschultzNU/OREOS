function [params, M, itc, s] = itcmp(y, M)
%ITCMP Information theoretic criteria and matrix pencil method
%	Fits a time series of (complex) data to a sum of complex exponentials
%	using the matrix pencil method.  The number of underlying signals can
%	be supplied directly or estimated based on information theoretic
%	criteria.  The returned model parameters describe the amplitudes,
%	frequencies, damping rates, and phases of the fitted exponentials.
%
% [params, M, itc, s] = itcmp(y, M)
%
% INPUT
%	y		(numeric) Time series vector of (complex) data points.
%	M		(real) Number of signal components or effective matrix rank.
%			For automatic selection based on estimated information content,
%			use M = -1 (AIC criterion) or M = -2 (MDL criterion).  For 
%			M >= 0, the user-defined value is used directly.
%
% OUTPUT
%	params	(real) M*4 matrix with columns for damping rate, frequency
%			amplitude, and phase
%	M		(real) Number of signal components used in the fit.
%	itc		(numeric) Vector of AIC or MDL function values.
%	s		(numeric) Complex frequencies
%
% NOTE
%	This function is adapted from work by Yung-Ya Lin (dated 4/15/96) 
%	published in "Lin, Y.-Y., Hodgkinson, P., Ernst, M. & Pines, A. A novel 
%	detection–estimation scheme for noisy NMR signals: applications to 
%	delayed acquisition data. Journal of Magnetic Resonance 128, 30–41 
%	(1997).  https://doi.org/10.1006/jmre.1997.1215"
%
%	See also itcmpEval, itcmpFilterOscillations, itcmpParamConversion,
%	svdecon.
%

%	Author: Austin P. Spencer
%	Email: austin.spencer (at) northwestern.edu
%	Last revision date: April 5, 2019
%
%	Copyright: Austin P. Spencer, 2019

	y = double(y(:));
	N = length(y);
	L = floor(N/3);						% pencil parameter
	Y = toeplitz(y(L+1:N),y(L+1:-1:1)); % YO = Y(:,2:L+1), Y1 = Y(:,1:L) Eq. [3]
	[U, S, V] = svdecon(Y(:,2:L+1));	% singular value decomposition
	S = diag(S);
	
	% Use information theoretic criterion to determine number of components
	itc = zeros(1,L);
	switch M
		case -1							% determining M by AIC
			for k = 0:L-1
				itc(k+1) = -2*N*sum(log(S(k+1:L))) ...
					+ 2*N*(L-k)*log((sum(S(k+1:L))/(L-k))) + 2*k*(2*L-k);
			end
			[~, tempI] = min(itc); M = tempI-1;
		case -2							% determining M by MDL Eq. [16]
			for k = 0:L-1
				itc(k+1) = -N*sum(log(S(k+1:L))) ...
					+ N*(L-k)*log((sum(S(k+1:L))/(L-k))) + k*(2*L-k)*log(N)/2;
			end
			[~, tempI] = min(itc); M = tempI-1;
		otherwise
			if M >= 0					% user-defined M
				
			else						% Invalid value for M.
				error('Invalid value of M: must be >= -2.');
			end
	end
	
	% Attempt to handle the case where M is set higher than the data can
	% support.  These cases otherwise result in erroneous fits.
	if M>(L/2)
		if any(itc)	% If using information theoretic criterion...
			% Truncate M at the point where itc goes from decreasing to
			% increasing.
% 			warning('M>(L/2). Truncating M to first local minimum of itc.');
			tempI = find(diff(itc)>0,1); M = tempI-1;
		else
			% No correction possible when not using information theortic
			% criterion.
% 			warning('M>(L/2). Fit may be unreliable.');
		end
	end
	
	s = log(eig(diag(1./S(1:M))*((U(:,1:M)'*Y(:,1:L))*V(:,1:M))));
	
	% signal pole z = exp(s)
	Z = zeros(N,M);
	for k = 1:M
		Z(:,k) = exp(s(k)).^(0:N-1).';
	end
	
	% Remove infinite values.
	s = s(~any(isinf(Z)));
	Z = Z(:,~any(isinf(Z)));
	
	% linear least squares analysis
	a = Z\y;
	
	% Extract component fit parameters
	params = [-real(s) imag(s)/2/pi abs(a) imag(log(a./abs(a)))];
end