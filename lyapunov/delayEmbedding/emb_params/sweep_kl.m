% Compute the Kosachenko-Leonenko differential entropy
% of phase spaces of increasing embedding dimension [m_min m_max]
% and a range of time lags TAU
%
% Usage: [H,M,TAU] = sweep_kl (X,m_min,m_max,TAU);

function [H,M,TAU] = sweep_kl (X,m_min,m_max,TAU);

if (nargin<1)
	error ('Provide at least an input X');
end
if (nargin<2)
	m_min = 2;
end
if (nargin<3)
	m_max = 10;
end
if (nargin<4)
	TAU = 1:10;
end

M = m_min:m_max;
H = sweep_kl_mex (X,m_min,m_max,TAU);


