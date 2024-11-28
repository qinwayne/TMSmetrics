% Generate a Henon Map
%
% Usage: X = generate_henon (pp, tau)

function x = generate_henon (pp, tau)

if (nargin<1)
	pp = 1000;
end
if (nargin<2)
	tau = 1;
end

a = 1.4;
b = 0.3;
x = randn(pp,1).*.01;
y = randn(pp,1).*.01;

for n = tau+1:pp,
	x(n) = 1 - a*x(n-tau)*x(n-tau) + y(n-tau);
	y(n) = b*x(n-tau);
end

