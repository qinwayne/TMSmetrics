% clear all
function [best_m,best_tau] = optim_m_tau(data)

rand ('seed', 1); randn ('seed', 1);

m1 = 2;		% Minimal m
m2 = 10;	% Maximal m
T = 1:10;	% Values for tau
n_surr = 5;	% Number of surrogates

s = 1;		% Resampling factor

		%%%%%%%%%%%%%
		% Read Data %
		%%%%%%%%%%%%%
X = data;
pp=500; 


		%%%%%%%%%%%%%%%%%%%%%%%%%%%
		% Standardise Time Series %
		%%%%%%%%%%%%%%%%%%%%%%%%%%%
X = X-mean(X);
X = X./std(X);
X = X(:);
pp = length(X);


		%%%%%%%%%%%%%%
		% Resampling %
		%%%%%%%%%%%%%%
if (s~=1)
	X = interp1(0:1/(pp-1):1,X,0:1/(s*pp-1):1,'linear');
end


		%%%%%%%%%%%%%%%%%%%%%%%%
		% Differential Entropy %
		%%%%%%%%%%%%%%%%%%%%%%%%
[H,M,T] = sweep_kl (X, m1, m2, T);
for s=1:n_surr
	Xs = generate_surrogate(X,0,0);
	Hs(:,:,s) = sweep_kl (Xs, m1, m2, T);
end
D = (H./mean(Hs,3))';
p1 = pp - max(T)*m2;
D =  D+log(p1)*repmat(m1:m2,length(T),1)./p1;	% MDL


		%%%%%%%%%%%%%%%%
		% Find Minimum %
		%%%%%%%%%%%%%%%%
[a b c] = minmin(D);
best_m = m1+c-1;
best_tau = T(b);
fprintf ('Min at m = %d, tau = %d (%.4f)\n', m1+c-1, T(b), D(b,c));

if (m1~=m2)
	Hp = mesh (M,T,D);
	set(Hp, 'EdgeColor','k');
	hold on
	plot3 (m1+c-1,T(b), D(b,c), 'or')
	hold off
	xlabel ('m')
	ylabel ('\tau');
	zlabel ('R_{ent}')
else
	plot (T,D);
	xlabel ('\tau');
end

