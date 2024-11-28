% Generate surrogate data with matching amplitude spectrum and
% amplitude distribution (Schreiber and Schmitz, 1996).
% Multivariate data as described in Schreiber and Schmitz (2000)
% Leakage effects have been compensated for (end2end) for 1D cases
%
% Usage: [Xs X] = generate_surrogate (X, endmode, specflag);
%	endmode=1	use end2end compensation
%	specflag	exact amplitude spectrum (1, default), otherise amp distr
%	X [pp x dim]

function [X2,X,E] = generate_surrogate (X, endmode, flag);


if (nargin<2)
	endmode = 1;
end
if (nargin<3)
	flag = 1;
end

max_it = 500;

		%%%%%%%%%%%%%%%%
		% End Matching %
		%%%%%%%%%%%%%%%%
[pp dim] = size(X);
if (pp==1)
	X = X';
	[pp dim] = size(X);
end
if ((dim==1) & (endmode) & (abs(X(1)-X(pp))>1e-3))
	l = floor(pp/25)+1;
	%l = floor(pp/10)+1;
	T1 = X(1:l);
	T2 = X(pp-l+1:pp);
	D_2 = (repmat(T1(:),1,l)-repmat(T2(:)',l,1)).^2;
	%D_2(find(eye(l))) = mmax(D_2);
	[v ind1 ind2] = minmin(D_2);
	X = X(ind1:pp-l+ind2);
	if (v>1e-1)
		%fprintf ('Unable to make ends meet\n');
		X2 = [];
		X = [];
	end
end

if (~isempty(X) & (dim==1))
	X = X(:);
	pp = length(X);

	Yamp = abs(fft(X));
	if (rem(pp,2)==0)
		l = pp/2-1;
		T1 = rand(1,l).*(2*pi)-pi;
		Yang = [0 T1 0 -T1(l:-1:1)];
	else
		l = (pp-1)/2;
		T1 = rand(1,l).*(2*pi)-pi;
		Yang = [0 T1 -T1(l:-1:1)];
	end
	%[N1 INDh] = hist(X, 50);

	% Initial Conditions
	rn = X(randperm(pp));
	Xs = sort(X);

	prev_err = 0;
	E = zeros(1,max_it);
	c = 1;
	prev_err = 1000000;
	err = prev_err - 1;
	while (prev_err>err) & (c<max_it)

		% Match Amplitude Spec
		Yrn = fft(rn);
		Yang = angle(Yrn);
		sn = real(ifft(Yamp.*(cos(Yang)+sqrt(-1).*sin(Yang))));

		% Scale to Original Amp Distr
		[sns INDs] = sort(sn);
		rn(INDs) = Xs;

		% Eval Convergence
			%[N2 xx] = hist(sn, INDh);
			%prev_err = err;
			%err = sum(abs(N1-N2))/pp;
		prev_err = err;
		A2 = abs(Yrn);
		%err = mean(mean(abs(A2-Yamp)));
		err = mean(abs(A2-Yamp));
		E(c) = err;

		c = c+1;
		%plot (E); pause
	end
	E = E(1:c-1);
	if (flag==1)
		X2 = sn;	% Exact Amp Spectrum
	else
		X2 = rn;	% Exact Amp Distribution
	end

	if (endmode==2)
		X2 = X2./h;
	end


end

if (dim>1)
	%Gref = randn(pp,dim);
	%Grefs = sort(Gref);
	Y = fft(X);		% fft every column
	Yamp = abs(Y);
	Porig = angle(Y);
	%[N1 INDh] = hist(X, 50);

	% Initial Conditions
	rn = zeros(size(X));
	for k=1:dim
		rn(:,k) = X(randperm(pp),k);
	end
	Xs = sort(X);

	prev_err = 1000000;
	err = prev_err - 1;
	%E = zeros(1,n_it); E(n_it)=1;
	c = 1;
	Pcurr = Porig;
	%while ((E(n_it)-E(n_it-1)>0) & (c<10))
	while (prev_err>err) & (c<max_it)
		% Match Amplitude Spec
		Prn = angle(fft(rn));
		goal = Prn - Porig;
		AUX1 = sum(cos(goal),2);
		alpha = repmat((AUX1~=0).*atan(sum(sin(goal),2)./sum(AUX1+(AUX1==0),2)),1,dim);
		alpha = alpha + repmat(pi.*(sum(cos(alpha-goal),2)<0),1,dim);
		Pcurr = Porig + alpha;
		sn = real(ifft(Yamp.*(cos(Pcurr)+sqrt(-1).*sin(Pcurr))));

		% Scale to Gaussian
		[sns INDs] = sort(sn);
		for k=1:dim
			rn(INDs(:,k),k) = Xs(:,k);
		end

		% Eval Convergence
			%[N2 xx] = hist(sn, INDh);
			%prev_err = err;
			%err = sum(sum(abs(N1-N2)))/pp/dim;
		prev_err = err;
		A2 = abs(fft(rn));
		err = mean(mean(abs(A2-Yamp)));
		E(c) = err;
		c = c+1;
	end
	%plot(E(1:c));drawnow
	if (flag==1)
		X2 = sn;	% Exact Amp Spectrum
	else
		X2 = rn;	% Exact Amp Distribution
	end

end
