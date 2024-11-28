function [LyaExp_b_sec,LyaExp_b_orb] = lyapunov_expEst(out)
close all;
% Mean Orbital Period from PSD
domFreq = 1;
meanPeriod = 1/domFreq;
% Output data fetout.txt from lyapunov.m code
% bitsPerSec_Data = importdata('fetout.txt'); 
% Estimate of lyapunov exponent for each increment
bitsPerSec_Exp = out(:,4); 
i = 1;
ind = 1;
while i <= length(bitsPerSec_Exp)
 if isnan(bitsPerSec_Exp(i)) ~= 1
 BPS(ind) = bitsPerSec_Exp(i);
 ind = ind + 1;
 end
 i = i + 1;
end
BPS = BPS';
bitsPerOrbit = BPS.*meanPeriod;
% plot(bitsPerOrbit,'linewidth',2)
% hold on; plot([1 length(BPS)],[0 0],'-k')
% grid on;
% axis([0 length(BPS) -inf inf])
% xlabel('Time \rightarrow')
% ylabel('Lyapunov Exponent (bits/orbit)')
% figureHandle = gcf;
% set(findall(figureHandle,'type','text'),'fontSize',18,'fontWeight','bold')
LyaExp_b_sec = BPS(end);
LyaExp_b_orb = bitsPerOrbit(end);
end