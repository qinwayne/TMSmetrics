% Determine min of 2D array and its position
%
% Usage: [M row col] = minmin(X)
%
function [M, row, col] = minmin(X)

[a b] = min(X);
[M col] = min(a);
row = b(col);

