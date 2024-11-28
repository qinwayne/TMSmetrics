% Read Binary Array File
% Extended to support N-dimensional arrays
%
% Usage: X = read_ar (name, end_type)
%	end_type = 'ieee-be' (default)

function X = read_ar (name, end_type)

if (nargin<2)
	end_type = 'ieee-be';
end

fid = fopen (name, 'r', end_type);	% Sun is ieee-be
if (fid~=-1)
	rows = fread (fid, 1, 'int');
	cols = fread (fid, 1, 'int');
	if (rows~=-1)
		X = fread (fid, [cols,rows], 'double')';
	else
		S = fread (fid, [1,cols], 'int');
		S = S(cols:-1:1);
		X = zeros(S);
		X(1:prod(S)) = fread (fid, [1,Inf], 'double');
		X = permute(X,cols:-1:1);
	end
	fclose (fid);
else
	X = [];
end

