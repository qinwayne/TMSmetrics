% Write a Binary Array File
% Extended to support N-dimensional arrays
%
% Usage: write_ar (X, name)

function write_ar (X, name)

ndim = length(size(X));
if (ndim==2)
	[rows cols] = size(X);
	fid = fopen (name, 'w', 'ieee-be');
	fwrite (fid, rows, 'int');
	fwrite (fid, cols, 'int');
	fwrite (fid, X', 'double')';
	fclose (fid);
else
	X = permute (X, ndim:-1:1);
	fid = fopen (name, 'w', 'ieee-be');
	fwrite (fid, [-1 ndim], 'int');
	SIZ = size(X);
	fwrite (fid, SIZ(ndim:-1:1), 'int');
	pp = prod(size(X));
	fwrite (fid, X(1:pp), 'double');
	fclose (fid);
end
