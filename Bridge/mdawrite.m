function mdawrite( X, filename, multichannel )
%MDAWRITE Write MATLAB array to MDA file
%
% syntax: MDAWRITE( X, filename, multichannel )
%
%	X: the array data (matrix)
%	filename: destination filename (including .mda suffix)
%	multichannel: bool. If true the final dimension of X is treated as the
%		number of channels. An MxNxZ array could either be an MxNxZ scalar
%		array, or an MxN array of Z channels.

% inspect the array to see what data type we're dealing with
S = whos( 'X' );
bytesPerElement = S.bytes / numel( X );
hasNegatives = min( X(:) ) < 0;

% select the appropriate MDA types based on the array data
switch bytesPerElement
	case 1
		if hasNegatives
			mdaFormat = 'byte';
			precision = 'schar';
		else
			mdaFormat = 'ubyte';
			precision = 'uchar';
		end
	case 2
		if hasNegatives
			mdaFormat = 'short';
			precision = 'short';
		else
			mdaFormat = 'ushort';
			precision = 'ushort';
		end
	case 4
		if isfloat( X )
			mdaFormat = 'float';
			precision = 'float32';
		else
			if hasNegatives
				mdaFormat = 'int';
				precision = 'int';
			else
				mdaFormat = 'uint';
				precision = 'uint';
			end
		end
	case 8
		mdaFormat = 'double';
		precision = 'float64';
end

% split dimensions and channels
dim = size( X );
if multichannel
	numChannels = dim(end);
else
	numChannels = 1;
end

% transpose the first two dimensions because matlab indexes arrays by
% (row,column) and we tend to index MDAs by (x,y) ie in the horizontal
% dimension first
%if (multichannel && length( dim ) > 2) || (~multichannel && length( dim ) > 1)
%	seq = 1:length( dim );
%	seq(1:2) = [2 1];
%	X = permute( X, seq );
%	dim(1:2) = dim(2:-1:1);
%end

% interleave channels
if multichannel
	dim = dim(1:end-1);
	numPixels = prod( dim );
	X = reshape( X, [numPixels numChannels] )';
end

% convert "[M N P]" into "M,N,P"
dimStr = mat2str( dim );
dimStr = strrep( dimStr(2:end-1), ' ', ',' );

% write the header
f = fopen( filename, 'w' );
fprintf( f, 'MDA 1.0\n' );
fprintf( f, 'Format:\t\t%s little\n', mdaFormat );
fprintf( f, 'Dimensions:\t%s\n', dimStr );
fprintf( f, 'Channels:\t%d\n', numChannels );
fprintf( f, '###\n' );
fclose( f );

% write the data (in little endian mode)
f = fopen( filename, 'a', 'ieee-le' );
fwrite( f, X, precision, 0, 'ieee-le' );
fclose( f );

