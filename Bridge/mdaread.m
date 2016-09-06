function [mda,succeed] = mdaread( filename )
%READMDA Load an MDA file as a Matlab matrix
%   This will convert an MxNxPx... MDA with C channels to a native Matlab
%   array of size MxNxPx...xC
%   Note that Matlab indexes rows as the first dimension, whereas we typically
%   use the first dim for columns in images. The caller may therefore need to
%   transpose the first two dimensions in the returned array.
%   Returns succeed=0 if ok, nonzero on error

	% dummy initial value
	mda = [];
	succeed = -1;
    
    % check if file exists
    if(~exist(filename, 'file'))
        error(sprintf('File does not exist: %s', filename));
        return;
    end

    % unzip .mda.gz files
    % check that filename is >7 chars, ends in .mda.gz
    % check if .mda version already exists, if so, use it
    isZip=false;
    if( length(filename)>7 && strcmpi(filename(end-6:end), '.mda.gz') )
        if( ~exist(filename(1:end-3), 'file') )  % .mda does not exist
            isZip=true;
            system(sprintf('gunzip %s', filename));
        end
        filename=filename(1:end-3);
    end

	% check that it's an MDA file
	f = fopen( filename, 'r' );
	signature = fgetl( f );
	if( ~strcmp(upper(signature(1:3)),'MDA') )
		fclose( f );
        if(isZip)
            system(sprintf('gzip %s', filename));
        end
		return;
	end
	
	% read header
	format = lower( fgetl( f ) );
	if( strfind(format,'ubyte') )
		datatype = 'uchar';
		datasize = 1;
	elseif( strfind(format,'byte') )
		datatype = 'schar';
		datasize = 1;
	end
	if( strfind(format,'ushort') )
		datatype = 'uint16';
		datasize = 2;
	elseif( strfind(format,'short') )
		datatype = 'int16';
		datasize = 2;
	end
	if( strfind(format,'uint') )
		datatype = 'uint32';
		datasize = 4;
	elseif( strfind(format,'int') )
		datatype = 'int32';
		datasize = 4;
	end
	if( strfind(format,'float') )
		datatype = 'float32';
		datasize = 4;
	end
	if( strfind(format,'double') )
		datatype = 'float64';
		datasize = 8;
	end
	if( strfind(format,'big') )
		endian = 'ieee-be';
	else
		endian = 'ieee-le';
	end
	dim = fgetl( f );
	dim = dim( 12:end );
	dim = sscanf( dim, '%d,' )';
	channels = sscanf( fgetl(f), 'Channels: %d' );
	dummy = fgetl( f );

	% seek to start of data chunk. it's necessary to seek from the end, because continuing to read
	% after parsing the header produces bugs when certain characters show up in the stream (they can
	% be mistakenly treated as EOF or EOL characters)
	nbytes = prod( dim ) * channels * datasize;
	status = fseek( f, -nbytes, 'eof' );
	if status ~= 0
		disp( ferror( f ) );
		succeed = 1;
        if(isZip)
            system(sprintf('gzip %s', filename));
        end
		return;
	end

	% NB: don't try to preallocate zeros for 'mda' because then they will be doubles and consume
	% lots of memory even if the input type is uint. Also, sometimes the fread returns an incorrect
	% number of elements, and then the array sizes won't match.
	% NB: do not attempt to read in channels one at a time from the disk because a scattered load
	% such as "mda(1:channels:end) = fread()" is very slow and uses lot of memory
	mda = fread( f, inf, strcat('*',datatype), 0, endian );

	% UPDATE: this was fixed by seeking to n bytes from the end of file marker, rather than trying to
	% read the next n bytes after the header. /ENDOFUPDATE
	% handle bad input (when headers are incorrect). this is a bug in mda/matlab/linux? which causes
	% fread to return a different number of bytes from the file than actually exist (or at least, which
	% i think do exist. according to linux the file contains the correct amount of data, but we can't
	% always guarantee it all gets read in here. failure rate is about 1 in 500 in my experience)
	%n = prod( dim ) * channels;
	%if numel(mda) ~= n
	%	disp( 'Warning: the number of elements read from the file does not' );
	%	disp( 'match number indicated by header.' );
	%	dif = abs( n - numel(mda) );
	%	if numel(mda) < n
	%		disp( sprintf( '%d zeros have been appended', dif ) );
	%		mda = [mda; zeros(dif,1)];
	%	else
	%		disp( sprintf( '%d elements have been truncated', dif ) );
	%		mda = mda(1:end-dif);
	%	end
	%	succeed = 0;
	%end

	% reshape the array into correct array dimensions....
	% i've also included some "bad" ways to do this, just to make sure the someone doesn't
	% make these mistakes again when modifying this code

%	% this is fast, uses no extra mem, but gives incorrect output on multichannel arrays
%	%mda = reshape( mda, [dim,channels] );

%	% this works generally, but uses extra mem since shiftdim is not in-place
%	%mda = shiftdim( reshape(mda,[channels,dim]), 1 ); 

	% this does the trick. it's not quite perfect though. uses about 35% of the total
	% array size in temp memory space. which is a big problem if you're trying to load
	% big video sequences with only 2Gb of RAM.
	mda = reshape( mda, [channels,prod(dim)] );
	mda = reshape( mda', [dim,channels] );

	% don't bother with transpose now, let caller decide whether to do that
	% ...



%	% this was the older version. it's very memory inefficient (lots of temp variables).
%	% but it does demonstrate a way to transpose the first two dimensions which could
%	% useful later on, so i'll leave the code here
%
%	% read input file as single column
%	data = fread( f, inf, strcat('*',datatype), 0, endian );
%
%	% handle bad input (when headers are incorrect)
%	n = prod( dim ) * channels;
%	if numel(data) ~= n
%		disp( 'Warning: the number of elements read from the file does not' );
%		disp( 'match number indicated by header.' );
%		dif = abs( n - numel(data) );
%		if numel(data) < n
%			disp( sprintf( '%d zeros have been appended', dif ) );
%			data = [data; zeros(dif,1)];
%		else
%			disp( sprintf( '%d elements have been truncated', dif ) );
%			data = data(1:end-dif);
%		end
%	end
%
%	% extract channels
%	for i = 1:channels
%		% reshape to match MDA dimensions
%		C = reshape( data(i:channels:end), dim );
%		% each 2D 'frame' must be transposed, because in matlab we index
%		% by (row,column) which is (y,x) and not (x,y). However we don't
%		% want to transpose the entire matrix (and the transpose of an
%		% nD matrix is undefined anyway) so just permute first 2 dimensions
%		if( length(dim) > 1 )
%			seq = 1:length(dim);
%			seq(1:2) = [2 1];
%			C = permute( C, seq );
%		end
%		% stack up the channels along the highest dimension
%		mda = squeeze( cat( length(dim)+1, mda, C ) );
%	end

	fclose( f );

	%if succeed == -1
	%	succeed = 1;
	%end
	n = prod( dim ) * channels;
	succeed = ( numel(mda) ~= n );

    % zip up .mda.gz file
    if(isZip)
        system(sprintf('gzip %s', filename));
    end

end
