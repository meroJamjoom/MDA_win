# ==========================================================================
# $Id: mda-info.csh 106 2007-06-06 20:49:42Z atcheson $
# extract and print MDA header info
# ==========================================================================
# License: Internal use at UBC only! External use is a copyright violation!
# ==========================================================================
# (C)opyright:
#
# 2007-, UBC
#
# Creator: heidrich (Wolfgang Heidrich)
# Email:   heidrich@cs.ubc.ca
# ==========================================================================

set full_command = '/MDA/ {print "MDA version", $2} /Format:/ {print "Data type:\t",$2,"\nEncoding:\t",$3,"endian"} /Dimensions:/ {print "Dimensions:\t",$2} /Channels:/ {print "Channels:\t",$2} /###/ {exit}'

set version_command = '/MDA/ {print $2}'
set type_command = '/Format:/ {print $2}'
set order_command = '/Format:/ {print $3}'
set dim_command = '/Dimensions:/ {print $2}'
set channel_command = '/Channels:/ {print $2}'
set numdims_command = '/Dimensions:/ {numdim = split($2,dims,","); if( length(dims)==1 && dims[1]=="1") {print "0"} else {print length(dims)} }'

set command = ""

foreach f ($*)
    switch ($f)
    case "-v":
    case "--version":
	set command = "$command $version_command"
	breaksw;
    case "-t":
    case "--type":
	set command = "$command $type_command"
	breaksw;
    case "-b":
    case "--byte-order":
	set command = "$command $order_command"
	breaksw;
    case "-d":
    case "--dimension":
	set command = "$command $dim_command"
	breaksw;
    case "-c":
    case "--channels":
	set command = "$command $channel_command"
	breaksw;
    case "-n":
    case "--numdimensions":
	set command = "$command $numdims_command"
	breaksw
    default:
	set echo_style sysv
	echo "Usage: $0 <options>\n"
        echo "Print information on an MDA file. Without options, a summary of"
	echo "the MDA properties is written in human readable form. The"
	echo "options restrict the information to a single piece of"
	echo "information that is easy to parse in shell scripts."
	echo "\nOptions:\n"
	echo "  -v | --version\n\tprint file version\n"
	echo "  -t | --type\n\tprint data type\n"
	echo "  -b | --byte-order\n\tprint byte order\n"
	echo "  -d | --dimension\n\tprint array dimensions\n"
	echo "  -c | --channels\n\tprint number of channels\n"
	echo "  -n | --numdimensions\n\tprint number of dimensions\n"
	exit 1
    endsw
end

if ( "$command" == "" ) then
    set command = "$full_command"
else
    set command = "$command  /###/ {exit}"
endif
    
awk "$command"

# consume the rest of the input to avoid non-zero return values if the
# input stream is a pipe
cat > /dev/null
