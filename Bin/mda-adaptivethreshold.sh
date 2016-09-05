# ==========================================================================
# $Id: mda-adaptivethreshold.sh 79 2011-07-22 05:59:14Z heidrich $
# adpative thresholding of an MDA array (single channel only)
# ==========================================================================
# License: Internal use at UBC only! External use is a copyright violation!
# ==========================================================================
# (C)opyright:
#
# 2011-, UBC
#
# Creator: heidrich (Wolfgang Heidrich)
# Email:   heidrich@cs.ubc.ca
# ==========================================================================


############################ default variable values #####

channel=0
radius=20
percentileThreshold=.5
noiseLevel=.1
background=0

############################ main code ###################

dir=${0%/*}

var=""
for f in "$@" ; do
    if [ "$var" != "" ]
    then
	eval "$var=\"$f\""
	var=""
    else
	case $f in
	'-r' | '--radius' )
	    var="radius";;
	'-c' | '--channel' )
	    var="channel";;
	'-p' | '--percentile-threshold')
	    var="percentileThreshold";;
	'-n' | '--noise-level')
	    var="noiseLevel";;
	'-b' | '--background')
	    var="background";;
	*)
	    echo 1>&2 "Usage: $0 <options>\n"
	    echo 1>&2 "Adaptive thresholding of a single channel."
	    echo 1>&2 "Options:"
	    echo 1>&2 "  -r | --radius\n\tradius of neighborhood ($radius)\n"
	    echo 1>&2 "  -c | --channel\n\tinput channel ($channel)\n"
	    echo 1>&2 "  -p | --percentile-threshold\n\tthreshold for each neighborhood ($percentileThreshold)\n"
	    echo 1>&2 "  -n | --noise-level\n\testimate of the noise level in the image ($noiseLevel)\n"
	    echo 1>&2 "  -b | --background\n\tbackground color ($background)\n"
	    echo 1>&2 "  -h | --help\n\tprint this message\n"
	    exit 1;;
	esac
    fi
done

# Method:
# - generate an array with three copies of the channel to be thresholded
# - erode channel 0 by the radus (i.e min over a local region)
# - dilate channel 1 by the radius (i.e. max over the region)
# - for each pixel
#     if min-max is < than noise level
#       flat image region - set to background
#     else
#       if pixel value is above local threshold (percetile between min and max)
#         set to 1
#       else
#         set to 0
${dir+$dir/}mda-swizzle -cl $channel,$channel,$channel | ${dir+$dir/}mda-filter -b mirror -r $radius -cl 0 erode -cl 2 dilate | ${dir+$dir/}mda-channelarith "#2-#0 < $noiseLevel ? $background :(#1< $percentileThreshold * #2+( 1-$percentileThreshold )*#0?0:1)"
