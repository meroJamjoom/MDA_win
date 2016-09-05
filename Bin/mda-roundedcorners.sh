# ==========================================================================
# $Id: mda-roundedcorners.sh 79 2011-07-22 05:59:14Z heidrich $
# add an alpha channel to a 2D MDA file (image) to generate rounded corners
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

radius=20
mode="alpha"

############################ main code ###################

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
	'--alpha' )
	    mode="alpha";;
	'--transparency' )
	    mode="transparency";;
	* )
            echo 1>&2 "Usage: $0 <options>\n"
            echo 1>&2 "Add an alpha channel with rounded courners to an image."
            echo 1>&2 "Options:"
            echo 1>&2 "  -r | --radius\n\tradius of neighborhood ($radius)\n"
	    echo 1>&2 "  --alpha | --transparency\n\twhether to generate an alpha or a transparency (1-alpha) channel.\n\t($mode).\n"
	    exit 1;;
	esac
    fi
done

dir=${0%/*}

# temp files
tmpTemplate="/tmp/tmp.XXXXXXX"
rgbFile=`mktemp $tmpTemplate`
alphaFile=`mktemp $tmpTemplate`

# get image size from standard input, while storing MDA stream to temp file
dimensions=`tee $rgbFile | mda-info -d`

# make arithmetic expression for pre-multiplying alpha
channels=`mda-info -c < $rgbFile`
case $mode in
'alpha' )
    premult="#$channels";;
* )
    premult="1-#$channels";;
esac
for (( i=$channels-1 ; $i>=0 ; i-- )) ; do
    premult="#$i*#$channels,$premult"
done

# create rounded corner single channel array with same dimensions
${dir+$dir/}mda-newmda -v "$dimensions,$radius" -p "(#2>%2&&#2<=%0-%2)||(#3>%2&&#3<=%1-%2)||pow(#2-%2+1,2)+pow(#3-%2+1,2)<=%2*%2||pow(#2-%0+%2,2)+pow(#3-%2+1,2)<=%2*%2||pow(#2-%2+1,2)+pow(#3-%1+%2,2)<=%2*%2||pow(#2-%0+%2,2)+pow(#3-%1+%2,2)<=%2*%2?1:0" $dimensions | ${dir+$dir/}mda-filter -b cyclic -s .5 gauss > $alphaFile

# merge the original channels and the mask, write to standard out
${dir+$dir/}mda-mergechannels $rgbFile $alphaFile | mda-channelarith $premult

# clean up
rm $rgbFile $alphaFile
