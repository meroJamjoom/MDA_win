# ==========================================================================
# $Id: mda-frompnm.csh 79 2007-05-16 05:59:14Z heidrich $
# convert PNM image into MDA
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

set awk_command = 'BEGIN { w = -1 ; h = -1 ; c = -1 } /#/ {} /P5/ { c = 1 } /P6/ { c = 3 } /^[0-9]* [0-9]*/ { w = $1 ; h= $2 } /^255/ { print "-t ubyte -c", c, "-co", w "," h ; exit } /^65535/ { print "-t ushort -c", c, "-co", w ',' h ; exit } /.*/ {print "Cannot read image: not a binary PNM!" > "/dev/stderr" ; exit 1}'

if (! -r $1) then
    echo "Cannot read $1\n"
else
    $0:h/mda-fromraw -off -1 `awk "$awk_command" $1` $1
endif
