# ==========================================================================
# $Id: mda-zoneplate.csh 79 2007-05-16 05:59:14Z heidrich $
# zone plate pattern
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


# radius of outer ring in normalized coordinates
set radius = .9

# scale of pattern
set scale = 17.5

# resolution
set res = 256

#mda-newmda -t float -v "sqrt(#0*#0+#1*#1)/$radius" -p "%0<1?((1+cos(2*pi*$scale*%0*%0))/2):0"  $res,$res > zoneplate2d.mda
mda-newmda -t float -v "sqrt(#0*#0+#1*#1+#2*#2)/$radius" -p "%0<1?((1+cos(2*pi*$scale*%0*%0))/2):0"  $res,$res,$res > zoneplate3d.mda
