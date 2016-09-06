
#define CALTAG_PS_CODE " \n\
 \n\
gsave \n\
 \n\
% hairlines for crop markers \n\
0.001 setlinewidth \n\
 \n\
% coordinate frame: \n\
% spacing is 1/8 inch (= 9 pt) * specified value \n\
9 spacing mul dup scale \n\
% move origin a bit to maximize space \n\
-2 -2 translate \n\
 \n\
% total number of markers \n\
/nummark { columns rows mul } bind def \n\
 \n\
% constant: 1/sqrt(2) \n\
/sqtwoinv { 1 2 sqrt div } bind def \n\
 \n\
% pixel at current position \n\
% no stack values consumed \n\
/pixel { \n\
    gsave \n\
    newpath \n\
    0 0 moveto 1 0 lineto 1 1 lineto 0 1 lineto \n\
    closepath fill \n\
    grestore \n\
} def \n\
 \n\
% row starting at current position, left to right \n\
% in stack: id \n\
% out stack: bit-shifted id \n\
/drawrow { \n\
    gsave \n\
    4 { \n\
	dup 1 and 1 eq {pixel} if \n\
	-1 0 translate \n\
	-1 bitshift \n\
    } repeat \n\
    grestore \n\
} def \n\
 \n\
% draw code for current marker \n\
% stack: id \n\
/drawcode { \n\
    5 2 translate \n\
    gsave \n\
    4 { \n\
	drawrow \n\
	0 1 translate \n\
    } repeat \n\
    pop \n\
    grestore \n\
} def \n\
 \n\
% draw marker \n\
% stack: codenumber x y \n\
/marker { \n\
    1 add exch 1 add exch \n\
     \n\
    layout 1 eq { \n\
	% checker-style unrotate pattern \n\
	gsave \n\
	2 copy \n\
	8 mul exch 8 mul exch translate \n\
	add 2 mod 1 eq { \n\
	    % black background for odd checkers \n\
	    0 setgray \n\
	    newpath \n\
	    0 0 moveto 8 0 lineto 8 8 lineto 0 8 lineto \n\
	    closepath fill \n\
	    % white boxes for 0 pixels \n\
	    1 setgray \n\
	    drawcode \n\
	} { \n\
	    % even checkers - black boxes for 1 pixels \n\
	    0 setgray \n\
	    not drawcode \n\
	} ifelse \n\
	grestore \n\
    } { \n\
	% rotated pattern \n\
	gsave \n\
	8 mul exch 8 mul 4 add exch translate \n\
	45 rotate \n\
	sqtwoinv dup scale \n\
	0 setgray \n\
	% black background octagon \n\
	newpath \n\
	1 0 moveto 7 0 lineto 8 1 lineto 8 7 lineto 7 8 lineto 1 8 lineto 0 7 lineto 0 1 lineto \n\
	closepath fill \n\
	% white squares for 0 pixels \n\
	1 setgray \n\
	drawcode \n\
	grestore \n\
    } ifelse \n\
} def \n\
 \n\
% draw a bowtie \n\
% stack: x y \n\
/bowtie { \n\
    gsave \n\
    8 mul exch 8 mul exch translate \n\
    0 setgray \n\
    newpath 0 0 moveto 0 sqtwoinv -4 mul dup 0 lineto lineto closepath fill \n\
    newpath 0 0 moveto 0 sqtwoinv 4 mul dup 0 lineto lineto closepath fill \n\
    grestore \n\
} def \n\
 \n\
% draw a horizontal boundary field for the checker layout \n\
% stack: x y \n\
/hBoundary { \n\
    gsave \n\
    8 mul exch 8 mul exch translate \n\
    0 setgray \n\
    newpath 0 0 moveto 8 0 lineto 8 4 lineto 0 4 lineto closepath fill \n\
    grestore \n\
} def \n\
 \n\
% draw a vertical boundary field for the checker layout \n\
% stack: x y \n\
/vBoundary { \n\
    gsave \n\
    8 mul exch 8 mul exch translate \n\
    0 setgray \n\
    newpath 0 0 moveto 0 8 lineto 4 8 lineto 4 0 lineto closepath fill \n\
    grestore \n\
} def \n\
 \n\
% draw a crop marker \n\
% stack: x y \n\
/cropmarker { \n\
    gsave \n\
    0 setgray \n\
    8 mul exch 8 mul exch translate \n\
    newpath -2 0 moveto 2 0 lineto stroke \n\
    newpath 0 -2 moveto 0 2  lineto stroke \n\
    newpath 0 0 0.2 0 360 arc closepath stroke \n\
    grestore \n\
} def \n\
 \n\
 \n\
% \n\
% actual drawing \n\
% \n\
 \n\
gsave \n\
 \n\
% draw markers (both layouts) \n\
0 1 rows 1 sub { \n\
    0 1 columns 1 sub { \n\
	2 copy 2 copy \n\
	exch columns mul add \n\
	markerids exch get 3 1 roll exch marker \n\
	pop \n\
    } for \n\
     \n\
    pop \n\
} for \n\
 \n\
% layout-specific elements \n\
layout 1 eq { \n\
    % draw all possible boundary elements in the checker \n\
     \n\
    % draw horizontal boundary fields \n\
    1 1 columns { \n\
	% draw odd-column boundary fields in row 0 \n\
	dup 2 mod 1 eq { \n\
	    dup 0.5 hBoundary \n\
	} if \n\
	% draw bundary fields with even (column+height) \n\
	dup rows add 2 mod 0 eq { \n\
	    dup rows 1 add hBoundary \n\
	} if \n\
    } for \n\
     \n\
    % draw vertical boundary fields \n\
    1 1 rows { \n\
	% draw odd-row boundary fields in row 0 \n\
	dup 2 mod 1 eq { \n\
	    dup	0.5 exch vBoundary \n\
	} if \n\
	% draw bundary fields with even (row+width) \n\
	dup columns add 2 mod 0 eq { \n\
	    dup	columns 1 add exch vBoundary \n\
	} if \n\
	1 add \n\
    } for \n\
     \n\
    % draw applicable corners \n\
    columns 2 mod 0 eq { \n\
	% bottom right corner \n\
	gsave \n\
	columns 1 add 8 mul 8 translate \n\
	newpath 0 0 moveto 4 0 lineto 4 -4 lineto 0 -4 lineto closepath fill \n\
	grestore \n\
    } if \n\
    rows 2 mod 0 eq { \n\
	% top left corner \n\
	gsave \n\
	8 rows 1 add 8 mul translate \n\
	newpath 0 0 moveto 0 4 lineto -4 4 lineto -4 0 lineto closepath fill \n\
	grestore \n\
    } if \n\
    rows columns add 2 mod 1 eq { \n\
	% top right corner \n\
	gsave \n\
	columns 1 add 8 mul rows 1 add 8 mul translate \n\
	newpath 0 0 moveto 0 4 lineto 4 4 lineto 4 0 lineto closepath fill \n\
	grestore \n\
    } if \n\
} { \n\
    % draw bowties for rotated grid \n\
    1 1 columns 1 add { \n\
	1 1 rows 1 add { \n\
	    2 copy bowtie \n\
	    pop \n\
	} for \n\
	pop \n\
    } for \n\
} ifelse \n\
 \n\
 \n\
% optional paper crop markers \n\
drawcropmarkers { \n\
    .5 rows 1.5 add cropmarker \n\
    columns 1.5 add rows 1.5 add cropmarker \n\
    .5 .5 cropmarker \n\
    columns 1.5 add .5 cropmarker \n\
 \n\
    % dot identifying lower left corner \n\
    newpath 5 5 0.2 0 360 arc closepath fill \n\
} if \n\
 \n\
grestore \n\
 \n\
grestore \n\
 \n\
showpage \n\
"
