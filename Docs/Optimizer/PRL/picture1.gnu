#!/usr/local/lib/gnuplot-3.8j/bin/gnuplot -persist
#
#    
#    	G N U P L O T
#    	Version 3.8j patchlevel 0
#    	last modified Wed Nov 27 20:49:08 GMT 2002
#    	System: Linux 2.4.20-30.9smp
#    
#    	Copyright(C) 1986 - 1993, 1999 - 2002
#    	Thomas Williams, Colin Kelley and many others
#    
#    	This is a pre-version of gnuplot 4.0. Please refer to the documentation
#    	for command syntax changes. The old syntax will be accepted throughout
#    	the 4.0 series, but all save files use the new syntax.
#    
#    	Type `help` to access the on-line reference manual
#    	The gnuplot FAQ is available from
#    		http://www.gnuplot.info/faq/
#    
#    	Send comments and requests for help to <info-gnuplot-beta@dartmouth.edu>
#    	Send bugs, suggestions and mods to <info-gnuplot-beta@dartmouth.edu>
#    
 set terminal postscript eps enhanced color blacktext \
   dashed dashlength 1.0 linewidth 2.2 defaultplex \
   palfuncparam 2000,0.003 \
   butt "Times Roman" 26
set output '1_7_6_Data.eps'
unset clip points
set clip one
unset clip two
set bar 1.000000
set border 31 lt -1 lw 1.000
set xdata
set ydata
set zdata
set x2data
set y2data
set boxwidth
set style fill empty border
set dummy x,y
set format x "% g"
set format y "% g"
set format x2 "% g"
set format y2 "% g"
set format z "% g"
set angles radians
unset grid
set key title ""
unset key
unset label
set label 1 "#1" at 119.162, 0.0127786, 0 centre norotate back nopoint
set label 2 "#2" at 117.103, 0.0249423, 0 centre norotate back nopoint
set label 3 "#3" at 115.045, 0.0268812, 0 centre norotate back nopoint
set label 4 "#4" at 110.927, 0.0230842, 0 centre norotate back nopoint
set label 5 "#5" at 106.715, 0.0167649, 0 centre norotate back nopoint
set label 6 "#6" at 103.403, 0.0107462, 0 centre norotate back nopoint
set label 7 "#7" at 101.853, 0.00705893, 0 centre norotate back nopoint
unset arrow
unset style line
unset style arrow
unset logscale
set offsets 0, 0, 0, 0
set pointsize 1
set encoding default
unset polar
unset parametric
unset decimalsign
set view 60, 30, 1, 1
set samples 100, 100
set isosamples 10, 10
set surface
unset contour
set clabel '%8.3g'
set mapping cartesian
set datafile separator whitespace
unset hidden3d
set cntrparam order 4
set cntrparam linear
set cntrparam levels auto 5
set cntrparam points 5
set size ratio 0 1,1
set origin 0,0
set style data points
set style function lines
set xzeroaxis lt -2 lw 1.000
set yzeroaxis lt -2 lw 1.000
set x2zeroaxis lt -2 lw 1.000
set y2zeroaxis lt -2 lw 1.000
set tics in
set ticslevel 0.5
set ticscale 1 0.5
set mxtics default
set mytics default
set mztics default
set mx2tics default
set my2tics default
set mcbtics default
#set xtics border mirror norotate 5
set xtics border mirror norotate 100,5
set ytics border mirror norotate autofreq 
set ztics border nomirror norotate autofreq 
set nox2tics
set noy2tics
set cbtics border mirror norotate autofreq 
set title "" 0.000000,0.000000  font ""
set timestamp "" bottom norotate 0.000000,0.000000  ""
set rrange [ * : * ] noreverse nowriteback  # (currently [0.00000:10.0000] )
set trange [ * : * ] noreverse nowriteback  # (currently [-5.00000:5.00000] )
set urange [ * : * ] noreverse nowriteback  # (currently [-5.00000:5.00000] )
set vrange [ * : * ] noreverse nowriteback  # (currently [-5.00000:5.00000] )
set timefmt x "%d/%m/%y,%H:%M"
set timefmt y "%d/%m/%y,%H:%M"
set timefmt z "%d/%m/%y,%H:%M"
set timefmt x2 "%d/%m/%y,%H:%M"
set timefmt y2 "%d/%m/%y,%H:%M"
set timefmt cb "%d/%m/%y,%H:%M"
set xlabel "angle (degree)" 0.000000,0.000000  font ""
set x2label "" 0.000000,0.000000  font ""
#set xrange [ * : * ] noreverse nowriteback  # (currently [-10.0000:10.0000] )
set xrange [98:120] noreverse nowriteback  # (currently [-10.0000:10.0000] )
set x2range [ * : * ] noreverse nowriteback  # (currently [-10.0000:10.0000] )
set ylabel "gradient (a.u.)" 0.000000,0.000000  font ""
set y2label "" 0.000000,0.000000  font ""
set yrange [ 0.00000 : 0.0300000 ] noreverse nowriteback
set y2range [ * : * ] noreverse nowriteback  # (currently [-10.0000:10.0000] )
set zlabel "" 0.000000,0.000000  font ""
set zrange [ * : * ] noreverse nowriteback  # (currently [-10.0000:10.0000] )
set cblabel "" 0.000000,0.000000  font ""
set cbrange [ * : * ] noreverse nowriteback  # (currently [-10.0000:10.0000] )
set zero 1e-08
set lmargin -1
set bmargin -1
set rmargin -1
set tmargin -1
set locale "C"
set pm3d scansautomatic flush begin noftriangles nohidden3d transparent implicit corners2color mean
unset pm3d
set palette positive nops_allcF maxcolors 0 gamma 1.5 color model RGB 
set palette rgbformulae 7, 5, 15
set colorbox default
set colorbox vertical origin 0.9,0.2 size 0.1,0.63 bdefault
set loadpath 
set fontpath 
plot "QUADR/1_7_6_Data" u 2:3 w points pointtype 4 pointsize 3,-.1781331503+.0017985874*x, "QUADR/1_7_6_Pred2" u 1:2 w points pointtype 3 pointsize 3
#    EOF
