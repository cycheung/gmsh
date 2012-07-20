#set terminal pdf font "Times-Roman,12" ; INTERACT = 0

set terminal aqua; INTERACT=-1
set terminal aqua 1

set style data line
set zeroaxis
set multiplot;          
 
set style function lines
set size 1.0, 1.0
set origin 0.0, 0.0

set multiplot
set grid

set title "Maximum skin temperature at different depths"
set size 0.5,0.5
set origin 0.0,0.5
set ylabel "Temperature [degC]"
set xlabel "Time [s]"
plot "temp.txt" u 1:($2)-273 t "",\
     "temp.txt" u 1:($9)-273 t "",\
     "temp.txt" u 1:($16)-273 t "",\
     "temp.txt" u 1:($23)-273 t "",\
     "temp.txt" u 1:($30)-273 t "", \
     320-273 t "threshold "

set title "Skin temperature (at t=Tlaser) at different depth "
set size 0.5,0.5
set origin 0.5,0.5  
set ylabel "Temperature [degC]"
set xlabel "Radial coord [mm]"
plot "templaser0.txt" u ($5)*1000:($8)-273 w l t "",\
     "templaser1.txt" u ($5)*1000:($8)-273 w l t "",\
     "templaser2.txt" u ($5)*1000:($8)-273 w l t "",\
     "templaser3.txt" u ($5)*1000:($8)-273 w l t "",\
     "templaser4.txt" u ($5)*1000:($8)-273 w l t "", \
     320-273 t "threshold "

skinWidth = (0.05+1.5)/1000

zsurf=0.00155;

set title "Maximum (in time) active surface"
set size 0.5,0.5
set origin 0.0,0.0
set xlabel "Skin Depth [mm]"
set ylabel "Active surface [mm^2]"
plot "activeMax.txt" u (zsurf-($6))*1000:($8)*10**6 w lp t ""

set title "Maximum (at x=0) duration at threshold "
set size 0.5,0.5
set origin 0.5,0.0
set xlabel "Skin Depth [mm]"
set ylabel "Duration [s]"
plot [0:0.2] "duration.txt" u (skinWidth-($6))*1000:8 w lp t "" 

unset multiplot





