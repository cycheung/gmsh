
INTERACT=-1

if (INTERACT==-1) set terminal aqua enhanced
if (INTERACT== 0) set terminal pdf enhanced font "Times-Roman,6"
if (INTERACT== 0) set output "plot.pdf"

set style data line
set zeroaxis
set multiplot;          
 
set style function lines
set size 1.0, 1.0
set origin 0.0, 0.0
set grid

if (INTERACT==-1) set title "Maximum temperature at different depths"
set size 0.5,0.5
set origin 0.0,0.5
set ylabel "Temperature [{\260}C]";
set xlabel "Time [s]"
plot "temp.txt" u 1:($2)-273 t "",\
     "temp.txt" u 1:($9)-273 t "",\
     "temp.txt" u 1:($16)-273 t "",\
     "temp.txt" u 1:($23)-273 t "",\
     "temp.txt" u 1:($30)-273 t "", \
     320-273 w l  lt rgb "black" t "threshold "

if (INTERACT==-1) set title "Temperature distribution at t_{laser} at different depths"
set size 0.5,0.5
set origin 0.5,0.5  
set xlabel "Radial coordinate [mm]"

nbfiles= OL.get(PostPro/ZSURF,choices.size())
filename(n) = sprintf("templaser%d.txt", n)
plot for [i=0:nbfiles-1] filename(i) u ($5)*1000:($8)-273 w l t "", \
     320-273 w l  lt rgb "black" t "threshold "

skinWidth = (OL.get(Parameters/Skin/EPIDERMIS)+OL.get(Parameters/Skin/DERMIS))/1000
zsurf=OL.get(PostPro/ZSURF);

if (INTERACT==-1) set title "Maximum (over time) activated surface"
set size 0.5,0.5
set origin 0.0,0.0
set xlabel "Depth [mm]"
set ylabel "Activated surface A_{A{/Symbol d}} [mm^2]"
plot "activeMax.txt" u (zsurf-($6))*1000:($8)*10**6 w lp t ""

if (INTERACT==-1) set title "Maximum (on axis) duration at threshold "
set size 0.5,0.5
set origin 0.5,0.0
set xlabel "Depth [mm]"
set ylabel "Duration {/Symbol D} t_{A{/Symbol d}} [s]"
plot [0:0.2] "duration.txt" u (skinWidth-($6))*1000:8 w lp t "" 

unset multiplot




