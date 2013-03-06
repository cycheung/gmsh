OL.block
OS.number(2,Configuration/); OS.setVisible(0);
OS.valueLabels(1,"Linux",2,"MAC-OS",3,"Windows");
OL.endblock


INTERACT=-1

OL.if(OL.get(OS) == 2)
if (INTERACT==-1) set terminal aqua enhanced
OL.else
if (INTERACT==-1) set terminal x11 persist
OL.endif

if (INTERACT== 0) set terminal pdf enhanced font "Times-Roman,6"
if (INTERACT== 0) set output "plot.pdf"

set style data line
set zeroaxis
set multiplot;          
 
set style function lines
set size 1.0, 1.0
set origin 0.0, 0.0
set grid

if (INTERACT==-1) set title "Temperature on axis at surface and interface"
set size 0.5,0.5
set origin 0.0,0.5
set ylabel "Temperature [{\260}C]";
set xlabel "Time [s]"
set ytics nomirror
set y2tics
set y2label "Activated volume [mm^3]"
plot "temp.txt" u 1:($2) lt 1 t "",\
     "temp.txt" u 1:($11) lt 2 t "",\
     OL.get(Parameters/Skin/OVERTEMP) w l  lt rgb "black" t "threshold ",\
     "volume.txt" u 2:8 axis x1y2 lt 8 lw 1 t "act. vol.",\
     "volderm.txt" u 2:8 axis x1y2 lt 8 lw 2 t ""

set ytics mirror
unset y2tics
unset y2label
 
if (INTERACT==-1) set title "Temperature distribution at t_{laser} at surface and interface"
set size 0.5,0.5
set origin 0.5,0.5  
set xlabel "Radial coordinate [mm]"

nbfiles= 2
filename(n) = sprintf("templaser%d.txt", n)
plot for [i=0:nbfiles-1] filename(i) u ($5)*1000:8 w l t "", \
     OL.get(Parameters/Skin/OVERTEMP) w l  lt rgb "black" t "threshold "

skinWidth = (OL.get(Parameters/Skin/2EPIDERMIS)+OL.get(Parameters/Skin/3DERMIS))/1000
dermis = OL.get(Parameters/Skin/2EPIDERMIS)

if (INTERACT==-1) set title "Activated area"
set size 0.5,0.5
set origin 0.0,0.0
set xlabel "Depth [mm]"
set ylabel "A_{A{/Symbol d}} [mm^2]"
set arrow from dermis,graph(0,0) to dermis,graph(1,1) nohead lt 8
plot "activeMax.txt" u (skinWidth-($6))*1000:($8)*10**6 w lp t ""

if (INTERACT==-1) set title "Duration of activation"
set size 0.5,0.5
set origin 0.5,0.0
set xlabel "Depth [mm]"
set ylabel "Duration {/Symbol D} t_{A{/Symbol d}} [s]"
plot [0:0.2] "duration.txt" u (skinWidth-($6))*1000:8 w lp t "" 

unset multiplot




