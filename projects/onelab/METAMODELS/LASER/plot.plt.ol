
set terminal pdf enhanced font "Times-Roman,6"
set output "plot.pdf"

set style data line
set zeroaxis
set multiplot;          
 
set style function lines
set size 1.0, 1.0
set origin 0.0, 0.0
set grid

set title "Temperature on axis"
set size 0.5,0.5
set origin 0.0,0.5
set ylabel "Temperature [{\260}C]";
set xlabel "Time [s]"
set ytics nomirror
set y2tics
set y2label "Volume above thresh. [mm^3]"
set key top left
plot "temp.txt" u 1:($2) lt 1 t "surf.",\
     "temp.txt" u 1:($11) lt 2 t "interf.",\
     OL.get(PostPro/5OVERTEMP) w l  lt rgb "black" t "",\
     "volume.txt" u 2:8 axis x1y2 lt 8 lw 1 t "vol.A.T.",\
     "volderm.txt" u 2:8 axis x1y2 lt 8 lw 2 t ""

set ytics mirror
unset y2tics
unset y2label

set title "Temperature distribution at OL.get(PostPro/PROBETIME) ms"
set size 0.5,0.5
set origin 0.5,0.5 
set key right
set xlabel "Radial coordinate [mm]"

#nbfiles= 3
#filename(n) = sprintf("templaser%d.txt", n)
#plot for [i=0:nbfiles-1] filename(i) u ($5)*1000:8 w l t "", \
#     OL.get(PostPro/5OVERTEMP) w l  lt rgb "black" t ""

plot "templaser0.txt" u ($5)*1000:8 w l t "surf.", \
     "templaser1.txt" u ($5)*1000:8 w l t "interf.", \
     "templaser2.txt" u ($5)*1000:8 w l t "depth=OL.get(PostPro/4ZSURF) mm", \
     OL.get(PostPro/5OVERTEMP) w l  lt rgb "black" t ""


skinWidth = OL.get(PostPro/4SKINWIDTH)*1e-3
dermis = OL.get(Parameters/Skin/2EPIDERMIS)

set title "Area above threshold at OL.get(PostPro/PROBETIME) ms"
set size 0.5,0.5
set origin 0.0,0.0
set xlabel "Depth [mm]"
set ylabel "A_{A{/Symbol d}} [mm^2]"
set arrow from dermis,graph(0,0) to dermis,graph(1,1) nohead lt 8
plot "activeMax.txt" u (skinWidth-($6))*1000:($8*1e6) w lp t ""

set title "Duration above threshold"
set size 0.5,0.5
set origin 0.5,0.0
set xlabel "Depth [mm]"
set ylabel "Duration {/Symbol D} t_{A{/Symbol d}} [s]"
plot [0:0.2] "duration.txt" u (skinWidth-($6))*1000:8 w lp t "" 

unset multiplot

set output "plotOL.get(Solution/0TAG).pdf"
replot




