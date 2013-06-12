set style data line
set style function lines
set zeroaxis
set grid

set title "Temperature on axis"
set ylabel "Temperature [{\260}C]";
set xlabel "Time [s]"
set ytics nomirror
set y2tics
set y2label "Volume above thresh. [mm^3]"
set key bottom right
plot "tempOL.get(TAGSIMU).txt" u 1:($2) lt 1 t "surf.",\
     "tempOL.get(TAGSIMU).txt" u 1:($11) lt 2 t "interf.",\
     OL.get(PostPro/5OVERTEMP) w l  lt rgb "black" t "",\
     "volume.txt" u 2:8 axis x1y2 lt 8 lw 1 t "vol.A.T.",\
     "volderm.txt" u 2:8 axis x1y2 lt 8 lw 2 t ""

set ytics mirror
unset y2tics
unset y2label

