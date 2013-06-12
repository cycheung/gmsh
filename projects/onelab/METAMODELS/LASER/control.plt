


set style data line
set style function lines
set key bottom

set title "Control"
set ylabel "Temperature [{\260}C]";
set xlabel "Time [s]"
set ytics nomirror
set y2tics
set y2label "Pin [W]"
plot "temp.txt" u 1:2 lt 1 lw 2 w linesp t "sensor",\
     "temp.txt" u 1:7 lt 2 t "imposed",\
     "temp.txt" u 1:(($6>0)?$6:0) axis x1y2 lt 0 t "Power"
set ytics mirror
unset y2tics
unset y2label
 
