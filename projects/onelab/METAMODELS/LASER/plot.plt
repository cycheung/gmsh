
#set terminal pdf font "Times-Roman,12" ; INTERACT = 0

set terminal aqua; INTERACT=-1

set style data line
set zeroaxis
set ylabel "Temperature [K]"
set xlabel "Time [s]"

plot "temp.txt" u 1:2 t"depth=0 mm",\
     "temp.txt" u 1:8 t"depth=0.05 mm",\
     "temp.txt" u 1:14 t"depth=0.10 mm",\
     "temp.txt" u 1:20 t"depth=0.15 mm",\
     "temp.txt" u 1:26 t"depth=0.20 mm", \
     320 t "seuil Ad"

set terminal aqua 1

set xlabel "coord [mm]"
plot "tempsurf.txt" u ($5)*1000:8 w linesp t "surface temp"


set terminal aqua 1

set xlabel "depth [mm]"
plot "activeMax0.txt" u ($6)*0+0:8 w p  t "",  "activeMax1.txt" u ($6)*0+0.05:8 w p  t "",  "activeMax2.txt" u ($6)*0+0.1:8 w p t "", "activeMax3.txt" u ($6)*0+0.15:8 w p t "", "activeMax4.txt" u ($6)*0+0.2:8 w p  t ""





