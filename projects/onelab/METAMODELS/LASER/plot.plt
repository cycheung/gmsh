
#set terminal pdf font "Times-Roman,12" ; INTERACT = 0

set terminal aqua; INTERACT=-1

set style data line
set zeroaxis
set ylabel "Temperature [K]"
set xlabel "time [s]"

plot "temp.txt" u 1:2 t"depth=0 mm",\
     "temp.txt" u 1:8 t"depth=0.05 mm",\
     "temp.txt" u 1:14 t"depth=0.10 mm",\
     "temp.txt" u 1:20 t"depth=0.15 mm",\
     "temp.txt" u 1:26 t"depth=0.020 mm"

set terminal aqua 1

set xlabel "coord [m]"

plot "temper.txt" u 5:8 w linesp



