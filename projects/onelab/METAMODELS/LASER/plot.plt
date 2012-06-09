
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
set size 1.0,0.5
set origin 0.0,0.5
set grid

set title "Skin Temperature of a laser CO2 for different depths"
set ylabel "Temperature [K]"
set xlabel "Time [s]"
plot "temp.txt" u 1:2 t"z0",\
     "temp.txt" u 1:8 t"z1",\
     "temp.txt" u 1:14 t"z2",\
     "temp.txt" u 1:20 t"z3",\
     "temp.txt" u 1:26 t"z4", \
     320 t "treshold Ad"


set title "Surface Temperature at Time=Tlaser"
set size 0.5,0.5
set origin 0.0,0.0  
set xlabel "coord [mm]"
plot "tempsurf.txt" u ($5)*1000:8 w linesp t "z0"

set title "Active surface "
set size 0.5,0.5
set origin 0.5,0.0
set ylabel "Active surface [mm^2]"
plot "activeMax.txt" u ($4)*0.05:($8)*10**6 w lp t "Fiber Ad"

unset multiplot


