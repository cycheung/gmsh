set zeroaxis
set style data linespoints

set multiplot;          

set size 1.0, 1.0
set origin 0.0, 0.0
set grid

set size 0.5,1
set origin 0,0
set xlabel "H/ref"
set ylabel "J/Jsat"
plot "results.txt" u 1:2 t""

set size 0.5,1
set origin 0.5,0
set xlabel "Pseudo time []"
set ylabel ""
plot "results.txt" u 1 t"H/Href",\
     "results.txt" u 2 t"J/Jsat"


