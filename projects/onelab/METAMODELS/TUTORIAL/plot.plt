set zeroaxis
set style data linespoints

set xlabel "H [A/m]"
set ylabel "J [Tesla]"
plot "results.txt" u 1:2 t""

