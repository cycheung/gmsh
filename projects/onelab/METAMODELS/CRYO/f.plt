
set style data linespoints

set xlabel "Time (h)"
set ylabel "Concentration (kg/m3)"
#plot [0:168] "f.txt" u ($2)/3600:8  w l t ''
plot "f.txt" u ($2)/3600:8  w l t ''

