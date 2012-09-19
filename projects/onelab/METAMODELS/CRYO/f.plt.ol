
set title "fmin=OL.get(Solution/fmin) | Optimum application time=OL.get(Solution/tmin)"
set style data linespoints

set xlabel "Time (s)"
set ylabel "Objective function J(t)"
plot "f.txt" u 2:8  w l t ''

