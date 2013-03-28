
#set logscale y

set multiplot;          

set size 0.5,0.5
set origin 0.0,0.0
set xlabel 'Time (s)'
set ylabel 'Omega'
plot 'omega.txt' u 1:3 w lp t 'surface', 'omega.txt' u 1:5 w lp t 'interface'

set size 0.5,0.5
set origin 0.5,0.0
set xlabel 'Time (s)'
set ylabel 'CEM'
plot 'CEM.txt' u 1:3 w lp t 'surface', 'CEM.txt' u 1:5 w lp t 'interface'