set xlabel 'Time (s)'
set ylabel 'Omega'
#set logscale y

plot 'omega.txt' u 1:3 w lp t 'surface', 'omega.txt' u 1:5 w lp t 'interface'