set style data line
set style function lines
set zeroaxis
set grid

set title "Duration above threshold"
set xlabel "Depth [mm]"
set ylabel "Duration {/Symbol D} t_{A{/Symbol d}} [s]"
set arrow from dermis,graph(0,0) to dermis,graph(1,1) nohead lt 8
plot [0:0.2] "duration.txt" u (skinWidth-($6))*1000:8 w lp t "" 
unset arrow