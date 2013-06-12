set style data line
set style function lines
set zeroaxis
set grid

skinWidth = OL.get(PostPro/4SKINWIDTH)*1e-3
dermis = OL.get(Parameters/Skin/2EPIDERMIS)

set title "Area above threshold at OL.get(PostPro/PROBETIME) ms"
set xlabel "Depth [mm]"
set ylabel "A_{A{/Symbol d}} [mm^2]"
set arrow from dermis,graph(0,0) to dermis,graph(1,1) nohead lt 8
plot "activeMax.txt" u (skinWidth-($6))*1000:($8*1e6) w lp t ""
unset arrow
