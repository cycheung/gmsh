set style data line
set style function lines
set zeroaxis
set grid


set title "Temperature distribution at OL.get(PostPro/PROBETIME) ms"
set key top right
set xlabel "Radial coordinate [mm]"
set xlabel "Time [s]"

#nbfiles= 3
#filename(n) = sprintf("templaser%d.txt", n)
#plot for [i=0:nbfiles-1] filename(i) u ($5)*1000:8 w l t "", \
#     OL.get(PostPro/5OVERTEMP) w l  lt rgb "black" t ""

plot "templaser0.txt" u ($5)*1000:8 w l t "surf.", \
     "templaser1.txt" u ($5)*1000:8 w l t "interf.", \
     "templaser2.txt" u ($5)*1000:8 w l t "depth=OL.get(PostPro/4ZSURF) mm", \
     OL.get(PostPro/5OVERTEMP) w l  lt rgb "black" t ""

