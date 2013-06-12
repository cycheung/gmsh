set terminal pdf enhanced font "Times-Roman,6"

set output "plotOL.get(TAGSIMU).pdf"

set multiplot     
set size 1.0, 1.0
set origin 0.0, 0.0

set size 0.5,0.5
set origin 0.0,0.5
load "plot1.plt"

set size 0.5,0.5
set origin 0.5,0.5 
load "plot2.plt"

set size 0.5,0.5
set origin 0.0,0.0
load "plot3.plt"

set size 0.5,0.5
set origin 0.5,0.0
load "plot4.plt"
unset multiplot


set output "plot.pdf"

set multiplot     
set size 1.0, 1.0
set origin 0.0, 0.0

set size 0.5,0.5
set origin 0.0,0.5
load "plot1.plt"

set size 0.5,0.5
set origin 0.5,0.5 
load "plot2.plt"

set size 0.5,0.5
set origin 0.0,0.0
load "plot3.plt"

set size 0.5,0.5
set origin 0.5,0.0
load "plot4.plt"
unset multiplot
