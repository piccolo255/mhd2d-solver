file = "data/rk3-order-scaling.dat"

outdir = "plots/"

#set term pngcairo size 800,600 lw 4 enhanced font "Arial,24"
#ext = ".png"
set term postscript eps color enhanced size 5,4 lw 2 font "Arial,30"
ext = ".eps"

set logscale xy
set ytics format "%.1e"
set xtics 4,4,1024
unset mxtics

set xrange [0.75*4:1.5*1024]

set title ""
set xlabel "number of steps, N"
set ylabel "error"
#set lmargin screen 0.12
#set rmargin screen 0.7

set yrange [1e-11:1e-1]
set ytics 1e-10,100,1e-2

set out outdir . file . "-rk3" . ext
plot file u 1:3 w lp t "RK3" \
   , x**-3 t "N^{-3}"

set yrange [1e-3:1e0]
set ytics autofreq

set out outdir . file . "-euler" . ext
plot file u 1:2 w lp t "Euler" \
   , x**-1 t "N^{-1}"

set out
