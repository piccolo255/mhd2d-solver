dir1 = "data/g16-t2/shlf-16-t2-beta0.2-u1.0.ini/"
dir2 = "data/g32-t2/shlf-32-t2-beta0.2-u1.0.ini/"
dir3 = "data/g64-t2/shlf-64-t2-beta0.2-u1.0.ini/"
file1 = "shlf-16-t2-beta0.2-u1.0-0080-02.5000000.dat"
file2 = "shlf-32-t2-beta0.2-u1.0-0080-02.5000000.dat"
file3 = "shlf-64-t2-beta0.2-u1.0-0080-02.5000000.dat"

g1 = 16
g2 = 32
g3 = 64
Lx = 32
Ly =  6
wherex = 2.0
x0 = -16.0

Nx1 = Lx * g1
Ny1 = Ly * g1
Nx2 = Lx * g2
Ny2 = Ly * g2
Nx3 = Lx * g3
Ny3 = Ly * g3
xmod1 = 1.0 / g1 / 2.0
ymod1 = 1.0 / g1 / 2.0
xmod2 = 1.0 / g2 / 2.0
ymod2 = 1.0 / g2 / 2.0
xmod3 = 1.0 / g3 / 2.0
ymod3 = 1.0 / g3 / 2.0

block_size = 11*4
where1 = (wherex - x0) * g1 * Ny1
where2 = (wherex - x0) * g2 * Ny2
where3 = (wherex - x0) * g3 * Ny3
offset1 = where1 * block_size
offset2 = where2 * block_size
offset3 = where3 * block_size

set xrange [-1:1]
set yrange [0:1.4]

#set term pngcairo size 1200,600 enhanced font "Arial,14"
#ext = ".png"
set term postscript eps color enhanced size 3,2 lw 2 font "Arial,15"
ext = ".eps"

set title ""
set xlabel "y_{}"
set ylabel "pressure"
set title offset 0,-0.5
set xlabel offset 0,0.0
set ylabel offset 2
set grid

#load "jet.pal"

set out "compare-16+32+64-" . file1 . "-p" . ext

plot dir1.file1 binary record = (Ny1) skip = offset1 format = "%11float" u ($2+ymod1):($7) w l t "grid density 16" lc "blue" lw 2 dt 6 \
   , dir2.file2 binary record = (Ny2) skip = offset2 format = "%11float" u ($2+ymod2):($7) w l t "grid density 32" lc "red" lw 2 dt 2 \
   , dir3.file3 binary record = (Ny3) skip = offset3 format = "%11float" u ($2+ymod3):($7) w l t "grid density 64" lc "black" lw 2 dt 1
set out
