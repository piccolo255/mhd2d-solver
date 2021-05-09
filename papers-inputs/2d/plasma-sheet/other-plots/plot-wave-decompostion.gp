dir = "data/g64-t2/shlf-64-t2-beta0.2-u1.0.ini/"
outdir = "plots/"
FILES = "\
waves-shlf-64-t2-beta0.2-u1.0-0128-04.0002441-vs_3.75.dat \
waves-shlf-64-t2-beta0.2-u1.0-0128-04.0002441-vs_3.94.dat \
"

g = 64
Lx = 32
Ly =  6

Nx = Lx * g
Ny = Ly * g
xmod = 1.0 / g / 2.0
ymod = 1.0 / g / 2.0

set xrange [-2:6]
set yrange [-2:2]

set term pngcairo size 1200,600 enhanced font "Arial,14"
ext = ".png"
#set term postscript eps color enhanced size 6,3 lw 2 font "Arial,30"
#ext = ".eps"

set title ""

do for [ file in FILES ] {
   set cbrange [0:*]
   load "greys.pal"
   
   set cblabel "component strength"
   
   # x
   set out outdir . file . "-x-1-neg-fast" . ext
   plot dir.file u ($1+xmod):($2+ymod):(abs($3)) w image notitle
   
   set out outdir . file . "-x-2-neg-alfven" . ext
   plot dir.file u ($1+xmod):($2+ymod):(abs($4)) w image notitle
   
   set out outdir . file . "-x-3-neg-slow" . ext
   plot dir.file u ($1+xmod):($2+ymod):(abs($5)) w image notitle
   
   set out outdir . file . "-x-4-entropy" . ext
   plot dir.file u ($1+xmod):($2+ymod):(abs($6)) w image notitle
   
   set out outdir . file . "-x-5-pos-slow" . ext
   plot dir.file u ($1+xmod):($2+ymod):(abs($7)) w image notitle
   
   set out outdir . file . "-x-6-pos-alfven" . ext
   plot dir.file u ($1+xmod):($2+ymod):(abs($8)) w image notitle
   
   set out outdir . file . "-x-7-pos-fast" . ext
   plot dir.file u ($1+xmod):($2+ymod):(abs($9)) w image notitle
   
   # y
   set out outdir . file . "-y-1-neg-fast" . ext
   plot dir.file u ($1+xmod):($2+ymod):(abs($11)) w image notitle
   
   set out outdir . file . "-y-2-neg-alfven" . ext
   plot dir.file u ($1+xmod):($2+ymod):(abs($12)) w image notitle
   
   set out outdir . file . "-y-3-neg-slow" . ext
   plot dir.file u ($1+xmod):($2+ymod):(abs($13)) w image notitle
   
   set out outdir . file . "-y-4-entropy" . ext
   plot dir.file u ($1+xmod):($2+ymod):(abs($14)) w image notitle
   
   set out outdir . file . "-y-5-pos-slow" . ext
   plot dir.file u ($1+xmod):($2+ymod):(abs($15)) w image notitle
   
   set out outdir . file . "-y-6-pos-alfven" . ext
   plot dir.file u ($1+xmod):($2+ymod):(abs($16)) w image notitle
   
   set out outdir . file . "-y-7-pos-fast" . ext
   plot dir.file u ($1+xmod):($2+ymod):(abs($17)) w image notitle
}

set out
