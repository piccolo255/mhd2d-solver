dir = "data/g64-t2/shlf-64-t2-beta0.2-u1.0.ini/"
outdir = "plots/detail-tokicon/"
FILES = "\
shlf-64-t2-beta0.2-u1.0-0160-05.0002441.dat \
"

g = 64
Lx = 32
Ly =  6

Nx = Lx * g
Ny = Ly * g
xmod = 1.0 / g / 2.0
ymod = 1.0 / g / 2.0

#set term pngcairo size 1200,600 enhanced font "Arial,14"
#ext = ".png"
set term postscript eps color enhanced size 6,3 lw 2 font "Arial,30"
ext = ".eps"

set title ""

do for [ file in FILES ] {
   set yrange [-1:1]
   
   ################################################################## pressure
   
   load "parula.pal"

   set cbrange [0:1.2]
   set cblabel "pressure"

   #set title "grid density = 64, {/Symbol t} = 2, {/Symbol b} = 0.2, u_{init} = 1.0; t = 1.00"

   set out outdir . file . "-left-p" . ext
   set xrange [-9:-5]
   plot dir.file binary record = (Nx,Ny) skip = 0 format = "%11float" \
        u ($1+xmod):($2+ymod):7 w image t ""
   
   set out outdir . file . "-right-p" . ext
   set xrange [0:5]
   replot
   
   ################################################################ y-velocity
   
   load "prgn.pal"

   set cbrange [-0.06:0.06]
   set cblabel "velocity v"

   #set title "grid density = 64, {/Symbol t} = 2, {/Symbol b} = 0.2, u_{init} = 1.0; t = 1.00"

   set out outdir . file . "-left-v" . ext
   set xrange [-9:-5]
   plot dir.file binary record = (Nx,Ny) skip = 0 format = "%11float" \
        u ($1+xmod):($2+ymod):5 w image t ""
   
   set out outdir . file . "-right-v" . ext
   set xrange [0:5]
   replot
}

set out
