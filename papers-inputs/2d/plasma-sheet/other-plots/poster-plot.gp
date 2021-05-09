dir = "g64-t2/shlf-64-t2-beta0.2-u1.0.ini/"

FILES = "\
shlf-64-t2-beta0.2-u1.0-0032-01.0000000.dat \
shlf-64-t2-beta0.2-u1.0-0080-02.5000000.dat \
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

#set term pngcairo size 1200,600 enhanced font "Arial,14"
#ext = ".png"
set term postscript eps color enhanced size 6,3 lw 2 font "Arial,30"
ext = ".eps"

set title ""

do for [ file in FILES ] {
   set out file . "-d" . ext

   set cbrange [0:1.2]
   load "moreland.pal"

   #set title "grid density = 64, {/Symbol t} = 2, {/Symbol b} = 0.2, u_{init} = 1.0; t = 1.00"
   set cblabel "density"

   plot dir.file binary record = (Nx,Ny) skip = 0 format = "%11float" \
        u ($1+xmod):($2+ymod):3 w image t ""

   set out file . "-p" . ext
   load "greys.pal"

   set cbrange [0:1.2]
   set cblabel "pressure"

   #set title "grid density = 64, {/Symbol t} = 2, {/Symbol b} = 0.2, u_{init} = 1.0; t = 1.00"

   plot dir.file binary record = (Nx,Ny) skip = 0 format = "%11float" \
        u ($1+xmod):($2+ymod):7 w image t ""
}

set out
