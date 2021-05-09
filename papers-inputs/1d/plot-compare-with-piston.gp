FILES = "\
data/sh1d-128-u0.5.dat \
data/sh1d-128-u1.0.dat \
data/sh1d-128-u1.5.dat \
data/sh1d-128-u2.0.dat \
"

step = 128

#set terminal pngcairo size 600,400 lw 2 enhanced
#ext = ".png"
set term postscript eps color enhanced size 6,4 lw 4 font "Arial,30"
ext = ".eps"

outdir = "plots/"

dt = 0.03125

p = 1.0
r = 1.0
g = 5.0/3.0
cs = sqrt(g*p/r)

x1(up,t) = (cs-up*(g+1.0)/2.0)*t
x2(up,t) = cs*t
a(t) = 2.0/((g+1.0)*t)
b(t) = -2.0*cs/(g+1.0)

vs(up,x,t) = x < x1(up,t) ? -up \
           : x < x2(up,t) ? a(t)*x+b(t) \
           : 0

rs(up,x,t) = r * (( 1.0 - 0.5*(g-1.0)*abs(vs(up,x,t))/cs )**(2.0/(g-1.0)))
ps(up,x,t) = p * (( 1.0 - 0.5*(g-1.0)*abs(vs(up,x,t))/cs )**(2.0*g/(g-1.0)))

set xlabel "x"
set key right top

set xrange [-10.0:10.0]
set yrange [-2.2:0.5]

set grid

set samples 1024
set colorsequence podo

do for [file in FILES] {
   # location-16-t1-beta0.2-u1.0-bx
   initustart = strstrt( file, "-u" )+2
   inituend   = strstrt( file, "-u" )+4
   
   initu = file[initustart:inituend]
   t = step*dt
   
   title = ""
   
   #title = sprintf( "Comparison, piston model and initial velocity model\nvelocity, t = %.5f", t )
   set title title
   set ylabel "velocity"
   set output outdir . file . "-" . step . "-v" . ext
   set yrange [-1.2:0.5]
   plot file index step u 2:4 w l lw 2 t sprintf( "initial velocity model, u_{init} = -%s", initu ) \
      , vs(initu/2.0,x,t) dt "-" lw 2 t sprintf( "piston model, u_{p} = -%.2f", initu/2.0 )
   
   #title = sprintf( "Comparison, piston model and initial velocity model\npressure, t = %.5f", t )
   set title title
   set ylabel "pressure"
   set output outdir . file . "-" . step . "-p" . ext
   set yrange [0.4:1.2]
   plot file index step u 2:5 w l lw 2 t sprintf( "initial velocity model, u_{init} = -%s", initu ) \
      , ps(initu/2.0,x,t) dt "-" lw 2 t sprintf( "piston model, u_{p} = -%.2f", initu/2.0 )
   
   #title = sprintf( "Comparison, piston model and initial velocity model\ndensity, t = %.5f", t )
   set title title
   set ylabel "density"
   set output outdir . file . "-" . step . "-d" . ext
   set yrange [0.3:1.2]
   plot file index step u 2:3 w l lw 2 t sprintf( "initial velocity model, u_{init} = -%s", initu ) \
      , rs(initu/2.0,x,t) dt "-" lw 2 t sprintf( "piston model, u_{p} = -%.2f", initu/2.0 )
}

set output
