up = 1
p0 = 1
r0 = 1
t = 1
cs = 1.29
gamma = 5.0/3.0

xl = (cs-1.0/2.0*(gamma+1.0)*up)*t
xr = cs * t

u(x) = x < xl ? -up : \
       x > xr ? 0 : \
                2.0*x/((gamma+1.0)*t) - 2*cs/(gamma+1.0)

alpha(x) = 1.0-1.0/2.0*(gamma-1.0)*abs(u(x))/cs
r(x) = r0 * (alpha(x)**(2.0/(gamma-1.0)))
p(x) = p0 * (alpha(x)**(2.0*gamma/(gamma-1.0)))

set term pngcairo size 800,400 enhanced lw 3 font "Arial,20"
ext = ".png"
set term postscript eps color enhanced size 6,3 lw 4 font "Arial,30"
ext = ".eps"

set xlabel "x"
set ylabel "(-u),p,{/Symbol r}"

set title "t = 1, u_p = 1, c_s = 1.29"
set title ""

set yrange [-0.1:1.1]
set xrange [-1.6:1.6]

set key top left box opaque

set grid

set samples 1024
set colorsequence podo

set out "1d-piston-model-exact" . ext
plot -u(x) t "velocity" w l lw 2 \
   , p(x) t "pressure"  w l lw 2 dt "-" \
   , r(x) t "density"   w l lt 4 lw 2 dt "." \
   , "dummy.txt" u (-1):(10) w impulses lw 2 lc 'red' notitle \
   , "" u (-1):(-10) w impulses lw 2 lc 'red' notitle

set out
