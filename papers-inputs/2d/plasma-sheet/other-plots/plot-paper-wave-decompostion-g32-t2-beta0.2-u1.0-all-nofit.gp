dir = "data/g32-t2/shlf-32-t2-beta0.2-u1.0.ini/decompostion"
outdir = "plots/waves/paper-fit-wave-shlf-g32-t2-beta0.2-u1.0/"
FILES = "\
waves-shlf-32-t2-beta0.2-u1.0-0004-00.1250000.wdat \
waves-shlf-32-t2-beta0.2-u1.0-0008-00.2500000.wdat \
waves-shlf-32-t2-beta0.2-u1.0-0016-00.5000000.wdat \
waves-shlf-32-t2-beta0.2-u1.0-0024-00.7500000.wdat \
waves-shlf-32-t2-beta0.2-u1.0-0032-01.0000000.wdat \
waves-shlf-32-t2-beta0.2-u1.0-0064-02.0000000.wdat \
waves-shlf-32-t2-beta0.2-u1.0-0128-04.0000000.wdat \
"

outdir_create = outdir
# for windows, replace '/' with '\':
# while( (i = strstrt( outdir_create, "/" )) > 0 ){
#    outdir_create = outdir_create[:i-1] . "\\" . outdir_create[i+1:]
# }
print "creating directory: " . outdir_create
system "mkdir " . outdir_create

g = 32
Lx = 32
Ly =  6

Nx = Lx * g
Ny = Ly * g
xmod = 1.0 / g / 2.0
ymod = 1.0 / g / 2.0

set xrange [-2.0:6.0]
set yrange [-2.0:2.0]
set cbrange [0:0.050]
#set cbrange [0:*]

#set term pngcairo size 1200,600 lw 4 enhanced font "Arial,24"
#ext = ".png"
set term postscript eps color enhanced size 6,3 lw 2 font "Arial,30"
ext = ".eps"

#set lmargin screen 0.12
#set rmargin screen 0.7

load "greys.pal"
set cbtics format "%.3f"
set cblabel "component strength"
set title

do for [ file in FILES ] {
   time = file[strstrt(file,".wdat")-10:]
   
   # replace dots with underscores
   file_under = file
   while( (i = strstrt( file_under, "." )) > 0 ){
      file_under = file_under[:i-1] . "_" . file_under[i+1:]
   }
   
   # x
   col = 3
   
   #set title time . "\nx neg. fast ms."
   set out outdir . file . "-x-1-neg-fast" . ext
   plot dir.file u ($1+xmod):($2+ymod):(abs(column(col  ))) w image notitle
   
   #set title time . "\nx neg. alfven"
   set out outdir . file . "-x-2-neg-alfven" . ext
   plot dir.file u ($1+xmod):($2+ymod):(abs(column(col+1))) w image notitle
   
   #set title time . "\nx neg. slow ms."
   set out outdir . file . "-x-3-neg-slow" . ext
   plot dir.file u ($1+xmod):($2+ymod):(abs(column(col+2))) w image notitle
   
   #set title time . "\nx entropy"
   set out outdir . file . "-x-4-entropy" . ext
   plot dir.file u ($1+xmod):($2+ymod):(abs(column(col+3))) w image notitle
   
   #set title time . "\nx pos. slow ms."
   set out outdir . file . "-x-5-pos-slow" . ext
   plot dir.file u ($1+xmod):($2+ymod):(abs(column(col+4))) w image notitle
   
   #set title time . "\nx pos. alfven"
   set out outdir . file . "-x-6-pos-alfven" . ext
   plot dir.file u ($1+xmod):($2+ymod):(abs(column(col+5))) w image notitle
   
   #set title time . "\nx pos. fast ms."
   set out outdir . file . "-x-7-pos-fast" . ext
   plot dir.file u ($1+xmod):($2+ymod):(abs(column(col+6))) w image notitle
   
   # y
   col = 11
   
   #set title time . "\ny neg. fast ms."
   set out outdir . file . "-y-1-neg-fast" . ext
   plot dir.file u ($1+xmod):($2+ymod):(abs(column(col  ))) w image notitle
   
   #set title time . "\ny neg. alfven"
   set out outdir . file . "-y-2-neg-alfven" . ext
   plot dir.file u ($1+xmod):($2+ymod):(abs(column(col+1))) w image notitle
   
   #set title time . "\ny neg. slow ms."
   set out outdir . file . "-y-3-neg-slow" . ext
   plot dir.file u ($1+xmod):($2+ymod):(abs(column(col+2))) w image notitle
   
   #set title time . "\ny entropy"
   set out outdir . file . "-y-4-entropy" . ext
   plot dir.file u ($1+xmod):($2+ymod):(abs(column(col+3))) w image notitle
   
   #set title time . "\ny pos. slow ms."
   set out outdir . file . "-y-5-pos-slow" . ext
   plot dir.file u ($1+xmod):($2+ymod):(abs(column(col+4))) w image notitle
   
   #set title time . "\ny pos. alfven"
   set out outdir . file . "-y-6-pos-alfven" . ext
   plot dir.file u ($1+xmod):($2+ymod):(abs(column(col+5))) w image notitle
   
   #set title time . "\ny pos. fast ms."
   set out outdir . file . "-y-7-pos-fast" . ext
   plot dir.file u ($1+xmod):($2+ymod):(abs(column(col+6))) w image notitle
}

set out
