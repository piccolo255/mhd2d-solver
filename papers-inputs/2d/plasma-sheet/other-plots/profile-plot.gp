FILES = "\
profile-16-t2-beta0.2-u1.0-bx.dat \
"
#FILES = system("dir /b profile-16*.dat")

floor_n(x, n) = floor(x*10**n)*10.0**(-n)

set term png size 600,400 enhanced

set title "Location of maximum sheet pressure"
set xlabel "x"
set ylabel "y"
set key right top

set xrange [0.0:12.0]
set yrange [-0.8:0.8]

set grid

do for [file in FILES] {
   # profile-16-t1-beta0.2-u1.0-bx
   tempstart = strstrt( file, "-t" )+2
   tempend   = strstrt( file, "-beta" )-1
   betastart = strstrt( file, "-beta" )+5
   betaend   = strstrt( file, "-u" )-1
   initustart = strstrt( file, "-u" )+2
   inituend   = strstrt( file, "-u" )+4
   
   temp = file[tempstart:tempend]
   beta = file[betastart:betaend]
   initu = file[initustart:inituend]
   
   info = sprintf( "T_{ratio} = %s, \beta = %s, u_{init} = %s", temp, beta, initu )
   main_title = "Plasma sheet profile, " . info
   
   dirname = file . "-plot"
   system( "mkdir ".dirname )
   
   stats file using 1
   
   file32 = file[1:strstrt(file,"_16-")] . "32" . file[strstrt(file,"_16-")+3:]
   file64 = file[1:strstrt(file,"_16-")] . "64" . file[strstrt(file,"_16-")+3:]
   
   locationfile = "location" . file[8:]
   locationfile32 = "location" . file32[8:]
   locationfile64 = "location" . file64[8:]
   
   do for [n = 0:STATS_blocks-1:10] {
   #do for [n = 50:50] {
      print sprintf( "Plotting record %d of %d", n, STATS_blocks-1 )
      set output "dummy"
      plot file index n u 1:(time=$1)
      
      time_line = sprintf( "t = %09.6f", time )
      set title main_title . "\n" . time_line
      
      record_id = sprintf( "-%04d-%010.7f", n, time )
      set output dirname . "/" . file . record_id . ".png"
      
      timef = floor_n(time,3)
      plot file index n u 2:4 w l lt 1 t "grid 16" \
         , file index n u 2:3 w l lt 1 notitle \
         , file32 index n u 2:4 w l lt 2 t "grid 32" \
         , file32 index n u 2:3 w l lt 2 notitle \
         , file64 index n u 2:4 w l lt 3 t "grid 64" \
         , file64 index n u 2:3 w l lt 3 notitle \
         , locationfile u 2:(floor_n($1,3)==timef?10:1/0)    w impulses lw 2 lt 1 t "grid 16 front" \
         , locationfile u 2:(floor_n($1,3)==timef?-10:1/0)   w impulses lw 2 lt 1 notitle \
         , locationfile32 u 2:(floor_n($1,3)==timef?10:1/0)  w impulses lw 2 lt 2 t "grid 32 front" \
         , locationfile32 u 2:(floor_n($1,3)==timef?-10:1/0) w impulses lw 2 lt 2 notitle \
         , locationfile64 u 2:(floor_n($1,3)==timef?10:1/0)  w impulses lw 2 lt 3 t "grid 64 front" \
         , locationfile64 u 2:(floor_n($1,3)==timef?-10:1/0) w impulses lw 2 lt 3 notitle
   }
}

set output
