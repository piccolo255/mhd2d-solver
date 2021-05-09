FILES = "\
profile-32-t2-beta0.2-u1.0-bx.dat \
profile-32-t5-beta1.0-u1.0-bx.dat \
profile-16-t5-beta1.0-u1.0-bx.dat \
"
# FILES = system("dir /b profile-32*.dat")

set term png size 600,400 enhanced lw 3 font ",14"

set xlabel "x" offset 0,0.5
set ylabel "y" offset 2
set key right top

set xrange [0.0:8.0]
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
   
   locationfile = "location" . file[8:]
   #info = sprintf( "T_{ratio} = %s, \beta = %s, u_{0} = %s", temp, beta, initu )
   info = sprintf( "T_S / T_L = %s, β_L = %s, u_{0} = %s", temp, beta, initu )
   main_title = "Plasma sheet profile, " . info
   
   dirname = "slides-" . file . "-plot"
   system( "mkdir ".dirname )
   
   stats file using 1
   
   #do for [n = 0:STATS_blocks-1:10] {
   do for [n = 160:320:160] {
      print sprintf( "Plotting record #%d of %d", n, STATS_blocks-1 )
      set output
      plot file index n u 1:(time=$1)
      
      time_line = sprintf( "t = %09.6f", time )
      set title main_title . "\n" . time_line offset 0,-0.5
      
      record_id = sprintf( "-%04d-%010.7f", n, time )
      
      set output dirname . "/" . file . record_id . ".png"
      plot file index n u 2:4 w l lt 1 t "sheet border" \
         , file index n u 2:3 w l lt 1 notitle \
         , locationfile u 2:($1==time?10:1/0) w impulses lw 2 lc 'red' notitle \
         , locationfile u 2:($1==time?-10:1/0) w impulses lw 2 lc 'red' t "detected front location"
   }
}

set output
