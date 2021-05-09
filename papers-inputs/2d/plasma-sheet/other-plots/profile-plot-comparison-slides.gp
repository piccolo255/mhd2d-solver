FILES1 = "\
profile-32-t2-beta0.2-u1.0-bx.dat \
profile-32-t2-beta0.2-u1.0-bx.dat \
profile-32-t2-beta0.2-u1.0-bx.dat \
"
FILES2 = "\
profile-32-t2-beta0.2-u0.5-bx.dat \
profile-32-t2-beta1.0-u1.0-bx.dat \
profile-32-t5-beta0.2-u1.0-bx.dat \
"
# FILES1 = system("dir /b profile-32*u1.0*.dat")
# FILES2 = system("dir /b profile-32*u0.5*.dat")

set term png size 600,400 enhanced lw 3 font ",14"

set xlabel "x" offset 0,0.5
set ylabel "y" offset 2
set key right top

set xrange [0.0:8.0]
set yrange [-0.8:0.8]

set grid

do for [ii = 1:words(FILES1)] {
   file1 = word(FILES1,ii)
   file2 = word(FILES2,ii)
   
   # profile-16-t1-beta0.2-u1.0-bx
   tempstart1 = strstrt( file1, "-t" )+2
   tempend1   = strstrt( file1, "-beta" )-1
   betastart1 = strstrt( file1, "-beta" )+5
   betaend1   = strstrt( file1, "-u" )-1
   initustart1 = strstrt( file1, "-u" )+2
   inituend1   = strstrt( file1, "-u" )+4
   tempstart2 = strstrt( file2, "-t" )+2
   tempend2   = strstrt( file2, "-beta" )-1
   betastart2 = strstrt( file2, "-beta" )+5
   betaend2   = strstrt( file2, "-u" )-1
   initustart2 = strstrt( file2, "-u" )+2
   inituend2   = strstrt( file2, "-u" )+4
   
   temp1 = file1[tempstart1:tempend1]
   beta1 = file1[betastart1:betaend1]
   initu1 = file1[initustart1:inituend1]
   temp2 = file2[tempstart2:tempend2]
   beta2 = file2[betastart2:betaend2]
   initu2 = file2[initustart2:inituend2]
   
   locationfile1 = "location" . file1[8:]
   locationfile2 = "location" . file2[8:]
   #info = sprintf( "T_{ratio} = %s, \beta = %s, u_{0} = %s", temp, beta, initu )
   info1 = sprintf( "T_S / T_L = %s, β_L = %s, u_{0} = %s", temp1, beta1, initu1 )
   info2 = sprintf( "T_S / T_L = %s, β_L = %s, u_{0} = %s", temp2, beta2, initu2 )
   main_title = "Plasma sheet profile comparison"
   
   dirname = "slides-comp-" . file1 . "-vs-" . file2 . "-plot"
   system( "mkdir ".dirname )
   
   stats file1 using 1
   
   #do for [n = 0:STATS_blocks-1:10] {
   do for [n = 160:320:160] {
      print sprintf( "Plotting record #%d of %d", n, STATS_blocks-1 )
      set output
      plot file1 index n u 1:(time=$1)
      
      time_line = sprintf( "t = %09.6f", time )
      set title main_title . "\n" . time_line offset 0,-0.5
      
      record_id = sprintf( "-%04d-%010.7f", n, time )
      
      #set output dirname . "/" . file1 . "-vs-" . file2 . record_id . ".png"
      set output dirname . "/" . file1 . record_id . ".png"
      plot file1 index n u 2:4 w l lt 1 t info1 \
         , file1 index n u 2:3 w l lt 1 notitle \
         , file2 index n u 2:4 w l lt 2 t info2 \
         , file2 index n u 2:3 w l lt 2 notitle
   }
}

set output
