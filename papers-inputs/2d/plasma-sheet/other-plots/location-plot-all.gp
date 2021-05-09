#FILES = "\
#location-16-t1-beta0.2-u1.0-bx.dat \
#location-16-t2-beta0.2-u1.0-bx.dat \
#location-16-t5-beta0.2-u1.0-bx.dat \
#location-16-t10-beta0.2-u1.0-bx.dat \
#"
FILES = system("dir /b location-*.dat")

get_alfven( tr, beta ) = \
   ( tr ==  1 ) ? ( \
      ( beta ==   0.2 ) ? 3.30200 : \
      ( beta ==   0.4 ) ? 2.26779 : \
      ( beta ==   0.7 ) ? 1.75023 : \
      ( beta ==   1.0 ) ? 1.41421 : \
      ( beta ==   2.6 ) ? 0.88465 : \
      ( beta ==   7.0 ) ? 0.53452 : \
      ( beta ==  31.0 ) ? 0.25400 : \
      ( beta == 199.0 ) ? 0.10025 : \
      1/0 \
   ) : ( tr ==  2 ) ? ( \
      ( beta ==   0.2 ) ? 2.33487 : \
      ( beta ==   0.4 ) ? 1.60357 : \
      ( beta ==   0.7 ) ? 1.23760 : \
      ( beta ==   1.0 ) ? 1.00000 : \
      ( beta ==   2.6 ) ? 0.62554 : \
      ( beta ==   7.0 ) ? 0.37796 : \
      ( beta ==  31.0 ) ? 0.17961 : \
      ( beta == 199.0 ) ? 0.07089 : \
      1/0 \
   ) : ( tr ==  5 ) ? ( \
      ( beta ==   0.2 ) ? 1.47670 : \
      ( beta ==   0.4 ) ? 1.01419 : \
      ( beta ==   0.7 ) ? 0.78272 : \
      ( beta ==   1.0 ) ? 0.63246 : \
      ( beta ==   2.6 ) ? 0.39563 : \
      ( beta ==   7.0 ) ? 0.23905 : \
      ( beta ==  31.0 ) ? 0.11359 : \
      ( beta == 199.0 ) ? 0.04483 : \
      1/0 \
   ) : ( tr == 10 ) ? ( \
      ( beta ==   0.2 ) ? 1.04419 : \
      ( beta ==   0.4 ) ? 0.71714 : \
      ( beta ==   0.7 ) ? 0.55347 : \
      ( beta ==   1.0 ) ? 0.44721 : \
      ( beta ==   2.6 ) ? 0.27975 : \
      ( beta ==   7.0 ) ? 0.16903 : \
      ( beta ==  31.0 ) ? 0.08032 : \
      ( beta == 199.0 ) ? 0.03170 : \
      1/0 \
   ) : ( 1/0 )

get_sound_lobe( tr, beta ) = \
   ( tr ==  1 ) ? ( 1.29099 ) : \
   ( tr ==  2 ) ? ( 0.91287 ) : \
   ( tr ==  5 ) ? ( 0.57735 ) : \
   ( tr == 10 ) ? ( 0.40825 ) : \
   ( 1/0 )

get_sound_sheet( tr, beta ) = 1.29099

set term png size 600,400 enhanced

set xlabel "time"
set ylabel "x"
set key right top

set xrange [2.0:10.0]
set yrange [0.0:12.0]

set grid

do for [file in FILES] {
   # example filename: location-16-t1-beta0.2-u1.0-bx
   tempstart = strstrt( file, "-t" )+2
   tempend   = strstrt( file, "-beta" )-1
   betastart = strstrt( file, "-beta" )+5
   betaend   = strstrt( file, "-u" )-1
   initustart = strstrt( file, "-u" )+2
   inituend   = strstrt( file, "-u" )+4
   
   temp = file[tempstart:tempend]
   beta = file[betastart:betaend]
   initu = file[initustart:inituend]
   
   alfven = get_alfven( temp, beta )
   soundl = get_sound_lobe( temp, beta )
   sounds = get_sound_sheet( temp, beta )
   
   thinv(x) = fa*x + fb
   fit thinv(x) file u 1:2 via fa,fb
   
   info = sprintf( "T_{ratio} = %s, \beta = %s, u_{init} = %s", temp, beta, initu )
   title = "Thinning front location\n" . info
   set title title
   set output file . ".png"
   
   plot file u 1:2 w l t "thinning front" \
      , thinv(x) t sprintf( "fit, v = %.3f", fa ) \
      , alfven*x t sprintf( "Lobe Alfven, v_A = %.2f", alfven ) \
      , soundl*x t sprintf( "Lobe sound, v_{s,L} = %.2f", soundl ) \
      , sounds*x t sprintf( "Sheet sound, v_{s,S} = %.2f", sounds ) lt 7
}

set output
