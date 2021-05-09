FILES = "\
velocity-32-sort-tau.dat \
velocity-32-sort-tau-split.dat \
"
#FILES = system("dir /b velocity-*.dat")

get_alfven( tr, beta ) = \
   ( tr ==  1 ) ? ( \
      ( beta == 0.2 ) ? 3.30200 : \
      ( beta == 0.4 ) ? 2.26779 : \
      ( beta == 0.7 ) ? 1.75023 : \
      ( beta == 1.0 ) ? 1.41421 : \
      1/0 \
   ) : ( tr ==  2 ) ? ( \
      ( beta == 0.2 ) ? 2.33487 : \
      ( beta == 0.4 ) ? 1.60357 : \
      ( beta == 0.7 ) ? 1.23760 : \
      ( beta == 1.0 ) ? 1.00000 : \
      1/0 \
   ) : ( tr ==  5 ) ? ( \
      ( beta == 0.2 ) ? 1.47670 : \
      ( beta == 0.4 ) ? 1.01419 : \
      ( beta == 0.7 ) ? 0.78272 : \
      ( beta == 1.0 ) ? 0.63246 : \
      1/0 \
   ) : ( tr == 10 ) ? ( \
      ( beta == 0.2 ) ? 1.04419 : \
      ( beta == 0.4 ) ? 0.71714 : \
      ( beta == 0.7 ) ? 0.55347 : \
      ( beta == 1.0 ) ? 0.44721 : \
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

set grid

do for [file in FILES] {
   # ********************************** velocity dependence on beta
   title = "Thinning front velocity\nvs plasma beta"
   set title title
   
   set output file . "-dependency-beta.png"
   set xlabel "beta"
   set ylabel "front velocity"
   set xrange [0.1:200.0]
   set yrange [*:*]
   
   plot file u 3:($2==1.0?$7:1/0) w lp t "tau = 1" \
      , file u 3:($2==2.0?$7:1/0) w lp t "tau = 2" \
      , file u 3:($2==5.0?$7:1/0) w lp t "tau = 5"
   
   # ********************************** velocity dependence on alfven
   title = "Thinning front velocity\nvs lobe Alfven velocity"
   set title title
   
   set output file . "-dependency-alfven.png"
   set xlabel "alfven velocity"
   set ylabel "front velocity"
   set xrange [0.02:4.0]
   set yrange [*:*]
   set key right bottom
   
   plot file u 5:($2==1.0?$7:1/0) w lp t "tau = 1" \
      , file u 5:($2==2.0?$7:1/0) w lp t "tau = 2" \
      , file u 5:($2==5.0?$7:1/0) w lp t "tau = 5"
   
   # ********************************** velocity dependence on alfven, logscale
   title = "Thinning front velocity\nvs lobe Alfven velocity"
   set title title
   
   set output file . "-dependency-alfven-logscale.png"
   set xlabel "alfven velocity"
   set ylabel "front velocity"
   set xrange [0.02:4.0]
   set yrange [*:*]
   set key right bottom
   set logscale x
   
   plot file u 5:($2==1.0?$7:1/0) w lp t "tau = 1" \
      , file u 5:($2==2.0?$7:1/0) w lp t "tau = 2" \
      , file u 5:($2==5.0?$7:1/0) w lp t "tau = 5"
   unset logscale x
   
   # ********************************** velocity dependence on bx
   title = "Thinning front velocity\nvs B_x"
   set title title
   
   set output file . "-dependency-bx.png"
   set xlabel "B_x"
   set ylabel "front velocity"
   set xrange [0.05:1.35]
   set yrange [*:*]
   set key right bottom
   
   plot file u 4:($2==1.0?$7:1/0) w lp t "tau = 1" \
      , file u 4:($2==2.0?$7:1/0) w lp t "tau = 2" \
      , file u 4:($2==5.0?$7:1/0) w lp t "tau = 5"
   
   # ********************************** velocity dependence on tau
   title = "Thinning front velocity\nvs tau"
   set title title
   
   set output file . "-dependency-tau.png"
   set xlabel "tau"
   set ylabel "front velocity"
   set xrange [0.5:7.5]
   set yrange [*:*]
   set key right bottom
   
   plot file u 2:($3==199.0?$7:1/0) w lp t "beta = 199.0" \
      , file u 2:($3==31.0?$7:1/0) w lp t "beta = 31.0" \
      , file u 2:($3==7.0?$7:1/0) w lp t "beta = 7.0" \
      , file u 2:($3==1.0?$7:1/0) w lp t "beta = 1.0" \
      , file u 2:($3==0.7?$7:1/0) w lp t "beta = 0.7" \
      , file u 2:($3==0.4?$7:1/0) w lp t "beta = 0.4" \
      , file u 2:($3==0.2?$7:1/0) w lp t "beta = 0.2"
   
   set logscale x
   set output file . "-dependency-logscale-tau.png"
   replot
   unset logscale
   
   # ********************************** tau-beta scatter plot
   title = "Thinning front velocity\ntau vs beta"
   set title title
   
   set output file . "-dependency-tau-beta-scatter.png"
   set xlabel "tau"
   set ylabel "beta"
   set cblabel "front velocity"
   set xrange [0.5:5.5]
   set yrange [0.1:200.0]
   set key right bottom
   
   load "rdbu.pal"
   
   plot file u 2:3:7 w p ps 1 lt 4 lw 12 lc palette notitle
   
   set xrange [0.5:5.5]
   set logscale x
   set output file . "-dependency-logscale-tau-beta-scatter.png"
   replot
   unset logscale
   
   # ********************************** tau-beta surface plot, gaussian
   title = "Thinning front velocity\ntau vs beta"
   set title title
   
   set output file . "-dependency-tau-beta-surface.png"
   set xlabel "tau"
   set ylabel "beta"
   set cblabel "front velocity"
   set xrange [*:*]
   set yrange [*:*]
   set key right bottom
   
   load "rdbu.pal"
   
   set dgrid3d 10,10 gauss 10e-2,10e-2
   set view map
   
   splot file u 2:3:7 w pm3d notitle
   
   set logscale x
   set output file . "-dependency-logscale-tau-beta-surface.png"
   replot
   unset logscale
   
   unset view
   unset dgrid3d
   unset logscale
   
   # ********************************** tau-beta contour plot, gaussian
   title = "Thinning front velocity\ntau vs beta"
   set title title
   
   set output file . "-dependency-tau-beta-surface-contour.png"
   set xlabel "tau"
   set ylabel "beta"
   set cblabel "front velocity"
   set xrange [*:*]
   set yrange [*:*]
   set key right bottom
   
   load "rdbu.pal"
   
   set dgrid3d 10,10 gauss 10e-2,10e-2
   set view map
   unset surface
   set contour
   set cntrparam levels incremental 0.35,0.10,0.85
   
   splot file u 2:3:7 w pm3d notitle
   
   set logscale x
   set output file . "-dependency-logscale-tau-beta-surface-contour.png"
   replot
   unset logscale
   
   unset view
   unset dgrid3d
   unset contour
   set surface implicit
   
   # ********************************** sound-alfven scatter plot
   title = "Thinning front velocity\nlobe sound velocity vs lobe Alfven velocity"
   set title title
   
   set output file . "-dependency-sound-alfven-scatter.png"
   set xlabel "lobe sound velocity"
   set ylabel "lobe Alfven velocity"
   set cblabel "front velocity"
   set xrange [*:*]
   set yrange [*:*]
   set key right bottom
   
   load "rdbu.pal"
   
   plot file u 6:5:7 w p ps 1 lt 4 lw 12 lc palette notitle
   
   set logscale x
   set output file . "-dependency-logscale-sound-alfven-scatter.png"
   replot
   unset logscale
   
   # ********************************** sound-alfven surface plot, gaussian
   title = "Thinning front velocity\nlobe sound velocity vs lobe Alfven velocity"
   set title title
   
   set output file . "-dependency-sound-alfven-surface.png"
   set xlabel "lobe sound velocity"
   set ylabel "lobe Alfven velocity"
   set cblabel "front velocity"
   set xrange [*:*]
   set yrange [*:*]
   set key right bottom
   
   load "rdbu.pal"
   
   set dgrid3d 10,10 gauss 10e-2,10e-2
   set view map
   
   splot file u 6:5:7 w pm3d notitle
   
   set logscale x
   set output file . "-dependency-logscale-sound-alfven-surface.png"
   replot
   unset logscale
   
   unset view
   unset dgrid3d
   unset logscale
   
   # ********************************** sound-alfven contour plot, gaussian
   title = "Thinning front velocity\nlobe sound velocity vs lobe Alfven velocity"
   set title title
   
   set output file . "-dependency-sound-alfven-surface-contour.png"
   set xlabel "lobe sound velocity"
   set ylabel "lobe Alfven velocity"
   set cblabel "front velocity"
   set xrange [*:*]
   set yrange [*:*]
   set key right bottom
   
   load "rdbu.pal"
   
   set dgrid3d 10,10 gauss 10e-2,10e-2
   set view map
   unset surface
   set contour
   set cntrparam levels incremental 0.35,0.10,0.85
   
   splot file u 6:5:7 w pm3d notitle
   
   set logscale x
   set output file . "-dependency-logscale-sound-alfven-surface-contour.png"
   replot
   unset logscale
   
   unset view
   unset dgrid3d
   unset contour
   set surface implicit
}

set output
