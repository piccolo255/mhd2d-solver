#!/usr/bin/env perl

use feature ':5.10';
use warnings;
use strict;

use File::Path qw(make_path);
use Getopt::Long qw(GetOptions);

Getopt::Long::Configure qw(gnu_getopt);

# debugging
# use Time::HiRes qw( usleep );

# ******************************************************************* CONSTANTS

# assumed simulation settings:
#     [output grid]
#     natural      = true
#     conservation = false
#     single file  = false
#     binary       = true
my %plotexprs = ( x    => "( \$1)",
                  y    => "( \$2)",
                  d    => "( \$3)",
                  u    => "( \$4)",
                  v    => "( \$5)",
                  w    => "( \$6)",
                  p    => "( \$7)",
                  bx   => "( \$8)",
                  by   => "( \$9)",
                  bz   => "(\$10)",
                  divb => "(\$11)",
                  b    => "(sqrt(\$8**2+\$9**2+\$10**2))",
                  vtot => "(sqrt(\$4**2+\$5**2+\$6**2))",
                  ptot => "(\$7+0.5*(\$8**2+\$9**2+\$10**2))",
                  vvec => "(scalex*\$4):(scaley*\$5)",
                  bvec => "(scalex*\$8):(scaley*\$9)",
                  mu   => "(\$3*\$4)",
                  mv   => "(\$3*\$5)" );

my %plotlabels = ( x    => "x",
                   y    => "y",
                   d    => "density",
                   u    => "u",
                   v    => "v",
                   w    => "w",
                   p    => "pressure",
                   bx   => "B_x",
                   by   => "B_y",
                   bz   => "B_z",
                   divb => "div B",
                   b    => "B",
                   vtot => "velocity",
                   ptot => "total pressure",
                   vvec => "velocity vector",
                   bvec => "magnetic vector",
                   mu   => "mass flow in x direction, rho u",
                   mv   => "mass flow in y direction, rho v" );

# data field types
use constant {
   TYPE_UNSET  => 0,
   TYPE_FLOAT  => 1,
   TYPE_DOUBLE => 2,
};

# ******************************************************************* FUNCTIONS

# trim functions from https://perlmaven.com/trim
sub ltrim { my $s = shift; $s =~ s/^\s+//;       return $s };
sub rtrim { my $s = shift; $s =~ s/\s+$//;       return $s };
sub  trim { my $s = shift; $s =~ s/^\s+|\s+$//g; return $s };

# ************************************************************* READ PARAMETERS

my $dump_plot_script;
GetOptions(
    'dump|d' => \$dump_plot_script,
) or die "Usage: $0 [--dump|-d] config_file\n";

# ***************************************************************** READ CONFIG

my %config;
while( <> ){
   # skip comment lines and empty lines
   next if( $_ =~ /^#.*|^\s*$/ );
   # stop processing when DATA segment is reached
   last if( $_ =~ /^DATA\r?$/ );
   
   # read config into hashmap
   my ( $var, $content ) = split( /=/ );
   die "config section needs to be of form 'key = value'" unless defined $var && length $var > 0;
   die "config section needs to be of form 'key = value'" unless defined $content && length $content > 0;
   $config{trim($var)} = trim($content);
}

my $data_dir_mask          = $config{"data directory"};
my $data_files_regexp_mask = $config{"data files regexp"};
my $plot_dir_mask          = $config{"plot directory"};

# data field type
my $data_type = TYPE_UNSET;
for( $config{"data type"} ){ # float | double
   if    ( not defined  ) { $data_type = TYPE_FLOAT;   } # old default
   elsif ( /^float$/    ) { $data_type = TYPE_FLOAT;   }
   elsif ( /^double$/   ) { $data_type = TYPE_DOUBLE;  }
   else                   { say "CONFIG ERROR: 'data type' must be one of: float, double." and die; }
}
my $data_type_name   = "";
my $data_type_size   = 0;
for( $data_type ){
   if( TYPE_FLOAT == $_ ){
      $data_type_name   = "float";
      $data_type_size   = 4;
   } elsif ( TYPE_DOUBLE == $_ ){
      $data_type_name   = "double";
      $data_type_size   = 8;
   } else {
      die;
   }
}

# grid
my $Lx = $config{"Lx"};
my $Ly = $config{"Ly"};

# ********************************************************************* PROCESS

while( <> ){
   # skip comment lines and empty lines
   next if( $_ =~ /^#|^\s*$/ );
   
   # tokenize
   my ( $grid_nx, $grid_ny, $lobe_tau, $beta_name, $u_name, $variable, $cbmin, $cbmax, $xmin, $xmax, $ymin, $ymax ) = split;
   
   # default values -> autoscale
   $cbmin = "*" unless defined $cbmin && length $cbmin > 0;
   $cbmax = "*" unless defined $cbmax && length $cbmax > 0;
   $xmin = "*" unless defined $xmin && length $xmin > 0;
   $xmax = "*" unless defined $xmax && length $xmax > 0;
   $ymin = "*" unless defined $ymin && length $ymin > 0;
   $ymax = "*" unless defined $ymax && length $ymax > 0;
   
   # process masks
   my $data_dir = $data_dir_mask;
   $data_dir =~ s/%gx/$grid_nx/g;
   $data_dir =~ s/%gy/$grid_ny/g;
   $data_dir =~ s/%r/$lobe_tau/g;
   $data_dir =~ s/%b/$beta_name/g;
   $data_dir =~ s/%u/$u_name/g;
   $data_dir =~ s/%v/$variable/g;
   my $plot_dir = $plot_dir_mask;
   $plot_dir =~ s/%gx/$grid_nx/g;
   $plot_dir =~ s/%gy/$grid_ny/g;
   $plot_dir =~ s/%r/$lobe_tau/g;
   $plot_dir =~ s/%b/$beta_name/g;
   $plot_dir =~ s/%u/$u_name/g;
   $plot_dir =~ s/%v/$variable/g;
   
   say "Processing $data_dir/$data_files_regexp_mask, variable: $variable";
   
   opendir DATADIR, $data_dir or die "$!: $data_dir";
   my @data_files = grep { /$data_files_regexp_mask/ && -f "$data_dir/$_" } readdir( DATADIR );
   # say "Creating directory $plot_dir...";
   make_path $plot_dir unless( -d $plot_dir );
   
   my $Nx = $Lx*$grid_nx;
   my $Ny = $Ly*$grid_ny;
   my $xmod = 1.0/$grid_nx/2.0;
   my $ymod = 1.0/$grid_ny/2.0;
   
   my $GNUPLOT;
   if( not $dump_plot_script ){
      open $GNUPLOT, '|-', 'gnuplot' or die "Couldn't pipe to gnuplot: $!";
   }
   
   foreach my $filename ( @data_files ){
      my ( $time ) = ( $filename =~ /.*-(\d+\.\d+)\.dat/ );
      next unless defined $time;
      
      if( $dump_plot_script ){
         open $GNUPLOT, '>', "$plot_dir/$filename-$variable.gp" or die "Couldn't open output plot script: $!";
      }
      
      say $GNUPLOT "set term png size 600,300 enhanced";
      #say $GNUPLOT "set term png size 2388,500 enhanced";
      say $GNUPLOT "set output \"$plot_dir/$filename-$variable.png\"";
      
      say $GNUPLOT "set xrange [$xmin:$xmax]";
      say $GNUPLOT "set yrange [$ymin:$ymax]";
      say $GNUPLOT "set cbrange [$cbmin:$cbmax]";
      #say $GNUPLOT "load \"jet.pal\"";
      
      say $GNUPLOT "set title \"{/*0.8 grid density = ($grid_nx,$grid_ny), T_S / T_L = $lobe_tau, β_L = $beta_name, u_0 = $u_name}\\nt = $time\"";
      say $GNUPLOT "set cblabel \"$plotlabels{$variable}\"";
      
      if( $variable =~ /.vec/ ){
         # vector plot
         say $GNUPLOT "set term png size 3000,1500 lw 3 font \",72\" enhanced";
         say $GNUPLOT "scalex = $cbmin";
         say $GNUPLOT "scaley = $cbmax";
         say $GNUPLOT "density = 1";
         say $GNUPLOT "plot \"$data_dir/$filename\" binary record=($Nx,$Ny) skip=0 format=\"%11$data_type_name\" \\";
         say $GNUPLOT " u (\$1+$xmod):(\$2+$ymod):$plotexprs{$variable} w vectors t \"\"";
      } else {
         #scalar plot
         say $GNUPLOT "plot \"$data_dir/$filename\" binary record=($Nx,$Ny) skip=0 format=\"%11$data_type_name\" \\";
         say $GNUPLOT " u (\$1+$xmod):(\$2+$ymod):$plotexprs{$variable} w image t \"\"";
      }
      
      say $GNUPLOT "set output";
      
      if( $dump_plot_script ){
         close $GNUPLOT;
      }
   }
   
   if( not $dump_plot_script ){
      close $GNUPLOT;
   }
   
   closedir DATADIR;
}
