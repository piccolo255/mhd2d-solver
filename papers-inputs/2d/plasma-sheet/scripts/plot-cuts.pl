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
                  ptot => "(\$7+0.5*(\$8**2+\$9**2+\$10**2))" );

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
                   ptot => "total pressure" );

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
my $startx = $config{"start x"};
my $starty = $config{"start y"};

# ********************************************************************* PROCESS

while( <> ){
   # skip comment lines and empty lines
   next if( $_ =~ /^#|^\s*$/ );
   
   # tokenize
   my ( $name, $varcount, $varx, $grid_nx, $grid_ny, $lobe_tau, $beta_name, $u_name, $xmin, $xmax, $ymin, $ymax ) = split;
   
   # default values
   $xmin = -4.0 unless defined $xmin && length $xmin > 0;
   $xmax =  4.0 unless defined $xmax && length $xmax > 0;
   $ymin =  0.0 unless defined $ymin && length $ymin > 0;
   $ymax =  1.0 unless defined $ymax && length $ymax > 0;
   
   # process masks
   my $data_dir = $data_dir_mask;
   $data_dir =~ s/%gx/$grid_nx/g;
   $data_dir =~ s/%gy/$grid_ny/g;
   $data_dir =~ s/%r/$lobe_tau/g;
   $data_dir =~ s/%b/$beta_name/g;
   $data_dir =~ s/%u/$u_name/g;
   $data_dir =~ s/%v/$name/g;
   my $plot_dir = $plot_dir_mask;
   $plot_dir =~ s/%gx/$grid_nx/g;
   $plot_dir =~ s/%gy/$grid_ny/g;
   $plot_dir =~ s/%r/$lobe_tau/g;
   $plot_dir =~ s/%b/$beta_name/g;
   $plot_dir =~ s/%u/$u_name/g;
   $plot_dir =~ s/%v/$name/g;
   
   say "Processing $data_dir/$data_files_regexp_mask, name: $name";
   
   opendir DATADIR, $data_dir or die "$!: $data_dir";
   my @data_files = grep { /$data_files_regexp_mask/ && -f "$data_dir/$_" } readdir( DATADIR );
   # say "Creating directory $plotdir...";
   make_path $plot_dir unless( -d $plot_dir );
   
   my $Nx = $Lx*$grid_nx;
   my $Ny = $Ly*$grid_ny;
   my $xmod = 1.0/$grid_nx/2.0;
   my $ymod = 1.0/$grid_ny/2.0;
   
   my ( $plotx, $plotstarty, $record, $format, $formatsize, $datavarcount, $where_multiplier, $where_n );
   if( $varx =~ /^x$/ ){
      $plotx  = "(\$1+$xmod)";
      $plotstarty = $starty;
      $record = "$Nx";
      $datavarcount = 11;
      $format = "%@{[$datavarcount*$Ny]}$data_type_name";
      $where_multiplier = 1;
      $where_n = $grid_ny;
   } elsif( $varx =~ /^y$/ ){
      $plotx = "(\$2+$ymod)";
      $plotstarty = $startx;
      $record = "$Ny";
      $datavarcount = 11;
      $format = "%$datavarcount$data_type_name";
      $where_multiplier = $Ny;
      $where_n = $grid_nx;
   } else {
      die "Error:$ARGV: unsupported variable for x axis: $varx";
   }
   
   # open my $GNUPLOT, '|-', 'gnuplot' or die "Couldn't pipe to gnuplot: $!";
   
   my @variable;
   my @posy;
   my @scaling;
   my @offsety;
   for my $ii ( 0..$varcount-1 ){
      my $line;
      while( <> ){
         next if( $_ =~ /^#|^\s*$/ );
         $line = $_;
         last;
      }
      ( $variable[$ii], $posy[$ii], $scaling[$ii], $offsety[$ii] ) = split ' ', $line;
      $scaling[$ii] = "1.0" unless defined $scaling[$ii] && length $scaling[$ii] > 0;
      $offsety[$ii] = "0.0" unless defined $offsety[$ii] && length $offsety[$ii] > 0;
      
   }
   
   # turn on auto flushing
   $| = 1;
   print "> ";
   my $empty = "";
   foreach my $filename ( @data_files ){
      my $GNUPLOT;
      if( $dump_plot_script ){
         open $GNUPLOT, '>', "$plot_dir/$filename-$name.gp" or die "Couldn't open output plot script: $!";
      } else {
         open $GNUPLOT, '|-', 'gnuplot' or die "Couldn't pipe to gnuplot: $!";
      }
      
      print "\r> $empty";
      print "\r> $filename";
      $empty = $filename;
      $empty =~ s/./ /g;
      my ( $time ) = ( $filename =~ /.*-(\d+\.\d+)\.dat/ );
      next if not defined $time;
      
      say $GNUPLOT "set term png size 600,300 enhanced";
      say $GNUPLOT "set output \"$plot_dir/$filename-$name.png\"";
      
      say $GNUPLOT "set xrange [$xmin:$xmax]";
      say $GNUPLOT "set yrange [$ymin:$ymax]";
      say $GNUPLOT "set grid";
      
      say $GNUPLOT "load \"jet.pal\"";
      
      say $GNUPLOT "set title \"{/*0.8 grid density = ($grid_nx,$grid_ny), T_S / T_L = $lobe_tau, β_L = $beta_name, u_0 = $u_name}\\nt = $time\"";
      
      say $GNUPLOT "plot \\";
      
      for my $ii ( 0..$#posy ){
         my $where = int( ($posy[$ii]-$plotstarty) * $where_n * $where_multiplier );
         
         my $block_size = $datavarcount*$data_type_size;
         my $offset = $block_size * $where;
         
         #my $title = "$offsety[$ii]+$scaling[$ii]*$variable[$ii] \@ $posy[$ii]";
         my $title = "$variable[$ii]";
         $title .= "*$scaling[$ii]" if( $scaling[$ii] != 1.0 );
         $title .= "+$offsety[$ii]" if( $offsety[$ii] != 0.0 );
         $title .= " \\\\\@ $posy[$ii]";
         say $GNUPLOT " \"$data_dir/$filename\" binary record=($record) skip=$offset format=\"$format\" \\";
         say $GNUPLOT " u $plotx:($plotexprs{$variable[$ii]}*$scaling[$ii]+$offsety[$ii]) w lines t \"$title\", \\";
      }
      say $GNUPLOT "";
      say $GNUPLOT "set output";
      
      close $GNUPLOT;
   }
   print "\r  $empty\r";
   
   # close $GNUPLOT;
   
   closedir DATADIR;
}

say "Done.";
