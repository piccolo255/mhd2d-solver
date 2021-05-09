#!/usr/bin/env perl

use feature ':5.10';
use warnings;
use strict;

use File::Path qw( make_path );
use File::Basename qw( dirname );

# debugging
# use Time::HiRes qw( usleep );

# ******************************************************************* CONSTANTS

# field indices
# assumed simulation settings:
#     [output grid]
#     natural      = true
#     conservation = false
#     single file  = false
#     binary       = true
my %field = (
   x     =>  0,
   y     =>  1,
   rho   =>  2,
   u     =>  3,
   v     =>  4,
   w     =>  5,
   p     =>  6,
   bx    =>  7,
   by    =>  8,
   bz    =>  9,
   divb  => 10,
);

# sheet-lobe transition types
use constant {
   TRANSITION_UNSET  => 0,
   TRANSITION_LINEAR => 1,
   TRANSITION_STEP   => 2,
};

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

# ***************************************************************** READ CONFIG

my %config;
while( <> ){
   # skip comment lines and empty lines
   next if( $_ =~ /^#|^\s*$/ );
   # stop processing when DATA segment is reached
   last if( $_ =~ /^DATA\r?$/ );
   
   # read config into hashmap
   my ( $var, $content ) = split( /=/ );
   $config{trim($var)} = trim($content);
}

my $data_dir_mask          = $config{"data directories"};
my $profile_file_mask      = $config{"thinning profile files"};
my $location_file_mask     = $config{"thinning location files"};
my $data_files_regexp_mask = $config{"data files regexp"};

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
my $data_type_unpack = "";
for( $data_type ){
   if( TYPE_FLOAT == $_ ){
      $data_type_name   = "float";
      $data_type_size   = 4;
      $data_type_unpack = "f*";
   } elsif ( TYPE_DOUBLE == $_ ){
      $data_type_name   = "double";
      $data_type_size   = 8;
      $data_type_unpack = "d*";
   } else {
      die;
   }
}

# grid density
my $Lx = $config{"Lx"};
my $Ly = $config{"Ly"};
# time range
my $time_min = $config{"time min"};
my $time_max = $config{"time max"};
# coordinate range
my $x_max = $config{"x max"};
my $x_min = $config{"x min"};
my $y_max = $config{"y max"};
my $y_min = $config{"y min"};

# sheet condition
my $watch_value               = $config{"watch value"};
my $sheet_condition           = $config{"sheet condition"};
   $sheet_condition           =~ s/%v/\$var/g;
my $sheet_condition_eval      = "\$is_sheet = $sheet_condition";
my $sheet_thickness_treshhold = $config{"thinning width max"};

my $transition_type = TRANSITION_UNSET;
for( $config{"transition type"} ){ # linear | step
   if    ( not defined  ) { $transition_type = TRANSITION_STEP;    } # old default
   elsif ( /^linear$/   ) { $transition_type = TRANSITION_LINEAR;  }
   elsif ( /^step$/     ) { $transition_type = TRANSITION_STEP;    }
   else                   { say "CONFIG ERROR: 'transition type' must be one of: linear, step." and die; }
}

# my $whereaxis;
# $whereaxis = 0 if $axis =~ "x";
# $whereaxis = 1 if $axis =~ "y";

# ********************************************************************* PROCESS

while( <> ){
   # skip comment lines and empty lines
   next if( $_ =~ /^#|^\s*$/ );
   
   # tokenize
   my ( $grid_n, $lobe_tau, $beta_name, $u_name ) = split;
   
   # process masks
   my $data_dir = $data_dir_mask;
   $data_dir =~ s/%g/$grid_n/g;
   $data_dir =~ s/%r/$lobe_tau/g;
   $data_dir =~ s/%b/$beta_name/g;
   $data_dir =~ s/%u/$u_name/g;
   my $profile_file = $profile_file_mask;
   $profile_file =~ s/%g/$grid_n/g;
   $profile_file =~ s/%r/$lobe_tau/g;
   $profile_file =~ s/%b/$beta_name/g;
   $profile_file =~ s/%u/$u_name/g;
   $profile_file =~ s/%v/$watch_value/g;
   my $location_file = $location_file_mask;
   $location_file =~ s/%g/$grid_n/g;
   $location_file =~ s/%r/$lobe_tau/g;
   $location_file =~ s/%b/$beta_name/g;
   $location_file =~ s/%u/$u_name/g;
   $location_file =~ s/%v/$watch_value/g;
   
   # my $filemask = "shlf-$grid_n-t$lobe_temp-beta$beta_name-u$u_name-.*\\.dat";
   # my $datadir = "$datadirbase/$datasubdir";
   
   say "Processing directory $data_dir";
   
   # get the list of data files to analyze
   opendir DATADIR, $data_dir or die "$!: $data_dir";
   my @data_files = sort grep { /$data_files_regexp_mask/ && -f "$data_dir/$_" } readdir( DATADIR );
   
   # prepare output file for sheet profile data
   my $profile_file_dir = dirname( $profile_file );
   say "Creating profile directory $profile_file_dir...";
   make_path $profile_file_dir unless( -d $profile_file_dir );
   open my $PROFILE_FILE, '>', "$profile_file" or die "Cannot open file '$profile_file' for writing";
   say $PROFILE_FILE "# 1) time  2) x  3) sheet lower y  4) sheet upper y";
   
   # prepare output file for sheet thinning location data
   my $location_file_dir = dirname( $location_file );
   say "Creating location directory $location_file_dir...";
   make_path $location_file_dir unless( -d $location_file_dir );
   open my $LOCATION_FILE, '>', "$location_file" or die "Cannot open file '$location_file' for writing";
   say $LOCATION_FILE "# 1) time  2) x";
   
   # my $Nx = $Lx*$grid_n;
   # my $Ny = $Ly*$grid_n;
   # my $xmod = 1.0/$grid_n/2.0;
   # my $ymod = 1.0/$grid_n/2.0;
   
   # record size in bytes; assume 11 fields
   my $record_size = 11 * $data_type_size;
   
   # turn on auto flushing
   $| = 1;
   print "> ";
   my $clear_line = "";
   my $is_first_record = 1;
   foreach my $filename ( @data_files ){
      # get simulation time from the filename
      my ( $time ) = ( $filename =~ /.*-(\d+\.\d+)\.dat/ );
      
      # skip if there was an error, or time is out of range
      if( !defined $time ){
         print "\r  \r";
         say "Warning: unable to extract time from file: $filename";
         print "> ";
         next;
      }
      next if $time < $time_min;
      next if $time > $time_max;
      
      # print status to terminal
      my $message_line = "$filename: t = $time";
      print "\r> $clear_line";
      print "\r> $message_line";
      $clear_line = $message_line;
      $clear_line =~ s/./ /g;
      
      # separate records
      if( $is_first_record ){
         $is_first_record = 0;
      } else {
         print $PROFILE_FILE "\n\n";
      }
      
      open my $DATAFILE, '<:raw', "$data_dir/$filename" or die "Cannot open file '$data_dir/$filename' for reading";
      
      my $record;
      my $record_count;
      my $sheet_y_min;
      my $sheet_y_max;
      my $is_sheet_thinning = 1;
      my $thinning_front_x = $x_min;
      my $x_prev = $x_min - 1.0;
      my $y_prev = $y_min;
      my $var_prev = 0;
      my $sheet_thickness_prev = 0;
      while( my $nbytes = read( $DATAFILE, $record, $record_size ) ){
         # unpack record
         my @values = unpack $data_type_unpack, $record;
         $record_count++;
         
         # read values from record
         my $x   = $values[$field{x}];
         my $y   = $values[$field{y}];
         my $var = $values[$field{$watch_value}];
         
         # check if record is in range under consideration
         next if $x < $x_min;
         next if $x > $x_max;
         next if $y < $y_min;
         next if $y > $y_max;
         
         # prepare on new scanline
         if( $x != $x_prev ){
            # print "\n\r                      sheet is from $sheet_y_min to $sheet_y_max";
            # print "\rat x = $x_prev";
            if( $x_prev >= $x_min && $x_prev <= $x_max ){
               printf $PROFILE_FILE "%.7e\t%+.7e\t%+.7e\t%+.7e\n", $time, $x_prev, $sheet_y_min, $sheet_y_max;
               
               if( $is_sheet_thinning ){
                  my $sheet_thickness;
                  for( $transition_type ){
                     if( TRANSITION_STEP == $_ ){
                        $sheet_thickness = $sheet_y_max - $sheet_y_min;
                        if( $sheet_thickness <= $sheet_thickness_treshhold ){
                           $thinning_front_x = $x;
                        } else {
                           $is_sheet_thinning = 0;
                        }
                     } elsif( TRANSITION_LINEAR == $_ ){
                        $sheet_thickness = abs(2.0*$sheet_y_min);
                        if( $sheet_thickness <= $sheet_thickness_treshhold ){
                           $thinning_front_x = $x;
                        } else {
                           my $before_threshold = $sheet_thickness_treshhold - $sheet_thickness_prev;
                           my $after_threshold  = $sheet_thickness - $sheet_thickness_treshhold;
                           $thinning_front_x = ( $before_threshold*$x + $after_threshold*$x_prev ) / ($before_threshold+$after_threshold);
                           
                           $is_sheet_thinning = 0;
                        }
                     } else {
                        die;
                     }
                  }
                  $sheet_thickness_prev = $sheet_thickness;
               }
            }
            $sheet_y_min = $y_max;
            $sheet_y_max = $y_min;
            $x_prev = $x;
         }

         # detect sheet
         my $is_sheet;
         eval $sheet_condition_eval;
         if( $is_sheet ){
            for( $transition_type ){
               if( TRANSITION_STEP == $_ ){
                  $sheet_y_min = $y if $sheet_y_min > $y;
                  $sheet_y_max = $y+(1.0/$grid_n);
               } elsif( TRANSITION_LINEAR == $_ ){
                  if( $sheet_y_min > $y ){ # entered the sheet
                     my $realvar = $var; # push backup
                     
                     my $varmin = $var_prev;
                     my $varmax = $var;
                     my $ymin = $y-(1.0/$grid_n);
                     my $ymax = $y;
                     
                     for my $i (0..9) { # binary search
                        $var = ($varmin+$varmax) / 2.0;
                        eval $sheet_condition_eval;
                        if( $is_sheet ){
                           # midpoint IS in the sheet, so move max towards min
                           $varmax = $var;
                           $ymax = ($ymin+$ymax) / 2.0;
                        } else {
                           # midpoint IS NOT in the sheet, so move min towards max
                           $varmin = $var;
                           $ymin = ($ymin+$ymax) / 2.0;
                        }
                     }
                     
                     $sheet_y_min = ($ymin+$ymax) / 2.0;
                     $var = $realvar; # pop backup
                  }
                  $sheet_y_max = $y+(1.0/$grid_n); # TODO
               } else {
                  die;
               }
            }
            # print "\nsheet, var = $var, x = $x, y = $y";
         } else {
            # print "\nlobe,  var = $var, x = $x, y = $y";
         }
         
         $record = '';
         
         $y_prev     = $y;
         $var_prev   = $var;
      }
      
      # say "\n$record_count records found.";
      printf $LOCATION_FILE "%.7e\t%+.7e\n", $time, $thinning_front_x;
      close $DATAFILE;
   }
   print "\r  $clear_line\r";
   
   close $PROFILE_FILE;
   close $LOCATION_FILE;
   closedir DATADIR;
}
