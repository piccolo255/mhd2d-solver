#!/usr/bin/env perl

use feature ':5.10';
use warnings;
use strict;

use Time::HiRes qw( time );
use File::Basename qw( dirname );
use File::Basename qw( basename );
use Getopt::Long qw( GetOptions );

Getopt::Long::Configure qw( gnu_getopt );

# debugging
# use Time::HiRes qw( usleep );

# ******************************************************************* CONSTANTS

# field indices
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

my $gamma   = 5.0/3.0;
my $EPS     = 1e-20;

# ******************************************************************* FUNCTIONS

# trim functions from https://perlmaven.com/trim
sub ltrim { my $s = shift; $s =~ s/^\s+//;       return $s };
sub rtrim { my $s = shift; $s =~ s/\s+$//;       return $s };
sub  trim { my $s = shift; $s =~ s/^\s+|\s+$//g; return $s };

sub EPS_EQUAL { my $a = shift; my $b = shift; return abs($a-$b)<$EPS; };
sub EPS_ZERO  { my $x = shift; return abs($x)<$EPS; };

sub array_to_string_f {
   my $array_ref = shift;
   my @a = @$array_ref;
   
   my $s = sprintf( "%+.12e\t" x @a, @a );
   
   return $s;
};

sub matrix_to_string_f {
   my $matrix_ref = shift;
   my @m = @$matrix_ref;
   
   my $s = "";
   for my $array_ref( @m ){
      my @a = @$array_ref;
      $s .= sprintf( "%+.8e\t" x @a, @a );
   }
   
   return $s;
};

sub array_to_string_mathematica {
   my $array_ref = shift;
   my @a = @$array_ref;
   
   my $s = "{";
   $s .= sprintf( "%+.8e," x @a, @a );
   $s .= "};";
   
   $s =~ s/e/\*10\^/g;
   $s =~ s/,}/}/g;
   
   return $s;
};

sub matrix_to_string_mathematica {
   my $matrix_ref = shift;
   my @m = @$matrix_ref;
   
   my $s = "{";
   for my $array_ref( @m ){
      my @a = @$array_ref;
      $s .= "{";
      $s .= sprintf( "%+.8e," x @a, @a );
      $s .= "},";
   }
   $s .= "};";
   
   $s =~ s/e/\*10\^/g;
   $s =~ s/,}/}/g;
   
   return $s;
};

sub m_dot_v {
   my $matrix_ref = shift;
   my $array_ref = shift;
   my @m = @$matrix_ref;
   my @a = @$array_ref;
   
   my @r;
   #my $s = "";
   for( my $i = 0; $i <= 7; $i++ ){
      $r[$i] = 0.0;
      #$s .= sprintf( "r$i = 0" );
      for( my $j = 0; $j <= 7; $j++ ){
         $r[$i] += $m[$i][$j]*$a[$j];
         #$s .= sprintf( " + %.2e*%.2e", $m[$i][$j], $a[$j] );
      }
      #$s .= sprintf( " = %.2e\n", $r[$i] );
   }
   
   #$s =~ s/e/\*10\^/g;
   #printf $s;
   
   return \@r;
};

sub v_sub_v {
   my $arrayl_ref = shift;
   my $arrayr_ref = shift;
   my @al = @$arrayl_ref;
   my @ar = @$arrayr_ref;
   
   my @r;
   for( my $i = 0; $i <= 7; $i++ ){
      $r[$i] = $al[$i] - $ar[$i];
   }
   
   return \@r;
};

sub transposed_m_dot_v {
   my $matrix_ref = shift;
   my $array_ref = shift;
   my @m = @$matrix_ref;
   my @a = @$array_ref;
   
   my @r;
   #my $s = "";
   for( my $i = 0; $i <= 7; $i++ ){
      $r[$i] = 0.0;
      #$s .= sprintf( "r$i = 0" );
      for( my $j = 0; $j <= 7; $j++ ){
         $r[$i] += $m[$j][$i]*$a[$j];
         #$s .= sprintf( " + %.2e*%.2e", $m[$j][$i], $a[$j] );
      }
      #$s .= sprintf( " = %.2e\n", $r[$i] );
   }
   
   #$s =~ s/e/\*10\^/g;
   #printf $s;
   
   return \@r;
};

sub to_conservation_val {
   my %data = @_;
   
   # natural values
   my $r  = $data{rho};
   my $u  = $data{u  };
   my $v  = $data{v  };
   my $w  = $data{w  };
   my $bx = $data{bx };
   my $by = $data{by };
   my $bz = $data{bz };
   my $p  = $data{p  };
   
   my $uu =  $u*$u   + $v*$v   + $w*$w;
   my $bb = $bx*$bx + $by*$by + $bz*$bz;
   my @U;
   
   # conservation values
   $U[0] = $r;
   $U[1] = $r*$u;
   $U[2] = $r*$v;
   $U[3] = $r*$w;
   $U[4] = $bx;
   $U[5] = $by;
   $U[6] = $bz;
   $U[7] = $p/($gamma-1.0) + 0.5*$r*$uu + 0.5*$bb;
   
   return \@U;
};

sub flux_F {
   my $data_ref = shift;
   my %data = %$data_ref;
   
   my $r  = $data{rho};
   my $u  = $data{u  };
   my $v  = $data{v  };
   my $w  = $data{w  };
   my $bx = $data{bx };
   my $by = $data{by };
   my $bz = $data{bz };
   my $p  = $data{p  };
   
   my $uu =  $u*$u  +  $v*$v  +  $w*$w;
   my $ub =  $u*$bx +  $v*$by +  $w*$bz;
   my $bb = $bx*$bx + $by*$by + $bz*$bz;
   
   my $mx = $r*$u;
   my $my = $r*$v;
   my $mz = $r*$w;
   my $e  = $p/($gamma-1.0) + 0.5*$r*$uu + 0.5*$bb;
   my $ptot = $p + 0.5*$bb;
   
   my @F;
   $F[0] = $mx;                        # rho
   $F[1] = $mx*$u - $bx*$bx + $ptot;   # mx
   $F[2] = $my*$u - $bx*$by;           # my
   $F[3] = $mz*$u - $bx*$bz;           # mz
   $F[4] = 0;                          # bx
   $F[5] = $by*$u - $bx*$v;            # by
   $F[6] = $bz*$u - $bx*$w;            # bz
   $F[7] = ($e+$ptot)*$u - $bx*$ub;    # e
   
   return \@F;
};

sub flux_G {
   my $data_ref = shift;
   my %data = %$data_ref;
   
   my $r  = $data{rho};
   my $u  = $data{u  };
   my $v  = $data{v  };
   my $w  = $data{w  };
   my $bx = $data{bx };
   my $by = $data{by };
   my $bz = $data{bz };
   my $p  = $data{p  };
   
   my $uu =  $u*$u  +  $v*$v  +  $w*$w;
   my $ub =  $u*$bx +  $v*$by +  $w*$bz;
   my $bb = $bx*$bx + $by*$by + $bz*$bz;
   
   my $mx = $r*$u;
   my $my = $r*$v;
   my $mz = $r*$w;
   my $e  = $p/($gamma-1.0) + 0.5*$r*$uu + 0.5*$bb;
   my $ptot = $p + 0.5*$bb;
   
   my @G;
   $G[0] = $my;                        # rho
   $G[1] = $mx*$v - $by*$bx;           # mx
   $G[2] = $my*$v - $by*$by + $ptot;   # my
   $G[3] = $mz*$v - $by*$bz;           # mz
   $G[4] = $bx*$v - $by*$u;            # bx
   $G[5] = 0;                          # by
   $G[6] = $bz*$v - $by*$w;            # bz
   $G[7] = ($e+$ptot)*$v - $by*$ub;    # e
   
   return \@G;
};

sub eigens_at_midpoint {
   my $datal_ref = shift;
   my $datar_ref = shift;
   my %datal = %$datal_ref;
   my %datar = %$datar_ref;
   
   # half-point values
   my $r  = 0.5 * ( $datal{rho} + $datar{rho} );
   my $u  = 0.5 * ( $datal{u  } + $datar{u  } );
   my $v  = 0.5 * ( $datal{v  } + $datar{v  } );
   my $w  = 0.5 * ( $datal{w  } + $datar{w  } );
   my $bx = 0.5 * ( $datal{bx } + $datar{bx } );
   my $by = 0.5 * ( $datal{by } + $datar{by } );
   my $bz = 0.5 * ( $datal{bz } + $datar{bz } );

   # vector inner products
   my $uu =  $u*$u  +  $v*$v  +  $w*$w;
   my $bb = $bx*$bx + $by*$by + $bz*$bz;

   # half-point pressure
   my $mml = ($datal{rho}*$datal{u})**2 + ($datal{rho}*$datal{v})**2 + ($datal{rho}*$datal{w})**2;
   my $mmr = ($datar{rho}*$datar{u})**2 + ($datar{rho}*$datar{v})**2 + ($datar{rho}*$datar{w})**2;
   my $bbl = $datal{bx}**2 + $datal{by}**2 + $datal{bz}**2;
   my $bbr = $datar{bx}**2 + $datar{by}**2 + $datar{bz}**2;
   my $ptotl = $datal{p} + 0.5*$bbl;
   my $ptotr = $datar{p} + 0.5*$bbr;
   my $ptot  = 0.5 * ( $ptotl + $ptotr );
   my $p = $ptot - 0.5*$bb;
   if( $p < 0.0 ){
      $p = 0.0;
   }
   
   my %datam;
   $datam{rho} = $r;
   $datam{u  } = $u;
   $datam{v  } = $v;
   $datam{w  } = $w;
   $datam{bx } = $bx;
   $datam{by } = $by;
   $datam{bz } = $bz;
   $datam{p  } = $p;
   
   return eigens( %datam );
}

sub eigens {
   my %data = @_;
   
   # values
   my $r  = $data{rho};
   my $u  = $data{u  };
   my $v  = $data{v  };
   my $w  = $data{w  };
   my $bx = $data{bx };
   my $by = $data{by };
   my $bz = $data{bz };
   my $p  = $data{p  };
      
   if( $r <= 0.0 ){
      die "non-positive density";
   }
   if( $p <= 0.0 ){
      die "non-positive pressure";
   }
   
   # vector inner products
   my $uu =  $u*$u   + $v*$v   + $w*$w;
   my $bb = $bx*$bx + $by*$by + $bz*$bz;
   
   # speeds
   my $a2  = $gamma*$p / $r;
   if( $a2 < 0.0 ){ $a2 = 0.0; }
   my $a   = sqrt($a2);    # sound speed
   my $ca2 = $bx*$bx/$r;
   my $ca  = sqrt($ca2);   # Alfven speed
   my $c2  = 0.5*($bb/$r+$a2);
   my $ctemp = $c2*$c2-$a2*$ca2; if( $ctemp < 0.0 ){ $ctemp = 0.0; }
   my $cs2 = $c2 - sqrt($ctemp); if( $cs2 < 0.0 ){ $cs2 = $EPS; }
   my $cs  = sqrt($cs2);   # slow magnetosonic speed
   my $cf2 = $c2 + sqrt($ctemp); if( $cf2 < 0.0 ){ $cf2 = $EPS; }
   my $cf  = sqrt($cf2);   # fast magnetosonic speed

   # other
   my $rroot   = sqrt( $r );
   my $bperp2  = $by*$by + $bz*$bz;
   my $sgnbx   = ($bx<0.0) ? -1.0 : 1.0;
   my $alphaf2 = ( ($bperp2<$EPS) && EPS_EQUAL($ca2,$a2) ) ? 1.0 : ($cf2-$ca2) / ($cf2-$cs2);
      $alphaf2 = $alphaf2 < 0.0 ? 0.0 : $alphaf2; # hack for cf ~ ca
   my $alphaf  = sqrt( $alphaf2 );
   my $alphas2 = ( ($bperp2<$EPS) && EPS_EQUAL($ca2,$a2) ) ? 1.0 : ($cf2- $a2) / ($cf2-$cs2);
      $alphas2 = $alphas2 < 0.0 ? 0.0 : $alphas2; # hack for cf ~ a
   my $alphas  = sqrt( $alphas2 );
   my $betay   = $bperp2 > $EPS ? $by / sqrt( $bperp2 ) : 1.0/sqrt(2.0);
   my $betaz   = $bperp2 > $EPS ? $bz / sqrt( $bperp2 ) : 1.0/sqrt(2.0);
   my $theta1  = 1.0 / ( $alphaf2*$a2*($cf2-($gamma-2.0)/($gamma-1.0)*$a2) + $alphas2*$cf2*($cs2-($gamma-2.0)/($gamma-1.0)*$a2) );
   my $theta2  = 1.0 / ( $alphaf2*$cf*$a*$sgnbx + $alphas2*$cs*$ca*$sgnbx );
   
   my @lambda;
   my @rv;
   my @lv;
   
   ### eigenvalues
   $lambda[0] = $u - $cf;
   $lambda[1] = $u - $ca;
   $lambda[2] = $u - $cs;
   $lambda[3] = $u;
   $lambda[4] = $u + $cs;
   $lambda[5] = $u + $ca;
   $lambda[6] = $u + $cf;
   $lambda[7] = 0;
   
   ### right eigenvectors
   # a1 = u - cf (fast magnetosonic)
   $rv[0][0] =   $alphaf;
   $rv[0][1] =   $alphaf * ($u-$cf);
   $rv[0][2] =   $alphaf*$v + $alphas*$betay*$ca*$sgnbx;
   $rv[0][3] =   $alphaf*$w + $alphas*$betaz*$ca*$sgnbx;
   $rv[0][4] = 0;
   $rv[0][5] =   $alphas*$betay*$cf / $rroot;
   $rv[0][6] =   $alphas*$betaz*$cf / $rroot;
   $rv[0][7] = 0.5*$alphaf*$uu + $alphaf*$cf2/($gamma-1.0) - $alphaf*$cf*$u + $alphas*$ca*$sgnbx*($betay*$v+$betaz*$w) + ($gamma-2.0)/($gamma-1.0)*$alphaf*($cf2-$a2);
   
   # a2 = u - ca (Alfven)
   $rv[1][0] = 0;
   $rv[1][1] = 0;
   $rv[1][2] =   $betaz*$sgnbx;
   $rv[1][3] = - $betay*$sgnbx;
   $rv[1][4] = 0;
   $rv[1][5] =   $betaz/$rroot;
   $rv[1][6] = - $betay/$rroot;
   $rv[1][7] =   ($betaz*$v-$betay*$w)*$sgnbx;
   
   # a3 = u - cs (slow magnetosonic)
   $rv[2][0] =   $alphas;
   $rv[2][1] =   $alphas * ($u-$cs);
   $rv[2][2] =   $alphas*$v - $alphaf*$betay*$a*$sgnbx;
   $rv[2][3] =   $alphas*$w - $alphaf*$betaz*$a*$sgnbx;
   $rv[2][4] = 0;
   $rv[2][5] = - $alphaf*$betay*$a2 / ($cf*$rroot);
   $rv[2][6] = - $alphaf*$betaz*$a2 / ($cf*$rroot);
   $rv[2][7] = 0.5*$alphas*$uu + $alphas*$cs2/($gamma-1.0) - $alphas*$cs*$u - $alphaf*$a*$sgnbx*($betay*$v+$betaz*$w) + ($gamma-2.0)/($gamma-1.0)*$alphas*($cs2-$a2);
   
   # a4 = u (entropy)
   $rv[3][0] = 1;
   $rv[3][1] = $u;
   $rv[3][2] = $v;
   $rv[3][3] = $w;
   $rv[3][4] = 0;
   $rv[3][5] = 0;
   $rv[3][6] = 0;
   $rv[3][7] = 0.5*$uu;
   
   # a5 = u + cs (slow magnetosonic)
   $rv[4][0] =   $alphas;
   $rv[4][1] =   $alphas * ($u+$cs);
   $rv[4][2] =   $alphas*$v + $alphaf*$betay*$a*$sgnbx;
   $rv[4][3] =   $alphas*$w + $alphaf*$betaz*$a*$sgnbx;
   $rv[4][4] = 0;
   $rv[4][5] = - $alphaf*$betay*$a2 / ($cf*$rroot);
   $rv[4][6] = - $alphaf*$betaz*$a2 / ($cf*$rroot);
   $rv[4][7] = 0.5*$alphas*$uu + $alphas*$cs2/($gamma-1.0) + $alphas*$cs*$u + $alphaf*$a*$sgnbx*($betay*$v+$betaz*$w) + ($gamma-2.0)/($gamma-1.0)*$alphas*($cs2-$a2);
   
   # a6 = u + ca (Alfven)
   $rv[5][0] = 0;
   $rv[5][1] = 0;
   $rv[5][2] = - $betaz*$sgnbx;
   $rv[5][3] = + $betay*$sgnbx;
   $rv[5][4] = 0;
   $rv[5][5] =   $betaz/$rroot;
   $rv[5][6] = - $betay/$rroot;
   $rv[5][7] = - ($betaz*$v-$betay*$w)*$sgnbx;
   
   # a7 = u + cf (fast magnetosonic)
   $rv[6][0] =   $alphaf;
   $rv[6][1] =   $alphaf * ($u+$cf);
   $rv[6][2] =   $alphaf*$v - $alphas*$betay*$ca*$sgnbx;
   $rv[6][3] =   $alphaf*$w - $alphas*$betaz*$ca*$sgnbx;
   $rv[6][4] = 0;
   $rv[6][5] =   $alphas*$betay*$cf / $rroot;
   $rv[6][6] =   $alphas*$betaz*$cf / $rroot;
   $rv[6][7] = 0.5*$alphaf*$uu + $alphaf*$cf2/($gamma-1.0) + $alphaf*$cf*$u - $alphas*$ca*$sgnbx*($betay*$v+$betaz*$w) + ($gamma-2.0)/($gamma-1.0)*$alphaf*($cf2-$a2);
   
   # a8 = 0 (dummy)
   $rv[7][0] = 0;
   $rv[7][1] = 0;
   $rv[7][2] = 0;
   $rv[7][3] = 0;
   $rv[7][4] = 0;
   $rv[7][5] = 0;
   $rv[7][6] = 0;
   $rv[7][7] = 0;
   
   ### left eigenvectors
   # a1 = u - cf (fast magnetosonic)
   $lv[0][0] = 0.25*$theta1*$alphaf*$a2*$uu + 0.5*$theta2*( $alphaf*$a*$u*$sgnbx - $alphas*$cs*($betay*$v+$betaz*$w) );
   $lv[0][1] = -0.5*$theta1*$alphaf*$a2*$u - 0.5*$theta2*$alphaf*$a*$sgnbx;
   $lv[0][2] = -0.5*$theta1*$alphaf*$a2*$v + 0.5*$theta2*$alphas*$betay*$cs;
   $lv[0][3] = -0.5*$theta1*$alphaf*$a2*$w + 0.5*$theta2*$alphas*$betaz*$cs;
   $lv[0][4] = 0;
   $lv[0][5] =  0.5*$theta1*$alphas*$betay*$cf*$rroot*($cs2-($gamma-2.0)/($gamma-1.0)*$a2);
   $lv[0][6] =  0.5*$theta1*$alphas*$betaz*$cf*$rroot*($cs2-($gamma-2.0)/($gamma-1.0)*$a2);
   $lv[0][7] =  0.5*$theta1*$alphaf*$a2;
   
   # a2 = u - ca (Alfven)
   $lv[1][0] = -0.5*($betaz*$v-$betay*$w)*$sgnbx;
   $lv[1][1] = 0;
   $lv[1][2] =  0.5*$betaz*$sgnbx;
   $lv[1][3] = -0.5*$betay*$sgnbx;
   $lv[1][4] = 0;
   $lv[1][5] =  0.5*$betaz*$rroot;
   $lv[1][6] = -0.5*$betay*$rroot;
   $lv[1][7] = 0;
   
   # a3 = u - cs (slow magnetosonic)
   $lv[2][0] = 0.25*$theta1*$alphas*$cf2*$uu + 0.5*$theta2*( $alphas*$ca*$u*$sgnbx + $alphaf*$cf*($betay*$v+$betaz*$w) );
   $lv[2][1] = -0.5*$theta1*$alphas*$cf2*$u - 0.5*$theta2*$alphas*$ca*$sgnbx;
   $lv[2][2] = -0.5*$theta1*$alphas*$cf2*$v - 0.5*$theta2*$alphaf*$betay*$cf;
   $lv[2][3] = -0.5*$theta1*$alphas*$cf2*$w - 0.5*$theta2*$alphaf*$betaz*$cf;
   $lv[2][4] = 0;
   $lv[2][5] = -0.5*$theta1*$alphaf*$betay*$cf*$rroot*($cf2-($gamma-2.0)/($gamma-1.0)*$a2);
   $lv[2][6] = -0.5*$theta1*$alphaf*$betaz*$cf*$rroot*($cf2-($gamma-2.0)/($gamma-1.0)*$a2);
   $lv[2][7] =  0.5*$theta1*$alphas*$cf2;
   
   # a4 = u (entropy)
   $lv[3][0] = 1.0 - 0.5*$theta1*$uu*($alphaf2*$a2+$alphas2*$cf2);
   $lv[3][1] =  $theta1*($alphaf2*$a2+$alphas2*$cf2)*$u;
   $lv[3][2] =  $theta1*($alphaf2*$a2+$alphas2*$cf2)*$v;
   $lv[3][3] =  $theta1*($alphaf2*$a2+$alphas2*$cf2)*$w;
   $lv[3][4] = 0;
   $lv[3][5] =  $theta1*$alphaf*$alphas*$betay*$cf*$rroot*($cf2-$cs2);
   $lv[3][6] =  $theta1*$alphaf*$alphas*$betaz*$cf*$rroot*($cf2-$cs2);
   $lv[3][7] = -$theta1*($alphaf2*$a2+$alphas2*$cf2);
   
   # a5 = u + cs (slow magnetosonic)
   $lv[4][0] = 0.25*$theta1*$alphas*$cf2*$uu - 0.5*$theta2*( $alphas*$ca*$u*$sgnbx + $alphaf*$cf*($betay*$v+$betaz*$w) );
   $lv[4][1] = -0.5*$theta1*$alphas*$cf2*$u + 0.5*$theta2*$alphas*$ca*$sgnbx;
   $lv[4][2] = -0.5*$theta1*$alphas*$cf2*$v + 0.5*$theta2*$alphaf*$betay*$cf;
   $lv[4][3] = -0.5*$theta1*$alphas*$cf2*$w + 0.5*$theta2*$alphaf*$betaz*$cf;
   $lv[4][4] = 0;
   $lv[4][5] = -0.5*$theta1*$alphaf*$betay*$cf*$rroot*($cf2-($gamma-2.0)/($gamma-1.0)*$a2);
   $lv[4][6] = -0.5*$theta1*$alphaf*$betaz*$cf*$rroot*($cf2-($gamma-2.0)/($gamma-1.0)*$a2);
   $lv[4][7] =  0.5*$theta1*$alphas*$cf2;
   
   # a6 = u + ca (Alfven)
   $lv[5][0] = +0.5*($betaz*$v-$betay*$w)*$sgnbx;
   $lv[5][1] = 0;
   $lv[5][2] = -0.5*$betaz*$sgnbx;
   $lv[5][3] =  0.5*$betay*$sgnbx;
   $lv[5][4] = 0;
   $lv[5][5] =  0.5*$betaz*$rroot;
   $lv[5][6] = -0.5*$betay*$rroot;
   $lv[5][7] = 0;
   
   # a7 = u + cf (fast magnetosonic)
   $lv[6][0] = 0.25*$theta1*$alphaf*$a2*$uu - 0.5*$theta2*( $alphaf*$a*$u*$sgnbx - $alphas*$cs*($betay*$v+$betaz*$w) );
   $lv[6][1] = -0.5*$theta1*$alphaf*$a2*$u + 0.5*$theta2*$alphaf*$a*$sgnbx;
   $lv[6][2] = -0.5*$theta1*$alphaf*$a2*$v - 0.5*$theta2*$alphas*$betay*$cs;
   $lv[6][3] = -0.5*$theta1*$alphaf*$a2*$w - 0.5*$theta2*$alphas*$betaz*$cs;
   $lv[6][4] = 0;
   $lv[6][5] =  0.5*$theta1*$alphas*$betay*$cf*$rroot*($cs2-($gamma-2.0)/($gamma-1.0)*$a2);
   $lv[6][6] =  0.5*$theta1*$alphas*$betaz*$cf*$rroot*($cs2-($gamma-2.0)/($gamma-1.0)*$a2);
   $lv[6][7] =  0.5*$theta1*$alphaf*$a2;
   
   # a8 = 0 (dummy)
   $lv[7][0] = 0;
   $lv[7][1] = 0;
   $lv[7][2] = 0;
   $lv[7][3] = 0;
   $lv[7][4] = 0;
   $lv[7][5] = 0;
   $lv[7][6] = 0;
   $lv[7][7] = 0;
   
   ### make vectors continuous
   my $sgnbt   = ( ($by>0.0) || ( EPS_ZERO($by) && ($bz>0.0) ) ) ? 1.0 : -1.0;
   my $csfactor = $a2 > $ca2 ? $sgnbt : 1.0;
   my $cffactor = $a2 < $ca2 ? $sgnbt : 1.0;
   
   for( my $k = 0; $k < 7; $k++ ){
      $rv[0][$k] *= $cffactor;
      $rv[6][$k] *= $cffactor;
      $lv[0][$k] *= $cffactor;
      $lv[6][$k] *= $cffactor;
      $rv[2][$k] *= $csfactor;
      $rv[4][$k] *= $csfactor;
      $lv[2][$k] *= $csfactor;
      $lv[4][$k] *= $csfactor;
   }
   
   return \@lambda, \@rv, \@lv;
};

# ************************************************************* READ PARAMETERS

my $targ_file;
my $result_file;
my $use_double = 0;
GetOptions(
    'input|i=s'   => \$targ_file,
    'output|o=s'  => \$result_file,
    'double|d'    => \$use_double,
) or die "Usage: $0 (--input=input_file|-i input_file) (--output=output_file|-o output_file) [--double|-d]\n";

# ********************************************************************* PROCESS

my $start_time = time();

# record: 4 byte float or 8 byte double * 11 vars
my $field_size    = $use_double ? 8 : 4;
my $field_format  = $use_double ? "d*" : "f*";
my $field_count   = 11;
my $record_size   = $field_count*$field_size;

my $targ_dir   = dirname( $targ_file );
my $targ_name  = basename( $targ_file );
$result_file = "$targ_dir/waves-$targ_name" unless $result_file;

open my $DATAFILE, '<:raw', "$targ_file" or die "Cannot open file '$targ_file' for reading";
open my $RESULTFILE, '>', "$result_file" or die "Cannot open file '$result_file' for writing";

# turn on auto flushing
$| = 1;
# prepare progress report
print "Processing: $result_file\n";
print "> ";
my $progress_message = "";

# header
my $r = 0;
$r++; printf $RESULTFILE "# $r) x";
$r++; printf $RESULTFILE "  $r) y";
printf $RESULTFILE "\n";
$r++; printf $RESULTFILE "# $r) x-1-neg-fast";
$r++; printf $RESULTFILE "  $r) x-2-neg-alfven";
$r++; printf $RESULTFILE "  $r) x-3-neg-slow";
$r++; printf $RESULTFILE "  $r) x-4-entropy";
$r++; printf $RESULTFILE "  $r) x-5-pos-slow";
$r++; printf $RESULTFILE "  $r) x-6-pos-alfven";
$r++; printf $RESULTFILE "  $r) x-7-pos-fast";
$r++; printf $RESULTFILE "  $r) x-8-dummy";
printf $RESULTFILE "\n";
$r++; printf $RESULTFILE "# $r) y-1-neg-fast";
$r++; printf $RESULTFILE "  $r) y-2-neg-alfven";
$r++; printf $RESULTFILE "  $r) y-3-neg-slow";
$r++; printf $RESULTFILE "  $r) y-4-entropy";
$r++; printf $RESULTFILE "  $r) y-5-pos-slow";
$r++; printf $RESULTFILE "  $r) y-6-pos-alfven";
$r++; printf $RESULTFILE "  $r) y-7-pos-fast";
$r++; printf $RESULTFILE "  $r) y-8-dummy";
printf $RESULTFILE "\n";

my $record;
my $record_count;
my $is_first_record = 1;
my $x_prev = 1e10;
my $i = 0;
my $j = 0;
my @targx_all;
my @targy_all;
my @targx_Uall;
my @targy_Uall;
while( my $nbytes = read( $DATAFILE, $record, $record_size ) ){
   # unpack record
   my @targ_values = unpack $field_format, $record;
   $record_count++;
   
   # read values from record
   my %targx;
   $targx{x  } = $targ_values[$field{x}];
   $targx{y  } = $targ_values[$field{y}];
   $targx{rho} = $targ_values[$field{rho}];
   $targx{u  } = $targ_values[$field{u}];
   $targx{v  } = $targ_values[$field{v}];
   $targx{w  } = $targ_values[$field{w}];
   $targx{bx } = $targ_values[$field{bx}];
   $targx{by } = $targ_values[$field{by}];
   $targx{bz } = $targ_values[$field{bz}];
   $targx{p  } = $targ_values[$field{p}];
   
   my %targy;
   $targy{x  } = $targ_values[$field{x}];
   $targy{y  } = $targ_values[$field{y}];
   $targy{rho} = $targ_values[$field{rho}];
   $targy{u  } = $targ_values[$field{v}];
   $targy{v  } = -$targ_values[$field{u}];
   $targy{w  } = $targ_values[$field{w}];
   $targy{bx } = $targ_values[$field{by}];
   $targy{by } = -$targ_values[$field{bx}];
   $targy{bz } = $targ_values[$field{bz}];
   $targy{p  } = $targ_values[$field{p}];
   
   # new line
   if( ($targx{x} != $x_prev) && (not $is_first_record) ){
      $i++;
      $j = 0;
   }
   
   #say "( $i, $j ) ( $targx{x}, $targx{y} )";
   $targx_all[$i][$j] = {
      x   => $targx{x  },
      y   => $targx{y  },
      rho => $targx{rho},
      u   => $targx{u  },
      v   => $targx{v  },
      w   => $targx{w  },
      bx  => $targx{bx },
      by  => $targx{by },
      bz  => $targx{bz },
      p   => $targx{p  },
   };
   $targy_all[$i][$j] = {
      x   => $targy{x  },
      y   => $targy{y  },
      rho => $targy{rho},
      u   => $targy{u  },
      v   => $targy{v  },
      w   => $targy{w  },
      bx  => $targy{bx },
      by  => $targy{by },
      bz  => $targy{bz },
      p   => $targy{p  },
   };
   my %targx_prev;
   $targx_prev{x  } = $i > 0 ? $targx_all[$i-1][$j]{x  } : $targx_all[$i][$j]{x  };
   $targx_prev{y  } = $i > 0 ? $targx_all[$i-1][$j]{y  } : $targx_all[$i][$j]{y  };
   $targx_prev{rho} = $i > 0 ? $targx_all[$i-1][$j]{rho} : $targx_all[$i][$j]{rho};
   $targx_prev{u  } = $i > 0 ? $targx_all[$i-1][$j]{u  } : $targx_all[$i][$j]{u  };
   $targx_prev{v  } = $i > 0 ? $targx_all[$i-1][$j]{v  } : $targx_all[$i][$j]{v  };
   $targx_prev{w  } = $i > 0 ? $targx_all[$i-1][$j]{w  } : $targx_all[$i][$j]{w  };
   $targx_prev{bx } = $i > 0 ? $targx_all[$i-1][$j]{bx } : $targx_all[$i][$j]{bx };
   $targx_prev{by } = $i > 0 ? $targx_all[$i-1][$j]{by } : $targx_all[$i][$j]{by };
   $targx_prev{bz } = $i > 0 ? $targx_all[$i-1][$j]{bz } : $targx_all[$i][$j]{bz };
   $targx_prev{p  } = $i > 0 ? $targx_all[$i-1][$j]{p  } : $targx_all[$i][$j]{p  };
   my %targy_prev;
   $targy_prev{x  } = $j > 0 ? $targy_all[$i][$j-1]{x  } : $targy_all[$i][$j]{x  };
   $targy_prev{y  } = $j > 0 ? $targy_all[$i][$j-1]{y  } : $targy_all[$i][$j]{y  };
   $targy_prev{rho} = $j > 0 ? $targy_all[$i][$j-1]{rho} : $targy_all[$i][$j]{rho};
   $targy_prev{u  } = $j > 0 ? $targy_all[$i][$j-1]{u  } : $targy_all[$i][$j]{u  };
   $targy_prev{v  } = $j > 0 ? $targy_all[$i][$j-1]{v  } : $targy_all[$i][$j]{v  };
   $targy_prev{w  } = $j > 0 ? $targy_all[$i][$j-1]{w  } : $targy_all[$i][$j]{w  };
   $targy_prev{bx } = $j > 0 ? $targy_all[$i][$j-1]{bx } : $targy_all[$i][$j]{bx };
   $targy_prev{by } = $j > 0 ? $targy_all[$i][$j-1]{by } : $targy_all[$i][$j]{by };
   $targy_prev{bz } = $j > 0 ? $targy_all[$i][$j-1]{bz } : $targy_all[$i][$j]{bz };
   $targy_prev{p  } = $j > 0 ? $targy_all[$i][$j-1]{p  } : $targy_all[$i][$j]{p  };
   
   # L_{i-1/2}
   my ( $targx_mid_lambda_ref, $targx_mid_rv_ref, $targx_mid_lv_ref ) = eigens_at_midpoint( \%targx_prev, \%targx );
   my ( $targy_mid_lambda_ref, $targy_mid_rv_ref, $targy_mid_lv_ref ) = eigens_at_midpoint( \%targy_prev, \%targy );
   # F_{i-1}
   my $targx_mid_prevF_ref = flux_F( \%targx_prev );
   my $targy_mid_prevG_ref = flux_F( \%targy_prev );
   # F_{i}
   my $targx_mid_nextF_ref = flux_F( \%targx );
   my $targy_mid_nextG_ref = flux_F( \%targy );
   # L_{i-1/2} F_{i-1}
   # L_{i-1/2} F_{i}
   my $targx_mid_prevLF_ref = m_dot_v( $targx_mid_lv_ref, $targx_mid_prevF_ref );
   my $targy_mid_prevLG_ref = m_dot_v( $targy_mid_lv_ref, $targy_mid_prevG_ref );
   my $targx_mid_nextLF_ref = m_dot_v( $targx_mid_lv_ref, $targx_mid_nextF_ref );
   my $targy_mid_nextLG_ref = m_dot_v( $targy_mid_lv_ref, $targy_mid_nextG_ref );
   # deref
   my @targx_mid_prevLF = @$targx_mid_prevLF_ref;
   my @targy_mid_prevLG = @$targy_mid_prevLG_ref;
   my @targx_mid_nextLF = @$targx_mid_nextLF_ref;
   my @targy_mid_nextLG = @$targy_mid_nextLG_ref;
   
   # L_{i-1/2} ( F_{i} - F_{i-1} )
   my @targx_mid_dLF;
   my @targy_mid_dLG;
   for( my $k = 0; $k <= 7; $k++ ){
      $targx_mid_dLF[$k] = $targx_mid_nextLF[$k] - $targx_mid_prevLF[$k];
      $targy_mid_dLG[$k] = $targy_mid_nextLG[$k] - $targy_mid_prevLG[$k];
   }
   
   $j++;
   if( $targx{x} != $x_prev && not $is_first_record ){
      printf $RESULTFILE "\n";
      
      # progress report
      print "\r> $progress_message";
      $progress_message = $targx{x};
      print "\r> $progress_message";
      $progress_message =~ s/./ /g;
   }
   $is_first_record = 0;
   $x_prev = $targx{x};
   
   printf $RESULTFILE array_to_string_f( [$targx{x},$targx{y}] );
   printf $RESULTFILE array_to_string_f( \@targx_mid_dLF );
   printf $RESULTFILE array_to_string_f( \@targy_mid_dLG );
   
   #printf $RESULTFILE array_to_string_f( [$targx{rho},$targx{u},$targx{v},$targx{w},$targx{bx},$targx{by},$targx{bz},$targx{p}] );
   #printf $RESULTFILE array_to_string_f( [$targx_prev{rho},$targx_prev{u},$targx_prev{v},$targx_prev{w},$targx_prev{bx},$targx_prev{by},$targx_prev{bz},$targx_prev{p}] );
   #printf $RESULTFILE array_to_string_f( [$targx{rho},$targx{u},$targx{v},$targx{w},$targx{bx},$targx{by},$targx{bz},$targx{p}] );
   #printf $RESULTFILE array_to_string_f( [$targy_prev{rho},$targy_prev{u},$targy_prev{v},$targy_prev{w},$targy_prev{bx},$targy_prev{by},$targy_prev{bz},$targy_prev{p}] );
   #printf $RESULTFILE array_to_string_f( [$targy{rho},$targy{u},$targy{v},$targy{w},$targy{bx},$targy{by},$targy{bz},$targy{p}] );
   #printf $RESULTFILE array_to_string_f( $targx_mid_prevF_ref );
   #printf $RESULTFILE array_to_string_f( $targx_mid_nextF_ref );
   #printf $RESULTFILE array_to_string_f( $targy_mid_prevG_ref );
   #printf $RESULTFILE array_to_string_f( $targy_mid_nextG_ref );
   
   printf $RESULTFILE "\n";
}

print "\r  $progress_message\r";

close $DATAFILE;
close $RESULTFILE;

my $end_time = time();
printf "Finished in %.2f seconds.\n", $end_time - $start_time;
