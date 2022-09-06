#!/usr/bin/perl
### rmg2grd.pl

$] >= 5.004 or die "Perl version must be >= 5.004 (Currently $]).\n";

use Env qw(INT_BIN INT_SCR MY_SCR);
use lib "$INT_SCR";  #### Location of Generic.pm
use Generic;

###Usage info/check
sub Usage{

`$INT_SCR/pod2man.pl  $MY_SCR/rmg2grd.pl`;
exit 1;
}
@ARGV >= 1 or Usage();
@args = @ARGV;

$in_rmg          = shift;
$proj            = shift;
$type            = shift;
$incid             = shift; # incidence angle for slant range

$proj or $proj = "LL";  # default geocoded latlong
$type or $type = "unw";  # default unwrapped phase
$incid or $incid = 23.; # default to ERS and Envisat IS2 incidence

if ($type eq "unw")
  {
    $in_file = "$in_rmg.unw";
    $in_mag = "$in_rmg.mag";
    $in_phs = "$in_rmg.phs";
    $mag_grd = "$in_mag.grd";
    $phs_grd = "$in_phs.grd";
    $los_grd = "$in_rmg.los.grd";
    @Outfiles = ("$mag_grd","$phs_grd","$los_grd");
    $phs_unit = "radians";
    $phs_var = "unwrapped_phase";
  }
elsif ($type eq "cor")
  {
    $in_file = "$in_rmg.cor";
    $in_mag = "$in_rmg.amp";
    $in_phs = "$in_rmg.corr";
    $mag_grd = "$in_mag.grd";
    $phs_grd = "$in_phs.grd";
    @Outfiles = ("$mag_grd","$phs_grd");
    $phs_unit = "=";
    $phs_var = "correlation";
  }
else
  {
    print "$type not recognized\n" and Usage();
  }

$Pi =  4 * atan2( 1, 1 );

#################
Message "Checking I/O";
#################
@Infiles  = ("$in_file.rsc");
&IOcheck(\@Infiles, \@Outfiles);

##############################################
Message "Reading resource file: $in_file.rsc";
##############################################
$length             = Use_rsc "$in_file read FILE_LENGTH";
$width              = Use_rsc "$in_file read WIDTH";
$wavelength         = Use_rsc "$in_file read WAVELENGTH";

if ($proj eq "LL" or $proj eq "UTM")
  {
    $first_E            = Use_rsc "$in_file read X_FIRST";
    $pixel_size_E       = Use_rsc "$in_file read X_STEP";
    $first_N            = Use_rsc "$in_file read Y_FIRST";
    $pixel_size_N       = Use_rsc "$in_file read Y_STEP";
    $unit_E             = Use_rsc "$in_file read X_UNIT";
    $unit_N             = Use_rsc "$in_file read Y_UNIT";

    # GMT uses opposite convention for pixel "north" and "down"
    $pixel_size_N = -1*$pixel_size_N;

    $right_E = $first_E + $width * $pixel_size_E;
    $bottom_N = $first_N - $length * $pixel_size_N; # assume first at top (N)

# first pixel at NW (top left) for geocoded
    $order = "-ZTLf";

  }
elsif ($proj eq "SR")
  {
    $incid_r = $incid*$Pi/180.; # conver to radians

    $first_E            = 0.0;
    $rng_pix_size       = Use_rsc "$in_file read RANGE_PIXEL_SIZE";
    $pixel_size_E       = $rng_pix_size/sin($incid_r);
    $bottom_N            = 0.0;
    $pixel_size_N       = Use_rsc "$in_file read AZIMUTH_PIXEL_SIZE";
    $unit_E             = "approx_meters";
    $unit_N             = "approx_meters";

    $right_E = $first_E + $width * $pixel_size_E;
    $first_N = $bottom_N + $length * $pixel_size_N; # assume first at top (N)

# first pixel at bottom left for raw
    $order = "-ZBLf";

    print "non-geocoded using $incid incidence so ground range set to $pixel_size_E\n";
  }
elsif ($proj eq "raw")
  {
    $first_E            = 0.0;
    $pixel_size_E       = 1.0;
    $bottom_N           = 0.0;
    $pixel_size_N       = 1.0;
    $unit_E             = "pixels";
    $unit_N             = "pixels";

# use postive along-track (but upward) coordinates 2006/1/19
    $right_E = $first_E + $width * $pixel_size_E;
    $first_N = $bottom_N + $length * $pixel_size_N;

# first pixel at bottom left for raw
    $order = "-ZBLf";

    print "raw uses original pixel units, positive 'down'\n";
  }
else
  {
    print "projection $proj not recognized" and Usage();
  }

$phs_fact = $wavelength / (4 * $Pi);

print "$Pi, $wavelength\n";
print "Convert phase to los factor: $phs_fact\n";

############################################
Message "splitting to $in_mag and $in_phs ";
############################################

system "$INT_BIN/rmg2mag_phs $in_file $in_mag $in_phs $width";

############################################
Message "converting to GMT grids $mag_grd and $phs_grd ";
############################################

$reg = "-R$first_E/$right_E/$bottom_N/$first_N -r"; # force pixel reg.
$incr = "-I$pixel_size_E/$pixel_size_N";

print "region is $reg, increments are $incr\n";

$mag_D =  "-D$unit_E/$unit_N/=/=/=/SAR_amplitude/uncalibrated";
$phs_D =  "-D$unit_E/$unit_N/$phs_unit/=/=/$phs_var/";

print "D: $phs_D";

system "gmt xyz2grd $in_mag $reg $incr $order -di0.0 $mag_D -G$mag_grd -V";
system "gmt xyz2grd $in_phs $reg $incr $order -di0.0 $phs_D -G$phs_grd -V";


if ($type eq "unw")
  {
    # convert to LOS meters
    system "gmt grdmath $phs_grd $mag_grd OR $phs_fact MUL = $los_grd";
    system "gmt grdedit $los_grd -D=/=/meters/=/=/range_change/converted_wavelength_$wavelength";
  }

system "rm -f $in_mag $in_phs";

exit 0;
################################################################################

=pod

=head1 USAGE

B<rmg2grd.pl>  I<rmg_file_root> I<proj> I<type> I<[incid]>

I<rmg_file_root>    = root name of rmg file to be converted

I<proj> = file map projection: LL, UTM, SR (slant range), or raw (pixel units)
I<type> = file type: unw or cor
I<[incid]> = optional incidence angle for SR projection (default 23 degrees)

=head1 FUNCTION

Converts I<rmg_file_root.type> to two GMT grids

=head1 ROUTINES CALLED

=head1 CALLED BY

=head1 FILES USED

I<rmg_file_root.type>.rsc
I<rmg_file_root.type>

=head1 FILES CREATED

I<rmg_file_root>.mag.grd
I<rmg_file_root>.phs.grd
I<rmg_file_root>.los.grd

I<or>

I<rmg_file_root>.amp.grd
I<rmg_file_root>.corr.grd

=head1 HISTORY

Eric Fielding, 2004/2/4
Eric Fielding, 2004/2/7
modified "raw" and "SR" projections to BL order  Eric Fielding, 2006/1/19
changed to allow input of incidence angle for SR projection 2006/12/3

=head1 LAST UPDATE

Eric Fielding, 2006/12/3

=cut
