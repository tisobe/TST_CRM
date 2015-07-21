#!/opt/local/bin/perl -w

# read the Chandra ephemeris data points from STDIN, 
# find the SI-config proton attenuation factor 
# for the same times and print to STDOUT

# create a 100-year hash of seconds offsets from 1998

@ylen = (31622400,31536000,31536000,31536000);
$s = 0;
foreach (1970..2070) { $soy{$_} = $s - 883612800.0; $s += $ylen[$_ % 4] };

# read the basic SI config data

$scf = "siconfig";
open (SF, $scf) or die "Cannot read SI config file $scf\n";

$hdr = <SF>;    # read header line
@sc = <SF>;     # read data

print STDERR "  SI config file read completed\n";

@hdr = split " ",$hdr;
foreach $i (0..$#hdr) { 
    $tsc = $i if ($hdr[$i] =~ /3TSCPOS/);
    $hetg = $i if ($hdr[$i] =~ /4HPOSARO/);
    $letg = $i if ($hdr[$i] =~ /4LPOSARO/);
}

print STDERR "  3TSCPOS index = $tsc, 4HPOSARO index = $hetg, 4LPOSARO index = $letg\n";

foreach (@sc) {
    @af = split;
#    print STDERR "Incorrect number of fields at $af[0]\n" if (@hdr != @af);
    @t = split /:/, $af[0];
    print STDERR "Invalid timestamp: $af[0]\n" unless (@t == 5);
    $t = $soy{$t[0]} + ($t[1]-1)*86400 + $t[2]*3600 + $t[3]*60 + $t[4]; 
    print STDERR " Bad TSC position at $af[0]: $af[$tsc]\n" if (abs($af[$tsc]) > 1.1e5);
    print STDERR " Bad HETG position at $af[0]: $af[$hetg]\n" if ($af[$hetg] > 90 || $af[$hetg] < 0);
    print STDERR " Bad LETG position at $af[0]: $af[$letg]\n" if ($af[$letg] > 90 || $af[$letg] < 0);
    $af{$t} = 1.0;
    $af{$t} *= 0.0 if ($af[$tsc] < 0);
    $af{$t} *= 0.2 if ($af[$hetg] < 40);
    $af{$t} *= 0.5 if ($af[$letg] < 40);
}

print STDERR "  SI config hash build completed\n";

@at = sort {$a<=>$b} (keys %af);        # numerically sort attenuation timestamps

print STDERR "  SI config time sort completed\n";

# timestamp should be first item on each line read from STDIN

$i = -1;
while (<>) {
    @e = split;
    while ($at[$i+1] && abs($e[0]-$at[$i+1]) < abs($e[0]-$at[$i])) { $i++ }
    printf "%.1f %.1f %8.1f %.2f\n",$e[0],$at[$i],$e[0]-$at[$i],$af{$at[$i]};
}

