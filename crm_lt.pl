#! /opt/local/bin/perl -w

# Track Chandra mission fluence from ACE and CRM data

# Robert Cameron
# September 2001

$ephdat = "/proj/rac/ops/ephem/longterm/dephem.gsme";
$acedat = "/proj/rac/ops/ACE/longterm/p3.dat";
$crmdat = "/proj/rac/ops/CRM/longterm/CRM_p.dat";
$attdat = "/proj/rac/ops/CRM/longterm/atten.dat";
$fludat = "/proj/rac/ops/CRM/longterm/fluence.dat";
$flxdat = "/proj/rac/ops/CRM/longterm/flux.dat";

@sw_factor = (0,1,2,0);
@crm_factor = (0,0,1,1);

open EPH, $ephdat or die "Cannot open $ephdat\n";
open ACE, $acedat or die "Cannot open $acedat\n";
open CRM, $crmdat or die "Cannot open $crmdat\n";
open ATT, $attdat or die "Cannot open $attdat\n";

# combine CRM and ACE fluxes according to CRM region

$delt = 300;

open FLU, ">$fludat" or die "Cannot open $fludat\n";
open FLX, ">$flxdat" or die "Cannot open $flxdat\n";

while (<EPH>) {
    ($teph,$alt,$dum,$dum,$dum,$dum,$fy,$mon,$d,$h,$m,$s) = split;
    ($tace,$dum,$dum,$p3) = split ' ',<ACE>;
    ($tcrm,$reg,$fmn,$dum,$dum,$dum) = split ' ',<CRM>;
    ($tatt,$dum,$dum,$af) = split ' ',<ATT>;
    print STDERR "Times mismatch: $teph, $tace, $tcrm, $tatt\n" 
	if ($tace != $teph or $tcrm != $teph or $tatt != $teph);
    $delt = $teph - $tprev if ($tprev);
    $tprev = $teph;
    print STDERR "Unexpected time step: $delt at time $teph\n" if ($delt != 300);
    $flx = $crm_factor[$reg]*$fmn + $sw_factor[$reg]*$p3;
    $flxatt = $flx * $af;
    $flu += ($flx * $delt);
    $fluatt += ($flxatt * $delt);
    $p3att = $p3 * $af;
    $flup3 += ($p3 * $delt);
    $flup3att += ($p3att * $delt);
    printf FLU "%12.1f%9.1f%2d%12.5e%12.5e%12.5e%12.5e %s %02d %02d %02d %02d %02d\n",
      $teph,$alt,$reg,$flup3,$flup3att,$flu,$fluatt,$fy,$mon,$d,$h,$m,$s;
    printf FLX "%12.1f%9.1f%2d%12.5e%12.5e%12.5e%12.5e %s %02d %02d %02d %02d %02d\n",
      $teph,$alt,$reg,$p3,$p3att,$flx,$flxatt,$fy,$mon,$d,$h,$m,$s;
}
