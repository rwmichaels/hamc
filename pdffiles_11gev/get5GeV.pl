#!/usr/bin/perl


open(infile, "<cteqPDF.dat.5");
$_=<infile>;
$_=<infile>;
$_=<infile>;

while ($_ = <infile>) {
    ($x,$qmu2,$uv,$dv,$u,$d,$s,$c,$ubar,$dbar$d2,$d3,$d4,$d5,$d6,$d7) = split(" ",$_);
    if ($ifile<3){
    print "$countp,$countm\n";
    $sump+=$countp;
    $summ+=$countm;
}
    $ifile+=1;
}
$asym=($sump-$summ)/($sump+$summ);

print "sump=$sump, summ=$summ asym=$asym\n";
