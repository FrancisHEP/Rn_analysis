#!/usr/bin/perl -w

use strict;

my $cmd = 'condor_q -submitter mawenbo -totals | grep jobs';

my $i=0;
while (<>) {
  chomp;
  my ($r) = split /\s/;
  while(1) {
    my $out = qx($cmd);
    my $num = 0;
    if ($out =~ m/(\d+) jobs/) {
      $num = $1;
    }
    if ($num<100) {
      print $r, " at ", $i, "\n";
      system ("condor_submit <<\"EOF\"
universe = vanilla
executable = process_ana1to2.sh
arguments = $r
log = log/point_$i.log
output = log/point_$i.out
error = log/point_$i.err
queue
EOF
");
      $i++;
      last;
    } else {
      sleep 2;
    }
  }
}
