#!/bin/bash
#checks for the following cases and ajusts as necessary to avoid an NaN value:
#every tissue has a shared value of 1 for all samples (happens for housekeeping genes, e.g. GADPH)
#effectively, if it were't for multiple groups being potentially present, we could just dump a pvalue of 1
#w/o needing to actually run the test

cat /dev/stdin | egrep -v -e '	((rail_id)|(header))' | perl -ne 'chomp; $f=$_; ($g,$rid,$c,$t)=split(/\t/,$f);  if($t eq "tissue") { $header=$f; next; } push(@{$a{$g}},$f); if($c==1) { $cn{$g}++; } $h{$g}->{$t}++; END { print "$header\n"; for $grp (keys %cn) { if($cn{$grp} == 9662) { for $t (keys %{$h{$grp}}) { print "$grp\t-1\t1\t$t\n"; print "$g\t-1\t2\t$t\n"; } } else { for $a (@{$a{$grp}}) { print "$a\n"; } } } }'
