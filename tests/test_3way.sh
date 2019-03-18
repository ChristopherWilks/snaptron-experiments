python client/query_snaptron.py --datasrc gtex,rpc,encode1159 --bulk-query-file data/alk_alt_tss.hg38.snap.tsv --bulk-query-stdout | egrep -e 'chr2	29227673	29228883' | cut -f 14-17 > 3way.test

rm 3way.check
for f in gtex rpc encode1159; do python client/query_snaptron.py --datasrc $f --bulk-query-file data/alk_alt_tss.hg38.snap.tsv --bulk-query-stdout | egrep -e 'chr2	29227673	29228883' >> 3way.check ; done

cut -f 13 3way.check | perl -ne '$f=$_; $f=~s/^,//; print "$f";' | tr \\n \\, | perl -ne 'chomp; $f=$_; $f=~s/,$//; @s=split(/,/,$_); $s=0; $c=0; @v=sort { $a <=> $b } map { ($sid,$cn)=split(/:/,$_); $c++; $s+=$cn; $cn; } @s; $a=$s/$c; $m=($v[$c/2]+$v[($c/2)+1])/2; $sz=scalar(@v); print "$c\t$s\t"; printf("%.11f\t%.1f\n",$a,$m);' > 3way.stat1

diff 3way.test 3way.stat1 

#intersection
python client/query_snaptron.py --datasrc srav2,gtex,ct_h_s --bulk-query-file data/multi_intersection.snapin.tsv --bulk-query-stdout --func intersection | cut -f 1-12,14- > 3way.intersection.test

diff 3way.intersection.test tests/multi_intersection.snapin.tsv.expected
