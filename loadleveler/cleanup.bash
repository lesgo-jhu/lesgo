#!/usr/local/bin/bash
find ./ -name "tree_data.*.c[1-9]*" | xargs rm
cd output
list=$(find ./ -name "vel.xskip*.c0")
for f in $list; do
	# isolate the time
	tmp=${f/.out.c0}
	t=${tmp/.\/vel.xskip[0-9]*.}
	echo "zipping files at t=$t"
	# construct list of all files at this time
 	list2=$(find ./ -name "vel.[xy]skip*.$t.out.c*")
	tar -cvf - $list2 | gzip > velskip.$t.tgz
	rm $list2
done
