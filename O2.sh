#!/bin/bash

for k in `echo "bridge graph gtf meta rnacore scallop util"`
do
	echo $k
	cat $k/Makefile.am | sed 's/pg/O2/g' > tmp.am
	mv tmp.am $k/Makefile.am
done

cat Makefile.am | sed 's/pg/O2/g' > tmp.am
mv tmp.am Makefile.am
