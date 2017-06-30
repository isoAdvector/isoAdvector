#!/bin/bash

cd figs
mkdir cropped
for n in $(ls *.png ); do ( echo $n; convert $n -crop 560x560+150+70 cropped/${n} ); done
cd cropped
for f in *; do fn=`echo $f | sed 's/\(.*\)\.\([^.]*\)$/\1\n\2/;s/\./_/g;s/\n/./g'`; mv $f $fn; done
cd ../../
