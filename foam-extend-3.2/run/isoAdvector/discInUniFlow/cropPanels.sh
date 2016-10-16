#!/bin/bash

cd figs
mkdir cropped
for n in $(ls *hex*.png ); do ( echo $n; convert $n -crop 685x685+130+5 cropped/${n} ); done
for n in $(ls *poly*.png ); do ( echo $n; convert $n -crop 650x650+170+27 cropped/${n} ); done
for n in $(ls *tri*.png ); do ( echo $n; convert $n -crop 650x650+170+27 cropped/${n} ); done
cd cropped
for f in *; do fn=`echo $f | sed 's/\(.*\)\.\([^.]*\)$/\1\n\2/;s/\./_/g;s/\n/./g'`; mv $f $fn; done
cd ../../
