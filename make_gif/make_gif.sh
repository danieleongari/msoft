#!/bin/bash

rm *.png
rm ../README.gif

gnuplot psisq_gif.gnuplot

for i in $(seq 0 9)
do
mv psisq_${i}.png psisq_0${i}.png
done

convert -delay 5 -loop 0 *.png ../README.gif
