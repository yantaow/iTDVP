#!/bin/bash
cmd="plot 'imps_A.out' u 2:5 w l, 'mx_g0_1.5_g1_2.0.dat' u 1:2 w l"

echo $cmd > hello 
echo "pause -1" >> hello
gnuplot hello; pause -1
# echo "$cmd ; pause -1"
# gnuplot -e "$cmd ; pause -1"
