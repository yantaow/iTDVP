#!/bin/bash
cmd="plot 'imps_A.out' u 2:$1 w l, 'imps_Z.out' u 2:$1 w l"

echo $cmd > hello 
echo "pause -1" >> hello
gnuplot hello; pause -1
# echo "$cmd ; pause -1"
# gnuplot -e "$cmd ; pause -1"
