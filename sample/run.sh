#!/bin/bash -l

export OMP_NUM_THREADS=1

date
time \
  ./real_time \
  `#m:`           10   \
  `#dt0:`         0.1   \
  `#dt1:`         0.01 \
  `#n_sw0:`       1000  \
  `#n_sw1:`       1000 \
  `#J:`           -1.0 \
  `#g0:`          1.5  \
  `#g1:`          2.0
