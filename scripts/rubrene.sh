#!/usr/bin/bash

make spectral

./spectral        \
--L      5        \
--T      0.025    \
--kappa  0.0961     \
--w      0.006     \
--lambda 0.00932031  \
--NM     10         \
--dw     0.0002       \
--wmax   0.7        \
--wmin   0.00001    

#--nev   400      
#--check-eigsys     \


