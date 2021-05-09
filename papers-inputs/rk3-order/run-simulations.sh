#!/bin/bash

mkdir -p "data"
mkdir -p "plots"

./rk3-order    4 > data/rk3-order-n0004.dat
./rk3-order    8 > data/rk3-order-n0008.dat
./rk3-order   16 > data/rk3-order-n0016.dat
./rk3-order   32 > data/rk3-order-n0032.dat
./rk3-order   64 > data/rk3-order-n0064.dat
./rk3-order  128 > data/rk3-order-n0128.dat
./rk3-order  256 > data/rk3-order-n0256.dat
./rk3-order  512 > data/rk3-order-n0512.dat
./rk3-order 1024 > data/rk3-order-n1024.dat

echo "# nmax err_euler err_rk3" > "data/rk3-order-scaling.dat"
tail -qn1 data/rk3-order-n*.dat | sed "s/# //" | sort -n >> "data/rk3-order-scaling.dat" 
