#!/bin/bash

mkdir -p "data"
mkdir -p "plots"

cd "data"

./euler_1d inputs/sh1d-128-u0.5.ini
./euler_1d inputs/sh1d-128-u1.0.ini
./euler_1d inputs/sh1d-128-u1.5.ini
./euler_1d inputs/sh1d-128-u2.0.ini
