#!/bin/bash
rm -v cylinder_skew
ifort -I/usr/local/lib64 -r8 -g -c cylinder_skew.f90
ifort -L/usr/local/lib -L/usr/local/lib64 -I/usr/local/lib64 -r8 -g -o cylinder_skew cylinder_skew.o -lgeometry -lafnl


