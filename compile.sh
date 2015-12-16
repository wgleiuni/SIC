#!/bin/bash

icpc -O3 main.cpp DefClass.cpp DefClass.h -I//home/quantum/FEAST/3.0/include/ -c
icpc -o a.out DefClass.o main.o -L /home/quantum/FEAST/3.0/lib/x64/ -lfeast_sparse -lfeast -mkl -lm
#icpc -g -shared-intel -debug -O0 main.cpp DefClass.cpp DefClass.h -mkl
