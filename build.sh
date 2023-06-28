#!/bin/bash

gcc -O3 -march=native -mtune=native -fPIC -shared -o g2.so g2_functions.c