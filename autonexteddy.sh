#!/bin/bash

for i in {2..20..1}
    do 
       fname=$(($i + 1))
       val1=$(($i + 260-1))
       val0=$(($i+260-2))
       val2=$(($i + 260))
       val3=$(($val0*69120))
       val4=$(($val3+69120))
       auto_change '/mtidy.'$val0 '/mtidy.'$val1 'nexteddy'$i 'next_eddy'$fname
       auto_change $val0'-'$val1 $val1'-'$val2 'next_eddy'$fname 'next__eddy'$fname
       auto_change '*'$val3 '*'$val4 'next__eddy'$fname 'next___eddy'$fname
       auto_change '.'$val3 '.'$val4 'next___eddy'$fname 'nexteddy'$fname
    done
wait
rm -f next_eddy* next__eddy* next___eddy*