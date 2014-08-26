#!/bin/bash

for i in {17971200..19353600..69120}
    do 
       val=$(($i + 69120))
       auto_change 'nIter0='$i 'nIter0='$val 'data.'$i 'data.'$val
    done

