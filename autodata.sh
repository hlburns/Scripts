#!/bin/bash

for i in {0..16588800..691200}
    do 
       val=$(($i + 691200))
       auto_change 'nIter0='$i 'nIter0='$val 'data.'$i 'data.'$val
    done

