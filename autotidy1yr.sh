#!/bin/bash

for i in {260..280..1}
    do 
       val=$(($i + 1))
       val2=$(($val + 1))
       auto_change $i'-'$val $val'-'$val2 'mtidy.'$i 'mtidy.'$val
    done

