#!/bin/bash

for i in {0..260..20}
    do 
       val=$(($i + 20))
       val2=$(($val + 20))
       auto_change $i'-'$val $val'-'$val2 'mtidy.'$i 'mtidy.'$val
    done

