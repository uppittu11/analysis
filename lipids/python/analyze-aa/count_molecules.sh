#!/bin/bash

file=$1

first_line=$( grep -n type $file | awk -F":|\n" 'NR==1{ print $1 }' )
last_line=$( grep -n type $file | awk -F":|\n" 'NR==2{ print $1 }' )

mhead2=$( sed -n "$first_line,$last_line p" $file | grep -c mhead )
chead=$( sed -n "$first_line,$last_line p" $file | grep -c chead )
head=$( sed -n "$first_line,$last_line p" $file | grep -c head )
let "head=$head-$mhead2-$chead"
water=$( sed -n "$first_line,$last_line p" $file | grep -c water )

printf "$mhead2\tCeramides\n$chead\tCholesterols\n$head\tFree Fatty Acids\n$water\tWaters\n"
