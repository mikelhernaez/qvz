#!/bin/bash

idx=0
rm -f trash.txt
while [ $idx -lt 20 ]; do
	comp=$(echo "$idx*0.05" | bc -l)
	bin/qvz -f $comp -s codebook.txt $1 $2 >> trash.txt
	idx=$((idx+1))
done
awk '{print $2 $4 $6}' trash.txt > $3
rm -f trash.txt
