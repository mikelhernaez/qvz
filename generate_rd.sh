#!/bin/bash

idx=0

while [ $idx -lt 20 ]; do
	comp=$(echo "$idx*0.05" | bc -l)
	bin/qvz -f $comp -s codebook.txt $1 $2
	idx=$((idx+1))
done
