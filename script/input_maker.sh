#!/bin/bash

N=(500 1000 1500 2000 4000)
L=(250 500 750 1000 2000)
t=(1000000 1000000 1000000 2000000 10000000)
ids=({1..10})

fname="script/input2.dat"

echo $N

for i in "${!N[@]}"; do
	for j in "${ids[@]}"; do
		printf "%5d %d %d 10000 %02d\n" "${N[i]}" "${L[i]}" "${t[i]}" "${ids[j]}" >> $fname
	done;
done
