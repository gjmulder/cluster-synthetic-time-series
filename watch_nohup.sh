#!/bin/sh

while :
do
	date
	awk 'BEGIN {perms=81} /^.1. .targets/ {count+=1} END {printf("%d/perms, progress=%3d%%\n", count, 100*count/perms)}' nohup.out
	sleep 600
done
