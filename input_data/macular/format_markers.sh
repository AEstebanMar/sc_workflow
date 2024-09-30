#! /usr/bin/env bash


file=$1

types=`cut -f 7 $1 | sort -u | grep -v "Cell fate"`

echo -e "marker\ttype" > canon.tsv

for type in $types; do
	grep $type $file | sort -k4 -r | head -n 20 | cut -f 6,7 >> canon.tsv
done
