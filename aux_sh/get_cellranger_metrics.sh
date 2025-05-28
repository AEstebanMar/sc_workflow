#! /usr/bin/env bash


input=$1
sample=$2
output=$3

if [ "$input" == "" ]; then
	echo No input specified
	exit 1
fi

if [ "$sample" == "" ]; then
	echo No sample specified
	exit 1
fi

if [ "$output" == "" ]; then
	echo No output specified
	exit 1
fi

head -n 1 $input | tr "," "\n" > tmp.header
tail -n 1 $input | sed 's/\%,/\n/g' | sed 's/\",/\n/g' | tr -d "\"" > tmp.content
paste tmp.header tmp.content -d "\t" > tmp

rm $output
while IFS= read line; do
	echo -e "$sample\t$line" >> $output
done < tmp

rm *tmp*
