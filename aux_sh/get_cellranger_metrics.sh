#! /usr/bin/env bash


input=$1
sample=$2
output=$3
mode=$4

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

rm $output

if [ "$mode" == "count" ]; then
	head -n 1 $input | tr "," "\n" > tmp.header
	tail -n 1 $input | sed 's/\%,/\n/g' | sed 's/\",/\n/g' | tr -d "\"" > tmp.content
	paste tmp.header tmp.content -d "\t" > tmp

	while IFS= read line; do
		echo -e "$sample\t$line" >> $output
	done < tmp
	rm tmp.header tmp.content

elif [ "$mode" == "multi" ]; then
	cut -f 5- $input -d "," | grep -v "Sample ID" | grep -v "per probe barcode" | tr -d "\"" > tmp
	while IFS= read line; do
		new_line=`echo $line | sed 's/,/\t/1' | tr -d ","`
		echo -e $sample"\t""$new_line" >> $output
	done < tmp
else
	echo "Unrecognized mode value. Must be \"count\" or \"multi\", was \"$mode\""
	exit 1
fi

rm tmp
