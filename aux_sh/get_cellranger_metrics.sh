#! /usr/bin/env bash
cat $1 | perl -pe 's/(\d),(\d)/$1$2/g'| sed '1 s/ /_/g' | sed 's/%//g' | sed 's/"//g' | sed 's/ /\n/g' | sed 's/,/\t/g' | awk '
	{ 
	    for (i=1; i<=NF; i++)  {
	        a[NR,i] = $i
	    }
	}
	NF>p { p = NF }
	END {    
	    for(j=1; j<=p; j++) {
	        str=a[1,j]
	        for(i=2; i<=NR; i++){
	            str=str" "a[i,j];
	        }
	        print str
	    }
	}' | awk -v var="$2" 'BEGIN {FS=OFS="\t"} {print var, $0}' | sed 's/ /\t/g'
