# simple awk script that picks out a row 
# from a tsv file and prints
# the values in a
#
# column_name:\tvalue
#
# format.

# use -v on invocation to set id to the value of the first column that picks
# out the row.


BEGIN {
	FS="\t";
	OFS="\t";
}

NR == 1 {
	for(i=1;i<=NF;i++) col[i]=$i;
	next;
}

id == $1 {
	for(i=1;i<=NF;i++) {
		print col[i] ":", $i;
	}
	n++
}


END {
	if(n>1) {
		print "More than one instance of", id, "while using tsv_row2yaml.awk.  That's a fail! " > "/dev/stderr";
		exit(1);
	}
}
