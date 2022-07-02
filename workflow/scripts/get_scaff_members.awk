# simple awk script that picks out the members of a give
# scaff group.  If it picks out none, then it throws an error.

# use -v on invocation to set the value of sg to the scaffold_group id


BEGIN {
	FS="\t"
	OFS="\t"
}

# check that the first row has columns in the correct order
NR == 1 {
	if($1 != "id") {
		kill = 1;
		print "First column of scaff_group file not named id" > "/dev/stderr"
	}
	if($2 != "chrom") {
		kill = 1;
		print "Second column of scaff_group file not named chrom" > "/dev/stderr"
	}
	if($3 != "start") {
		kill = 1;
		print "Third column of scaff_group file not named start" > "/dev/stderr"
	}
	if($4 != "stop") {
		kill = 1;
		print "Fourth column of scaff_group file not named stop" > "/dev/stderr"
	}
	if(kill > 0) {
		exit(1);
	}
	next
}

# all others just print out the last three columns of the correct scaff group
$1 == sg {
	print $2, $3, $4;
	n++
}

END {
	if(n==0) {
		print "Didn't find any elements in scaff group", sg ". Bailing out!" > "/dev/stderr"
		exit(1);
	}
}