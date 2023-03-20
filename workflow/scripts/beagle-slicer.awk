
# this little awk script is designed to take a beagle file to which has
# been prepended a scaff group file that has scaff groups in the first
# column and angsd-chromosome names in the fifth column.  That scaff group
# file has a header line which will be skipped.  

# Note that the header line of the beagle file is the first one that has more
# than 6 columns.

# you must set a variable `path` with -v that gives the directory you
# want these to go into.

# and also set an `ext` variable with for example "post" or
# "glikes", etc.  So the results will be like:
# path/scaff_group_0001-beagle-post.gz or
# path/scaff_group_0001-beagle-glikes.gz, etc.


# Ensure tab delimited output
BEGIN {OFS = "\t"; IFS="\t"; gzcomm = " gzip - "}

# first, deal with the header line
NR == 1 {
	if($1 != "id" && $5 != "angsd_chrom") {
		print "Column 1 not named id and column five not named angsd_chrom.  Bailing! " > "/dev/stderr"
		exit(1)
	}
	next
}

NF > 6 && n==0 {go=1;  header=$0; n++; next}

# record the scaff group for each chromosome
go == 0 {
	sg[$5] = $1
	next
}

# deal with printing to different files
go == 1 {
	an = split($1, a, /_/)
	if(an != 2) print "Error! ANGSD beagle marker names must have one underscore!" > "/dev/stderr"

	filep = path "/" sg[a[1]] "-beagle-" ext ".gz"
	gzp = "gzip > " filep

	#print a[1], sg[a[1]], gzp

	if(f == 0) {
		op = filep
		print header | gzp
	} else {
		if(filep != op) {
			close(gzp)
			op = filep
			print header | gzp
		}
	}
	f++

	print $0 | gzp
	next
}
