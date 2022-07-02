# A simple script to spit out every thin_int-th line in a file
# starting from line thin_start.

# The user has to define thin_int and thin_start outside of
# the script with the -v option.

# Currently, it is required that thin_start < thin_int

BEGIN {
	both_equal = 0;
	OFS="\t";
	if(thin_int == 0 || thin_start == 0) {
		print "One or both of thin_start and thin_int are not set or are set to 0 in thin_positions.awk.  That's a fail! " > "/dev/stderr"
		exit(1)
	}
	if(thin_start > thin_int) {
		print "thin_start cannot be > thin_int in thin_positions.awk.  Bailing out! " > "/dev/stderr"
		exit(1)
	}
	if(thin_start == thin_int) {
		both_equal = 1
	}
}


(NR == thin_start) || 
	(NR % thin_int == thin_start) ||
	(both_equal && NR % thin_int == 0) {print $0}

