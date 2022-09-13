# simple awk script for subsetting individuals out of
# a beagle file in a single pass.  Probably still too
# slow for really large files (in which case I would
# scatter over chromosomes, etc.), but still could be
# useful.

# to use it, you must supply a file that has individuals
# and populations, whitespace delimited in it, like this:
# T199967 POP_1
# T199969 POP_4
# T199977 POP_3
# T199978 POP_2
# T199985 POP_1
# T199990 POP_4
# T199992 POP_3

# Let's call that pops.txt.
# then if you have a gzipped beagle file you can do like this
# (cat pops.txt; zcat beagle.gz) | awk -f beagle-subset.awk

# this produces separate beagle files with names like
# POP_1.beagle, POP_2.beagle, etc.


BEGIN{
	IFS="\t"
	SUBSEP="  "  # only while testing
}

# get pops and inds info
NF < 3 {
	pop_of[$1] = $2;
	indn[$1]++ 
	ind_of[$2, ++I[$2]] = $1;
	next
}


# get positions of the columns
got_head == 0 {

	###### first check that no ind is given twice
	wrongos=0
	for(i in indn) {
		if(indn[i] > 1) {
			print "Individual", i, " appears more than once in pops file.  Abort! " > "/dev/stderr";
			wrongos++
		}
	}
	if(wrongos > 0) exit 1;
	######


	###### then get info
	for(i=4;i<=NF;i++) {
		pos[$i,++n[$i]] = i; # for example pos[T199967, 1] = 4, pos[T199967, 2] = 5, etc
	}
	got_head = 1
	######



	###### error checking.  Make sure requested indivs exist
	wrongos = 0
	for(i in pop_of) {
		if(!(i in n)) {
			print "Error! Requested indiv", i, " not found in beagle header" > "/dev/stderr";
			wrongos++
		}
	}
	if(wrongos > 0) exit(1);
	#######


	##### Set up filenames.  In the future set an output
	##### directory here, too.   
	for(i in I) {
		filename[i] = i ".beagle"
	}
	#####



	##### print the header lines in each of the files
	for(i in I) {
		printf("marker\tallele1\tallele2") > filename[i]
		for(j=1;j<=I[i];j++) {
			for(k=1;k<=3;k++) printf("\t%s", ind_of[i, j]) > filename[i]
		}
		printf("\n") > filename[i]
	}
	#####
	next
}

# Now when processing each line, print the markers
# and the alleles, then cycle over the pops,
# and for each pop, print the appropriate columns

{
	for(i in I) { # cycle over pops
		printf("%s\t%s\t%s", $1, $2, $3) > filename[i]
		for(j=1;j<=I[i];j++) {  # cycle over number of indivs in pop
			for(k=1;k<=3;k++) { # cycle over the three genotypes
				printf("\t%s", $(pos[ind_of[i, j], k])) > filename[i]
			}
		}
		printf("\n") > filename[i]
	}
}


