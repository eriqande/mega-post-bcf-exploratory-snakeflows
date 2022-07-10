
# this is a simple script that takes lines that are have comma-separated
# Phred-scaled genotype liklelihoods starting in column three, like:
#
# CHROMPOS	REF	ALT	0,12,34	0,3,16 		
#
# (i.e., the CHROM and POS are smooshed into a single column.)
#
# and it turns each of those comma-containing columns into three columns
# of genotype likelihoods scaled to sum to one.

BEGIN {
	IFS="\t"
	base2int["A"] = 1
	base2int["C"] = 2
	base2int["G"] = 3
	base2int["T"] = 4
}

{
	printf("%s\t%s\t%s", $1, base2int[$2], base2int[$3]); 
	for(i=4;i<=NF;i++) {
		n=split($i,a,/,/); 
		if(n!=3) {
			print "Not 3 GL values at ", $1  > "/dev/stderr";
			exit;
		} 
		x=10^(-a[1]/10); 
		y=10^(-a[2]/10);
		z=10^(-a[3]/10); 
		sum=x+y+z; 
		printf("\t%e\t%e\t%e",x/sum, y/sum,z/sum); 
	} 
	printf("\n");
}
