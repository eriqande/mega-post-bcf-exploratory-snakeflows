# simple awk script.  Pass it a file with sample names
# one per line and it makes the Beagle3 header


BEGIN {printf("marker\tallele1\tallele2")} 

{for(i=1;i<=3;i++) printf("\t%s", $1)} 

END {printf("\n")}

