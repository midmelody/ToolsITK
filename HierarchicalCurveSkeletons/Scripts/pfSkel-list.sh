#! /bin/bash

##
## This script reads the volume file names from a text file
## and calls the pfSkel.sh script for every volume file that it can parse 
## Output is written to the specified outputDir
## 
## The volumelist text file should contain the volume full name (including path)
## of the volumes to be processed, one file name on a line
##   The filenames should be of the form: <some_name>.<sx>x<sy>x<sz>.vol
##	Exaple: /home/user/cow.47x100x56.vol
##


if [ $# -lt 5 ] 
then
	echo "Usage: $0 <volumelist> <distCharges> <fieldStrength> <percHighDiverg> <outputDir>"
	exit 1
fi


awk '
BEGIN {
	nrex = 0;
	nrer = 0;	
}

/./ {
	##
	## If line starts with #, ignore it
	##

	if(substr($0, 1, 1) == "#") {
		printf "Ignoring line <" $0 ">\n";
		next;
	}

	##
	## Split input line into 3 sections, separated by "."
	##  acording to the file name convention: <some_name>.<sx>x<sy>x<sz>.vol
	##  we should get: the first part of the filename, the size and the extension
	##
	n = split($0, args, "\.");
	if(n != 3) {
		printf "** Error processing line: <%s>.\n", $0;
		nrer++;
		next;
	}
	
	##
	## Split the first section, into n sections separated by "/"
	##  we should get the path and the file name, and we just need the file name
	##
	n = split(args[1], pathElms, "[/]")
	volFile = pathElms[n] "." args[2];
	
	##
	## Split the second section into 3 sections, separated by "x"
	##  acording to the file name convention: <some_name>.<sx>x<sy>x<sz>.vol
	##  we should get the 3 sizes of the volume
	##
	n = split(args[2], size, "x");
	if(n != 3) {
		printf "** Error processing line: <%s>.\n", $0;
		nrer++;
		next;
	}
	
	printf "\n** Executing: \n\tpfSkel.sh " $0 " " size[1] " " size[2] " " size[3] " \n\t" dc " " fs " \n\t->" od "/" volFile ".skel\n";
	nrex++;
	}
	
	system("./pfSkel.sh " $0 " " size[1] " " size[2] " " size[3] " " dc " " fs " " phd " " od "/" volFile ".skel");
END {
	printf "********\nProcessed lines: %d.\nErrors: %d.\n********", nrex, nrer;
	printf "\nDone.\n";
}
' dc=$2 fs=$3 phd=$4 od=$5 $1	
	
