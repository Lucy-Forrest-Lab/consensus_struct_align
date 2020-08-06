#!/usr/bin/perl -w

####################################
## matt2fa.pl 			  ##
## Program to convert MATT 	  ##
## alignments into fasta format	  ##
## date: 05.05.2015 		  ##
## by Marcus Stamm 		  ##
####################################

#This is a script to convert the output of a parwise structural alignment of MATT to a pairwise sequence alignment in fasta format. 
#The script requires as an input the structural alignment output of Matt and a specification of the filetype name. 


$inputfile = $ARGV[0];
$inputfile_ending = $ARGV[1];
$outputfile = $ARGV[2];

if (@ARGV<3)
{
print "	Converts alignment outputs from Matt to fasta files
	Requires: <matt_in> <file_ending> <outputfile>
	Example: /home/matt/1KPL_1OTS.txt .txt /home/matt/1KPL_1OTS.txt
";
exit(-1);
}

$firstheader = "true";


open(IN,'<'.$inputfile) or die "Cann't open file: $inputfile";
open(OUTPUT, ">".$outputfile);
while (<IN>)
{
	if (substr($_, 0, 1) eq ">"  &&	$firstheader ne "true")
	{
		print OUTPUT ">".substr($inputfile, length($inputfile)-length($inputfile_ending)-6, 6)."\n";
	}
	elsif (substr($_, 0, 1) eq ">"  &&	$firstheader eq "true")
	{
		print OUTPUT ">". substr($inputfile, length($inputfile)-length($inputfile_ending)-13, 6)."\n";
		$firstheader = "false";
	}
	else
	{
		print OUTPUT uc($_);	
	}
}

close (OUTPUT);
close (IN);
