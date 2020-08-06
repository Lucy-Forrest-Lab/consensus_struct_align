#!/usr/bin/perl -w

####################################
## factat2fa.pl 		  ##
## Program to convert FATCAT 	  ##
## alignments into fasta format	  ##
## date: 05.05.2015 		  ##
## by Marcus Stamm 		  ##
####################################

#This is a script to convert the output of a parwise structural alignment of FATCAT to a pairwise sequence alignment in fasta format. 
#The script requires as an input the structural alignment output of FATCAT.

#Please note: 
#FATCAT is missing terminal segments that are aligned against gaps in the output
#The script add_terminals has to be run after this script to ensure that the complete sequences are aligned! 



if (@ARGV<2)
{
	print "USAGE:
	Converts alignments outputs from FASTA to fasta files
	Please not that the script add_terminals.pl might be required to add missing terminals to this alignment. 
	Requires: <FATCAT_in> <FATCAT_fasta_out> 
	Example: FATCAT.log FATCAT.fa 
	";
	exit(-1);
}

#$ce = "/home/marcus/test/fatcat.log";
#$fastaoutput = "/home/marcus/test/fatcat.fa";

$ce = $ARGV[0];
$fastaoutput = $ARGV[1];


$getsequence1 = "false";
$getsequence2 = "false";

$fastaheader1 = ">";
$fastaheader2 =  ">";
$fastaseq1 = "";
$fastaseq2 = "";

open(CE, '<'.$ce);
while(<CE>)
{
	@array = split(' ',$_);
	if (@array > 1 && $array[0] eq "Chain" && $array[1] eq "1:")
	{
		if ($getsequence1 eq "false")
		{
			$fastaheader1 = $fastaheader1.$array[2];
		}
		else
		{
			$fastaseq1 = $fastaseq1.$array[3];
		}
		$getsequence1 = "true";
		
	}
	if (@array > 1 && $array[0] eq "Chain" && $array[1] eq "2:")
	{
		if ($getsequence2 eq "false")
		{
			$fastaheader2 = $fastaheader2.$array[2];
		}
		else
		{
			$fastaseq2 = $fastaseq2.$array[3];
		}
		$getsequence2 = "true";
	}
}
close (CE);
open (OUT, '>'.$fastaoutput);
print OUT $fastaheader1."\n".$fastaseq1."\n".$fastaheader2."\n".$fastaseq2."\n";
close (OUT);
