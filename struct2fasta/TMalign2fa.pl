#!/usr/bin/perl -w

####################################
## TMalign2fa.pl 		  ##
## Program to convert DaliLite 	  ##
## alignments into fasta format	  ##
## date: 05.05.2015 		  ##
## by Marcus Stamm 		  ##
####################################

#This is a script to convert the output of a parwise structural alignment of TM-align or FR-TM-align to a pairwise sequence alignment in fasta format. 
#The script requires as an input the structural alignment output of TM-align and the two names of the headers that should be written into the fasta file 

if (@ARGV<4)
{
	print "USAGE:
	Converts alignment outputs from TM-align to fasta files
	Requires: <TMalign_in> <id1> <id2> <fast_alignment_out> 
	Example: /home/TMalign/1KPL_1OTS.log 1KPL 1OTS /home/TMalign/1KPL_1OTS.fa 
	";
	exit(-1);
}



$TMfile = $ARGV[0];
$ID1 = $ARGV[1];
$ID2 = $ARGV[2];
$FASTA= $ARGV[3];

$first = "false";
$tmp = "false";
$second = "false";


open(TM, '<'.$TMfile);
open(FASTA, '>'.$FASTA);
while (<TM>)
{
	$line = $_;
	chomp($line);
	if ($second eq "true")
	{
		$second = "false";
		print FASTA ">".$ID2."\n".$line."\n";
	}
	if ($tmp eq "true")
	{
		$tmp = "false";
		$second = "true";
	}
	if ($first eq "true")
	{
		$first = "false";
		$tmp =  "true";
		print FASTA ">".$ID1."\n".$line."\n";
	}
	if ($line eq "(\":\" denotes the residue pairs of distance < 5.0 Angstrom)")
	{
		$first = "true";
		$tmp = "false";
		$second = "false";
	}

}
close (TM);
close (FASTA);
