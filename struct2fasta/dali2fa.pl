#!/usr/bin/perl -w

####################################
## dali2fa.pl 			  ##
## Program to convert DaliLite 	  ##
## alignments into fasta format	  ##
## date: 05.05.2015 		  ##
## by Marcus Stamm 		  ##
####################################

#This is a script to convert the output of a parwise structural alignment of DaliLite to a pairwise sequence alignment in fasta format. 
#The script requires as an input the structural alignment output of DaliLite and two fasta files with the sequences used in the alignment. 


if (@ARGV<4)
{
	print "USAGE:
	Converts alignment outputs from dalilite to fasta files
	Requires: <dali_in> <fasta1> <fasta2> <fast_alignment_out> 
	Example: /home/dalilite/1KPL_1OTS.log /home/fastas/1KPL.fa /home/fasts/1OTS.fa /home/dalilite/1KPL_1OTS.fa 
	";
	exit(-1);
}


$daliresult = $ARGV[0];
$fasta1 = $ARGV[1];
$fasta2 = $ARGV[2];
$fastaoutput = $ARGV[3];



$fasta1seq = "";
$fasta2seq = "";

$final1seq = "";
$final2seq = "";

open(FASTA1, '<'.$fasta1);
while (<FASTA1>)
{
	$line = $_;
	chomp($line);
	if (substr($line, 0 ,1) ne ">")
	{
		$fasta1seq = $fasta1seq.$line;
	}
	else
	{
		$fasta1header = $line;	
	}
}
close (FASTA1);

open(FASTA2, '<'.$fasta2);
while (<FASTA2>)
{
	$line = $_;
	chomp($line);
	if (substr($line, 0 ,1) ne ">")
	{
		$fasta2seq = $fasta2seq.$line;
	}
	else
	{
		$fasta2header = $line;	
	}
}
close (FASTA2);

$structures = "false";
$best = "nomodel";
$start = "false";

$start1 = 0;
$start2 = 0;
$end1 = 0;
$end2 = 0;
$oldend1 = -1000;
$oldend2 = -1000;

open(DALI, '<'.$daliresult);
while (<DALI>)
{

	@array = split(' ',$_);
	if (@array > 2 && $array[0] eq  "#" &&  $array[1] eq  "Structural" &&  $array[2] eq "equivalences")
	{
		$structures = "true";
	}

	if (@array > 1 && $array[0] eq "1:" && $structures eq "true")
	{
		$best = "model1";
	}
	else
	{
		if ($best eq "model1")
		{
			for ($x = $oldend1 + 1; $x <= length($fasta1seq); $x++)
			{
				$final1seq = $final1seq.substr($fasta1seq, $x-1, 1);
				$final2seq = $final2seq."-";
			}	  

			for ($x = $oldend2 + 1; $x <=  length($fasta2seq); $x++)
			{
				$final1seq = $final1seq."-";
				$final2seq = $final2seq.substr($fasta2seq, $x-1, 1);

			} 
			$best = "false";
		}
	}


	if ($structures eq "true" && $best eq "model1")
	{
		$start1 = substr($_, 21, 4);
		$end1 = substr($_, 27, 4);
		$diff1 = $end1 - $start1 +1;
		$start2 = substr($_, 36, 4);
		$end2 = substr($_, 42, 4);
		$diff2 = $end2 - $start2 + 1;
		
		if ($start eq "false")
		{
			for ($i = 1; $i < $start1; $i++)
			{
				$final1seq = $final1seq.substr($fasta1seq, $i-1, 1);	
				$final2seq = $final2seq."-";	
			}
			for ($i = 1; $i < $start2; $i++)
			{
				$final1seq = $final1seq."-";
				$final2seq = $final2seq.substr($fasta2seq, $i-1, 1);	
			}
			$final1seq = $final1seq.substr($fasta1seq, $start1-1, $diff1);
			$final2seq = $final2seq.substr($fasta2seq, $start2-1, $diff2);	

			$oldend1 = $end1;
			$oldend2 = $end2;

			$start = "true";
		}
		else
		{
			for ($x = $oldend1 + 1; $x < $start1; $x++)
			{
				$final1seq = $final1seq.substr($fasta1seq, $x-1, 1);
				$final2seq = $final2seq."-";
			}	  

			for ($x = $oldend2 + 1; $x < $start2; $x++)
			{
				$final1seq = $final1seq."-";
				$final2seq = $final2seq.substr($fasta2seq, $x-1, 1);

			}  
			$final1seq = $final1seq.substr($fasta1seq, $start1-1, $diff1);
			$final2seq = $final2seq.substr($fasta2seq, $start2-1, $diff2);	

			$oldend1 = $end1;
			$oldend2 = $end2;
		}

	}


	#$final1seq = $final1seq."-";
	#$final2seq = $final2seq.substr($fasta2seq, $i-1, 1);
	#$final1seq = $final1seq.substr($fasta1seq, $start1-1, $diff1);
	#$final2seq = $final2seq.substr($fasta2seq, $start2-1, $diff2);	
}
close (DALI);
open (OUT, '>'.$fastaoutput);
print OUT $fasta1header."\n".$final1seq."\n".$fasta2header."\n".$final2seq."\n";
close (OUT);
