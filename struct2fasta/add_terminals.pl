#!/usr/bin/perl -w

####################################
## add_terminals.pl 		  ##
## Program to convert DaliLite 	  ##
## alignments into fasta format	  ##
## date: 05.05.2015 		  ##
## by Marcus Stamm 		  ##
####################################

#Some structural alignment methods miss in their output terminal segments that are aligned against gaps. 
#This is a script to add those missing residues to back to the pairwise sequence alignment output.  
#The script requires as an input a pairwise sequence alignment in fasta format and the two original fasta files that correspond to that alignment



if (@ARGV<4)
{
	print "USAGE:
	Adds missing terminals to an alignment (e.g., to an alignment output of FATCAT) 
	Requires: <alignment_in> <fasta1_in> <fasta2_in> <outputfile>
	Example: /home/alignments/1KPL_1OTS.txt /home/fastas/1KPL.fa /home/fastas/1OTS.fa  /home/alignments_with_terminals/1KPL_1OTS.txt
	";
	exit(-1);
}

$hhalignfile = $ARGV[0];
$first_fasta = $ARGV[1];
$second_fasta =  $ARGV[2];
$outputfile=  $ARGV[3];

open(HHalign, '<'.$hhalignfile);
$id1 = <HHalign>;
chomp($id1);
$orig_hhseq1 = <HHalign>;
chomp($orig_hhseq1);
$id2 = <HHalign>;
chomp($id2);
$orig_hhseq2= <HHalign>;
chomp($orig_hhseq2);

$hhseq1 = "";
$hhseq2 = "";

for ($i = 0; $i < length($orig_hhseq1); $i++)
{
	if (substr($orig_hhseq1, $i, 1) ne "-")
	{
		$hhseq1 .= substr($orig_hhseq1, $i, 1);	
	}
}

for ($i = 0; $i < length($orig_hhseq2); $i++)
{
	if (substr($orig_hhseq2, $i, 1) ne "-")
	{
		$hhseq2 .= substr($orig_hhseq2, $i, 1);	
	}
}



$fastaseq1 = "";
&read_fasta ($first_fasta, $fastaseq1);
$fastaseq2 = "";
&read_fasta ($second_fasta, $fastaseq2);



print $fastaseq1."\n";
print $hhseq1."\n";
$start_gaps = index($fastaseq1, $hhseq1);
$end_gaps = length($fastaseq1) - length($hhseq1) - $start_gaps;
for ($i = 0; $i < $start_gaps;  $i++)
{
	$tempseq1 .= substr($fastaseq1, $i, 1);
	$tempseq2 .= "-";
}
$tempseq1 .= $orig_hhseq1;
$tempseq2 .= $orig_hhseq2;
for ($i = 0; $i < $end_gaps;  $i++)
{	
	$tempseq1 .=  substr($fastaseq1, length($fastaseq1)-$end_gaps+$i, 1);
	$tempseq2 .= "-";
}



$start_gaps = index($fastaseq2, $hhseq2);
$end_gaps = length($fastaseq2) - length($hhseq2) - $start_gaps;
for ($i = 0; $i < $start_gaps;  $i++)
{
	$finalseq1 .= "-";
	$finalseq2 .= substr($fastaseq2, $i, 1);
}
$finalseq2 .= $tempseq2;
$finalseq1 .= $tempseq1;
for ($i = 0; $i < $end_gaps;  $i++)
{	
	$finalseq1 .= "-";
	$finalseq2 .=  substr($fastaseq2, length($fastaseq2)-$end_gaps+$i, 1);
}




open(OUT, '>'.$outputfile);
print OUT $id1."\n".$finalseq1."\n".$id2."\n".$finalseq2."\n";
close (OUT);


sub read_fasta
{
	$input_fasta = $_[0];
	$fastaseq = $_[1];
	
	open (FASTA, '<'.$input_fasta);
	while(<FASTA>)
	{
		if (substr($_, 0,1) ne ">")
		{
			chomp($_);
			$fastaseq .= $_
		}
	}
	$_[1] = $fastaseq;
}




