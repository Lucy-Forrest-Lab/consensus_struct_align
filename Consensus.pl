#!usr/bin/perl

####################################
## Consensus.pl 		  ##
## Program to fuse alignments and ##
## calculate confidence values	  ##
## date: 05.05.2015 		  ##
## by Marcus Stamm 		  ##
####################################

#This is a script to which 4 different alignments are submitted
#These alignments are fused together and a confidence value for each position is provided
#The higher the value - the more confident a position is


use lib '.';
use FuseAlignments;

use warnings;
use POSIX;

my $list;
my $flag = "";
my $name_counter = 0;

my $ref = "reference";
my $obj = "object";


my $hash;
my ($aln, $id);
my $out; 

my $flag_for_inputs_used  = "false";
my $flag_for_ouput_used  = "false";

my  %confhash =();
$confhash{"9"} = 9;
$confhash{"6"} = 6;
$confhash{"4"} = 3;
$confhash{"1"} = 1;


### loop to read in and process all arguments that the user provided and checks if the correct flags are used 
for (my $i=0 ; $i < @ARGV; $i++)
{
	if ($flag eq "-help" || $flag eq "--help" )
	{
		####### add the following options: -ref [1KPL] -obj [1OTS] -m [DaliLite TMalign CE FATCAT]
		print "Usage: This script generates a consensus alignment\nRequired inputs: perl Consensus.pl -i [multiple inputs] -m [multiple identifiers] -o [output]\n";
		exit;
	}

	if (substr($ARGV[$i], 0, 1) eq "-")
	{
		$flag = $ARGV[$i]
	}
	else
	{
		if ($flag eq "")
		{
			print "Error: You did not provide any flag, please have a look at the manual or have a look at the help using perl Consensus.pl -help\n";	
			exit;	
		}
		elsif ($flag eq "-i")
		{
   			$alignment = FuseAlignments->new();
   			$alignment->id($ARGV[$i]);
			#print $alignment->id."\n";
			$alignment->alignments(get_alignments($ARGV[$i]));
			push @list, $alignment;
			$flag_for_inputs_used  = "true";
			#print $alignment->get_seq_ids(1)."\n";
			#print $alignment->get_aln(1)."\n";
		}

		elsif ($flag eq "-m")
		{
			if (@list < $name_counter)
			{
				print "Please provide first your input files and then the names of the inputs!\n";
				print "Usage:  perl Consensus.pl -i [multiple inputs] -m [multiple identifiers] -o [output]\n";
				exit;
			}
			$list[$name_counter]->name($ARGV[$i]);
			$name_counter++;
		}

		elsif ($flag eq "-obj")
		{
			$obj = $ARGV[$i];
		}
		elsif ($flag eq "-ref")
		{
			$ref = $ARGV[$i];
		}


		elsif ($flag eq "-o")
		{
			$out = $ARGV[$i];
			$flag_for_ouput_used = "true";
		}
		else
		{
			print "ERROR: The flag < ".$flag." > is not provided by this program!\n";
			exit;	
		}
	}
}

if ($flag_for_inputs_used ne "true")
{
	print "ERROR: You did not provide any inputfiles!\n";
	print "Usage: This script generates a consensus alignment\nRequired inputs:  perl Consensus.pl -i [multiple inputs] -m [multiple identifiers] -o [output]\n";
	exit;

}
if (@list != 4 ) 
{
	print "ERROR: The number of your submitted alignments is incorrect. You have to provide 4 sequence alignments in fasta format as an input for this program!\n";
	print "Required inputs:  perl Consensus.pl -i [multiple inputs] -m [multiple identifiers] -o [output]\n";
	exit;
}

#subroutine to check if the same sequences are used in the same order in the inputs of all alignments

my $oldseq1="";
my $oldseq2="";
my $newseq1="";
my $newseq2="";

for (my $i = 0; $i < @list; $i++)
{

	if ($i>0)
	{
		$val1= $i;
		$val2= $i+1;
		$newseq1=$list[$i]->get_aln("1");
		$newseq1 =~ s/-//g;
		$newseq2=$list[$i]->get_aln("2");
		$newseq2 =~ s/-//g;
		if ($newseq1 ne $oldseq1)
		{
			print "ERROR: The first sequences of your alignments $val1 and $val2 do not agree with each other. Please make sure that you submit 4 alignments in which the same sequences are aligned and they have to be in the same order in the fasta input file\n";
			exit(-1);
		}
		if ($newseq2 ne $oldseq2)
		{
			print "ERROR: The second sequences of your alignments $val1 and $val2 do not agree with each other. Please make sure that you submit 4 alignments in which the same sequences are aligned and they have to be in the same order in the fasta input file\n";
			exit(-1);
		}
	}
	$oldseq1=$list[$i]->get_aln("1");
	$oldseq1 =~ s/-//g;
	$oldseq2=$list[$i]->get_aln("2");
	$oldseq2 =~ s/-//g;
}




#### convert the sequence-based alignmen to a numerical alignment #####
######## add gaps into the template sequence #########
my @final_array = ();
for (my $i = 0; $i < @list; $i++)
{
	#print $list[$i]->id."\n";
	#$seq = $list[$i]->get_aln(1);
	$list[$i]->nr_seq_1(convert_seqaln_to_nraln($list[$i]->get_aln(1)));
	$list[$i]->nr_seq_2(convert_seqaln_to_nraln($list[$i]->get_aln(2)));

	if ($i > 0)
	{
		my @temp_array = ();
		my $pos = 0;
		my $final_pos = 0;
		my $maxpos = $list[$i]->nr_seq_1();
		while ($pos < $maxpos &&  $final_pos < @final_array)
		{
			my $recent_value =  $list[$i]->return_nrseq1_value($pos);
			my $final_value = $final_array[$final_pos];

			if($recent_value < $final_value)
			{
				push (@temp_array, $recent_value); 
				$pos++;
			}
			elsif ($recent_value == $final_value)
			{
				push (@temp_array, $recent_value); 
				$pos++;
				$final_pos++;
			}
			else
			{
				push (@temp_array, $final_value); 
				$final_pos++;
			}
				
		}
		while ($pos < $maxpos)
		{
			my $recent_value =  $list[$i]->return_nrseq1_value($pos);
			push (@temp_array, $recent_value); 
			$pos++;
		}
		while ($final_pos < @final_array)
		{
			my $final_value = $final_array[$final_pos];
			push (@temp_array, $final_array[$final_pos]); 
			$final_pos++;
		}
		@final_array = @temp_array;
			
	}
	else
	{
		for ($x =0; $x < $list[$i]->nr_seq_1(); $x++)
		{
			push  (@final_array, $list[$i]->return_nrseq1_value($x));
		}
	}
	#print "listpeers \n\n";
	#print "length of final array is: ".@final_array."\n";

	#@otherarray =  $list[$i]->nr_seq_1();

	#print "entry is ".$otherarray[0]."\n";
	#print "entry is ". $list[$i]->return_nrseq1_value(0)."\n";
	#print "\n\n";
}



######## get the final sequences #########
for ($i = 0; $i < @list; $i++)
{
	my $final_seq1 = "";
	my $final_seq2 = "";
	my $x = 0;
	my $pos=0;
	while ($x <  @final_array)
	{
		my$test="";
		$test = $list[$i]->return_nrseq1_value($pos);
		if (defined $test && $test ne '') {
		    # do something with $name
		}
		else
		{
			$test=-1;
		}
		if ($test != $final_array[$x])
		{
			$final_seq1 .= "-";
			$final_seq2 .= ".";
			$x++;
		}
		else	
		{
			$final_seq1 .= 	$list[$i]->return_seq1_value($pos);
			$final_seq2 .= 	$list[$i]->return_seq2_value($pos);
			$x++;
			$pos++;
		}
	}
	#print length($final_seq1)."\n";
	#while (length($final_seq2) < length($final_seq1))
	#{
	#	$final_seq2 .= ".";
	#}
	
	$list[$i]->set_final_alignments(1, $final_seq1);
	$list[$i]->set_final_alignments(2, $final_seq2);

	$list[$i]->set_final_alignments_numbers(1, $final_seq1);
	$list[$i]->set_final_alignments_numbers(2, $final_seq2);


	if ($i == 0)
	{
		#print $list[$i]->get_final_alignments(1)."\n";
		#print $list[$i]->get_final_alignments(2)."\n";
	}
	else
	{
		#print $list[$i]->get_final_alignments(2)."\n";		
	}
}

if ($flag_for_ouput_used eq "true")
{
	open(OUT, '>'.$out);
	print OUT "CLUSTAL W formatted alignment\n\n";
}
else
{
	print "CLUSTAL W formatted alignment\n\n";
}



# print out the numbered alignment 
#for ($x=0; $x <  @final_array; $x++)
#{
#	print $x."\t";
#	print $list[0]->return_finalnrseq1_value($x)."x";
#	print $list[0]->get_substr_of_final_alignments(1, $x,1 )."\t";
#	for ($i = 0; $i < @list; $i++)
#	{
#		print $list[$i]->return_finalnrseq2_value($x)."x";
#		print $list[$i]->get_substr_of_final_alignments(2, $x,1 )."\t";
#	}
#	print "\n";
#
#}


if (length($ref) < 15)
{
	for(my $i=length($ref); $i < 15; $i++)
	{
		$ref .= " ";	
	}
}
#print length($obj);


########### print out the alignment #########
my $confidence = "";
my $clustalw_id = 13;
my $clustalw_seq = 60;
my $total_size = length($list[0]->get_final_alignments(1));

$rows = floor($total_size/$clustalw_seq);

for ($r = 0; $r < $rows; $r++)
{
	for ($i = 0; $i < @list; $i++)
	{
		if ($i == 0)
		{
			if ($flag_for_ouput_used eq "true")
			{
				print  OUT $ref.$list[$i]->get_substr_of_final_alignments(1, $r*$clustalw_seq, $clustalw_seq)."\n";	
			}
			else
			{
				print  $ref.$list[$i]->get_substr_of_final_alignments(1, $r*$clustalw_seq, $clustalw_seq)."\n";	
			}
		}
		#$id = remove_extension($list[$i]->id);
		#if (length($id) <  $clustalw_id )
		#{
		#	for ($l = length($id); $l< $clustalw_id ; $l++)
		#	{	
		#		$id.= " ";
		#	}		
		#}
		
		my $name = "";
		if ($list[$i]->name() ne "")
		{
			$name = $list[$i]->name();
		}
		else 
		{ 
			$name = $obj."-".$i;
		}
		if (length($name) < 14)
		{
			for(my $j=length($name); $j < 14; $j++)
			{
				$name .= " ";	
			}
		}

		if ($flag_for_ouput_used eq "true")
		{
			#print OUT substr($id, length($id)-$clustalw_id , $clustalw_id )."  ".$list[$i]->get_substr_of_final_alignments(2, $r*$clustalw_seq, $clustalw_seq)."\n";
			print OUT substr($name, 0, 14)." ".$list[$i]->get_substr_of_final_alignments(2, $r*$clustalw_seq, $clustalw_seq)."\n";
		}
		else
		{
			#print substr($id, length($id)-$clustalw_id , $clustalw_id )."  ".$list[$i]->get_substr_of_final_alignments(2, $r*$clustalw_seq, $clustalw_seq)."\n";
			print substr($name, 0, 14)." ".$list[$i]->get_substr_of_final_alignments(2, $r*$clustalw_seq, $clustalw_seq)."\n";
		}
	}
	
	if ($flag_for_ouput_used eq "true")
	{
		print OUT "confidence     ";
	}
	else
	{
		print "confidence     ";
	}
	for ($confpos = $r*$clustalw_seq; $confpos < $r*$clustalw_seq+$clustalw_seq; $confpos++)
	{
		my %acids = ();
		for ($i = 0; $i < @list; $i++)
		{
			my $value = "x";
			### change here to check also neighboring residues or borders of gap or something like this. 
			### alternative: convert string to numbers with gaps being 0.5 values  (set_final_alignments)
			my $actual_acid = $list[$i]->return_finalnrseq2_value($confpos);
			my $gapcheck = $list[$i]->get_substr_of_final_alignments(2, $confpos,1 );
			if ($gapcheck ne ".")
			{
				if ($acids{$actual_acid})
				{
					$acids{$actual_acid} += 1; 
				}
				else
				{
					$acids{$actual_acid} = 1; 				
				}
			}
		}
		$max = 0;
		foreach my $key ( keys %acids )
		{
			#print $key."-".$acids{$key}."";
			if ($acids{$key} >  $max)
			{
				$max = $acids{$key};
			}
		}
		my $confvalue =  ($max/@list)*10-1;
		if ($flag_for_ouput_used eq "true")
		{
			print OUT $confhash{floor($confvalue)};
		}
		else
		{
			print $confhash{floor($confvalue)};
		}
	}		

	if ($flag_for_ouput_used eq "true")
	{
		print OUT "\n\n\n";
	}
	else
	{
		print  "\n\n\n";			
	}
	
}



########### print out the last block of the alignment (if the length of the sequence is smaller than the max. allowed length) #########

if ( $rows*$clustalw_seq < $total_size)
{
	for ($i = 0; $i < @list; $i++)
	{
		if ($i == 0)
		{
			if ($flag_for_ouput_used eq "true")
			{
				print OUT $ref.$list[$i]->get_substr_of_final_alignments(1, $rows*$clustalw_seq, $total_size-$rows*$clustalw_seq)."\n";	
			}
			else
			{
				print $ref.$list[$i]->get_substr_of_final_alignments(1, $rows*$clustalw_seq, $total_size-$rows*$clustalw_seq)."\n";	
			}
		}
		#$id = remove_extension($list[$i]->id);
		#if (length($id) <  $clustalw_id )
		#{
		#	for ($l = length($id); $l< $clustalw_id ; $l++)
		#	{	
		#		$id.= " ";
		#	}		
		#}

		my $name = "";
		if ($list[$i]->name() ne "")
		{
			$name = $list[$i]->name();
		}
		else 
		{ 
			$name = $obj."-".$i;
		}
		if (length($name) < 14)
		{
			for(my $j=length($name); $j < 14; $j++)
			{
				$name .= " ";	
			}
		}


		if ($flag_for_ouput_used eq "true")
		{
			#print OUT substr($id, length($id)-$clustalw_id , $clustalw_id )."  ".$list[$i]->get_substr_of_final_alignments(2, $rows*$clustalw_seq, $total_size-$rows*$clustalw_seq)."\n";
			print OUT substr($name, 0, 14)." ".$list[$i]->get_substr_of_final_alignments(2, $rows*$clustalw_seq, $total_size-$rows*$clustalw_seq)."\n";
		}
		else
		{
			#print  substr($id, length($id)-$clustalw_id , $clustalw_id )."  ".$list[$i]->get_substr_of_final_alignments(2, $rows*$clustalw_seq, $total_size-$rows*$clustalw_seq)."\n";
			print  substr($name, 0, 14)." ".$list[$i]->get_substr_of_final_alignments(2, $rows*$clustalw_seq, $total_size-$rows*$clustalw_seq)."\n";
		}
	}

	if ($flag_for_ouput_used eq "true")
	{
		print OUT "confidence     ";
	}
	else
	{
		print "confidence     ";
	}
	for ($confpos =  $rows*$clustalw_seq; $confpos <$total_size; $confpos++)
	{
		my %acids = ();
		for ($i = 0; $i < @list; $i++)
		{
			my $value = "x";
			my $actual_acid = $list[$i]->return_finalnrseq2_value($confpos);
			my $gapcheck = $list[$i]->get_substr_of_final_alignments(2, $confpos,1 );
			if ($gapcheck ne ".")
			{
				if ($acids{$actual_acid})
				{
					$acids{$actual_acid} += 1; 
				}
				else
				{
					$acids{$actual_acid} = 1; 				
				}
			}
		}
		$max = 0;
		foreach my $key ( keys %acids )
		{
			#print $key."-".$acids{$key}."";
			if ($acids{$key} >  $max)
			{
				$max = $acids{$key};
			}
		}
		my $confvalue =  ($max/@list)*10-1;
		if ($flag_for_ouput_used eq "true")
		{
			print OUT $confhash{floor($confvalue)};
		}
		else
		{
			print $confhash{floor($confvalue)};
		}
	}		

	if ($flag_for_ouput_used eq "true")
	{
		print OUT "\n\n\n";
	}
	else
	{
		print  "\n\n\n";			
	}
}


#subroutine to read in the alignments! 
sub get_alignments 
{
	my $id_hash;
	my $seq_hash;
	my $file = $_[0];
	my $i = 0;
	open (IN, '<'.$file) or die "Cann't open file: $file - please check if you have submitted correct filenames!";
	while (<IN>)
	{
		#print $_;
		chomp($_);
		if (substr($_, 0, 1) eq ">")
		{
			$i++;
			$id_hash->{$i} = substr($_, 1, 6);
			$seq_hash->{$i} = "";
		}
		else
		{
			$seq_hash->{$i} .= $_;	
		}

	}
	close(IN);

	# check which of the two sequences in the alignment has the most amino acids 
	# and use this sequence as a reference!  (COMMENTED OUT for structural alns 2013)
	my $l1 = $seq_hash->{1};
	$l1 =~ s/-//g;
	my $l2 = $seq_hash->{2};
	$l2 =~ s/-//g;
	#if (length($l2) > length($l1))
	#{
	#	my $tmpseq = $seq_hash->{1};
	#	$seq_hash->{1} = $seq_hash->{2};
	#	$seq_hash->{2} = $tmpseq;
#
#		my $tmpid = $id_hash->{1};
#		$id_hash->{1} = $id_hash->{2};
#		$id_hash->{2} = $tmpid;	
#	}
	#print $seq_hash->{1};
	return ($seq_hash, $id_hash); 
}


#subroutine to convert the sequence alignment to an alignment consisting of sequence numbers 
sub convert_seqaln_to_nraln
{

	my $rawsequence = $_[0];
	my @nr_array = ();
	my $seqpos = 0;
	for (my $pos = 0; $pos < length($rawsequence);$pos++)
	{

		if (substr($rawsequence, $pos, 1) ne "-")
		{
			
			if (($seqpos - floor($seqpos)) > 0)
			{
			    $seqpos = ceil($seqpos);
			}
			else
			{
				$seqpos++;
			} 
			#print "nogap ". $seqpos."\n";
			push (@nr_array, $seqpos); 
		}
		else
		{
			$seqpos+= 0.0001;
			#print  "gap ".$seqpos."\n";
			push (@nr_array, $seqpos); 
		}
	}
	#print "array length ".@nr_array."\n";
	#print "first entry ".$nr_array[0]."\n";
	return 	@nr_array;
}


#subroutine to remove extensions from filenames 
sub remove_extension {
    my $filename = shift @_;

    $filename =~ s/
                (.)             # matches any character
                \.              # the literal dot starting an extension
                [^.]+           # one or more NON-dots
                $               # end of the string
        /$1/x;

    return $filename;
}


