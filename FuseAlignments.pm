#!usr/bin/perl
####################################
## FuseAlignments.pm is part of	  ##
## Consensus.pl 		  ##
## Program to fuse alignments and ##
## calculate confidence values	  ##
## date: 05.05.2015 		  ##
## by Marcus Stamm 		  ##
####################################
 
    package FuseAlignments;
    use strict;

    ##################################################
    ## the object constructor (simplistic version)  ##
    ##################################################
    sub new {
        my $self  = {};
	my $array = [];
        $self->{FILEID}   = undef;
        $self->{NAME}   = "";
        $self->{ALN}    = {};
        $self->{SEQIDS}    = {};
        $self->{SEQ1}    = undef;
        $self->{SEQ2}    = undef;
        $self->{FINALSEQ1}    = undef;
        $self->{FINALSEQ2}    = undef;
        $self->{NRSEQ1}  = [];
        $self->{NRSEQ2}  = [];
        $self->{FINALNRSEQ1}  = [];
        $self->{FINALNRSEQ2}  = [];
        $self->{PEERS}  = [];
        bless($self);           # but see below
        return $self;
    }

    ##############################################
    ## methods to access per-object data        ##
    ##                                          ##
    ## With args, they set the value.  Without  ##
    ## any, they only retrieve it/them.         ##
    ##############################################

    sub id 
    {
        my $self = shift;
        if (@_)
	{
		 $self->{FILEID} = $_[0];
	}
        return $self->{FILEID};
    }

    sub name
    {
        my $self = shift;
        if (@_)
	{
		 $self->{NAME} = $_[0];
	}
        return $self->{NAME};
    }


    sub alignments 
    {
        my $self = shift;

        if (@_) 
	{ 
		$self->{ALN} = $_[0];
		$self->{SEQ1} = $_[0]->{1};
		$self->{SEQ2} = $_[0]->{2};
		$self->{SEQIDS} = $_[1];
		$self->{SEQID1} = $_[1]->{1};
		$self->{SEQID2} = $_[1]->{2};
	}
	
        return $self->{ALN};
    }


	sub set_final_alignments_numbers
	{
        	my $self = shift;
		my $finalseq = $_[1];
		my @aln_nr_array;
		my $gap_open = "false";
		my $acid_counter = 0;
		for (my $i=0; $i<length($finalseq); $i++)
		{

#			if (substr($finalseq, $i, 1) eq "-" || substr($finalseq, $i, 1) eq ".")
#			{
#				if (substr($finalseq, $i, 1) eq "-" )
#				{
#					$acid_counter = $acid_counter + 0.0001;
#				}
#				push (@aln_nr_array, $acid_counter); 
#				$gap_open = "true"
#	
#			}
#			else
#			{	
#				if ($gap_open eq "true")
#				{ 
#					$acid_counter = int($acid_counter + 1);
#					push (@aln_nr_array, $acid_counter); 
#					$gap_open = "false";
#				}
#				else
#				{
#					$acid_counter++;
#					push (@aln_nr_array, int($acid_counter)); 
#				}
#			}


			if (substr($finalseq, $i, 1) eq "-" || substr($finalseq, $i, 1) eq ".")
			{
				if ($gap_open eq "true")
				{
					push (@aln_nr_array, $acid_counter); 
				}
				else
				{
					$acid_counter = $acid_counter + 0.5;
					push (@aln_nr_array, $acid_counter); 
				}
				$gap_open = "true";
			}
			else
			{	
				if ($gap_open eq "true")
				{
					$acid_counter = $acid_counter + 0.5;
					push (@aln_nr_array, $acid_counter); 
					$gap_open = "false";
				}
				else
				{
					$acid_counter++;
					push (@aln_nr_array, $acid_counter); 
				}
			}
		}
			
        	if (@_[0] eq "1") 
		{ 
			@{$self->{FINALNRSEQ1}} = @aln_nr_array;
			return @{ $self->{FINALNRSEQ1} };
		}
		elsif (@_[0] eq "2")
		{
			@{$self->{FINALNRSEQ2}} =@aln_nr_array;
			return @{ $self->{FINALNRSEQ1} };
		} 
	}


    sub return_finalnrseq1_value {
        my $self = shift;
        if (@_) { return @{ $self->{FINALNRSEQ1} }[ @_[0]];}
    }

    sub return_finalnrseq2_value {
        my $self = shift;
        if (@_) { return @{ $self->{FINALNRSEQ2} }[ @_[0]];}
    }



	sub set_final_alignments
	{
        	my $self = shift;
        	if (@_[0] eq "1") 
		{ 
			$self->{FINALSEQ1} = $_[1];
		}
		elsif (@_[0] eq "2")
		{
			$self->{FINALSEQ2} = $_[1];
		} 
	}



	sub get_final_alignments
	{
        	my $self = shift;
        	if (@_[0] eq "1") 
		{ 
			return $self->{FINALSEQ1};
		}
		elsif (@_[0] eq "2")
		{
			return $self->{FINALSEQ2};
		} 
	}

	sub get_substr_of_final_alignments
	{
        	my $self = shift;
        	if (@_[0] eq "1") 
		{ 
			return substr($self->{FINALSEQ1}, @_[1], @_[2]);
		}
		elsif (@_[0] eq "2")
		{
			return substr($self->{FINALSEQ2}, @_[1], @_[2]);
		} 
	}



	sub get_aln
  	{
        	my $self = shift;
        	if (@_[0] eq "1") 
		{ 
			return $self->{SEQ1};
		}
		elsif (@_[0] eq "2")
		{
			return $self->{SEQ2};
		} 
	}




    sub return_seq1_value {
        my $self = shift;
        if (@_) { return substr($self->{SEQ1},  @_[0], 1);}
    }

    sub return_seq2_value {
        my $self = shift;
   	if (@_) { return substr($self->{SEQ2},  @_[0], 1);}
    }



	sub get_seq_ids
  	{
        	my $self = shift;
        	if (@_[0] eq "1") 
		{ 
			return $self->{SEQID1};
		}
		elsif (@_[0] eq "2")
		{
			return $self->{SEQID2};
		} 
	}
	

    sub nr_seq_1 {
        my $self = shift;
        if (@_) { @{ $self->{NRSEQ1} } = @_ }
        return @{ $self->{NRSEQ1} };
    }


    sub return_nrseq1_value {
        my $self = shift;
        if (@_) { return @{ $self->{NRSEQ1} }[ @_[0]];}
    }

    sub return_nrseq2_value {
        my $self = shift;
        if (@_) { return @{ $self->{NRSEQ2} }[ @_[0]];}
    }

    sub nr_seq_2 {
        my $self = shift;
        if (@_) { @{ $self->{NRSEQ2} } = @_ }
        return @{ $self->{NRSEQ2} };
    }


    1;  # so the require or use succeeds

