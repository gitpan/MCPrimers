package Bio::MCPrimers;

use strict;
use warnings;

my $VERSION = '1.00';

# Bio::MCPrimers.pm - generates molecular cloning PCR primers for pET-32a

##########################################################################
# This software comes with no guarantee of usefulness.
# Use at your own risk. Check any solutions you obtain.
# Stephen G. Lenk assumes no responsibility for the use of this software.
##########################################################################

#######################################################
###### See bottom for pod and other information    ####
#######################################################


BEGIN {
    
    if ($^O =~ /^MSW/) { # Microsoft

        # check PRIMER3DIR
        if (not (defined $ENV{PRIMER3DIR})) {
            $ENV{PRIMER3DIR} = '.'
        } 

        # Is Primer3 executable available?
        if (not (-x "$ENV{PRIMER3DIR}\\primer3.exe")) {

print STDERR "\nprimer3.exe not available at $ENV{PRIMER3DIR}\\primer3.exe\n\n";
print STDERR "Check your environment variable PRIMER3DIR\n\n";

            exit 1;
        }

    } else { # non-Microsoft

        # IPC::Open3 used for non-MSW
        require IPC::Open3;

        # check PRIMER3DIR      
        if (not (defined $ENV{PRIMER3DIR})) {
            $ENV{PRIMER3DIR} = '.'
        } 

        # Is Primer3 executable available?
        if (not (-x "$ENV{PRIMER3DIR}/primer3_core")) {
                 
print STDERR "\nprimer3_core not available at $ENV{PRIMER3DIR}/primer3_core\n\n";
print STDERR "Check your environment variable PRIMER3DIR\n\n";

            exit 1;
        }
    }
}


sub find_mc_primers;
sub _solver;
sub _get_re_patterns;
sub _get_re_matches;
sub _start_codon;
sub _stop_codon;
sub _primers_ok;
sub _primer_pairs;
sub _create_primers;
sub _create_left_primers;
sub _create_right_primers;
sub _handle_reply;
sub _sanity_check_gene;
sub _primer_pairs;
sub _number_dots;
sub _generate_re_patterns;
sub _too_many_substitutions_in_a_row;
sub _check_3prime_gc;

my $CODON_SIZE = 3;   # likely to stay a constant

my %left_bad;         # for tracking bad primers during calculation
my %right_bad;        # and check right as well

my $fh_write;         # for non-MSW pipes from IPC::Open3
my $fh_read;          # read and write handles used



#### front end solver sets up for OS needs ####
sub find_mc_primers {

    my ($gene,      # ATGC string
        $flag_hr,   # anonymous hash reference to flags
        @re         # array of restriction enzyme strings
       ) = @_;
      
    my $search_shift = 0;       # permit the search for left primer
                                # to shift to the left
    my $all          = 0;       # return all solutions found
    my $clamp        = 'both';  # GC clamp at both ends  
    
    my $pid;    # process id needed to harvest child
   
    # MSW uses intermediate files, otherwise pipes
    if ( not $^O =~ /^MSW/) {
    
        # pipe on non-MSW
        $pid = IPC::Open3::open3($fh_write, 
                                 $fh_read, 
                                 $fh_read, 
               "$ENV{PRIMER3DIR}/primer3_core -format_output");
    }

    # digest flags into local variables
    if (defined $flag_hr->{search_shift}) { 
        $search_shift = $flag_hr->{search_shift} 
    } 
    if (defined $flag_hr->{all}) {
        $all = 1
    }
    if (defined $flag_hr->{clamp}) {

        if ($flag_hr->{clamp} eq 'both' or 
            $flag_hr->{clamp} eq '3prime') {

            $clamp = $flag_hr->{clamp}
        }
    }
    
    # invoke solver
    my $answer = _solver($gene, $search_shift, $all, $clamp, @re);

    # clean up intermediates for MSW
    if ($^O =~ /^MSW/) { 
    
        # intermediate files in MSW
        if (-e 'p3in.txt') {        
            system("del p3in.txt")
        }
        if (-e 'p3out.txt') {        
            system("del p3out.txt")
        }
        
    } else {
    
        # clean up pipe and child process
        close $fh_write;  
        close $fh_read;
        waitpid $pid, 0;
    }

    return $answer;
}



### solver function ####
sub _solver {

    my ($gene,         # ATGC string
        $search_shift, # extra search into start of gene for left primer
        $all,          # provide all solutions
        $clamp,        # both or 3prime
        @re            # array of RE strings
       ) = @_;

    my $MAX_CHANGES  = 3;       # changes to gene for primer
    my $solution_ar  = [];      # solution array reference

    my @re_matches;   # array of anonymous arrays of RE matches
    for (my $i = 0; $i < @re; $i++) {  $re_matches[$i] = []  }

    my $gene_start = _start_codon($gene); 
    my $gene_stop  = _stop_codon($gene, $gene_start); 

    if (not defined $gene_stop) { return undef }

    # patterns and matches for RE
    my @patterns = _get_re_patterns($re[0], $MAX_CHANGES);
    my @matches  = _get_re_matches($gene, @patterns);

    # generate RE matches
    for (my $i=0; $i < @re; $i += 1) {

        @patterns   = _get_re_patterns($re[$i], $MAX_CHANGES); 
        my @matches = _get_re_matches($gene, @patterns);

        foreach my $match (@matches) { 
            push @{$re_matches[$i]}, $match;
        }
    }
    
    # keep pulling out a potential left primer until exhausted or solved
    while (@re > 1) {

        my $re_l         = shift @re;
        my $re_pos_l_ref = shift @re_matches;

        # position match list for head enzyme
        foreach my $re_pos_l (@{$re_pos_l_ref}) {

            # control where left primer is placed
            if ($re_pos_l > (10 + $gene_start + $search_shift)) { next }
         
            # left primers
            my $modified_gene = $gene;
            substr($modified_gene, $re_pos_l, length $re_l, $re_l);
            my @left = _create_left_primers($modified_gene, 
                                            $re_pos_l, 
                                            $re_l,
                                            $clamp);
    
            # loop across rest of sites
            my $num_sites_left = @re;
            my $site = 0;
            
            # rest of enzymes left after head has been unshifted
            while ($site < $num_sites_left) {
        
                my $re_r = $re[$site];
                my $re_pos_r_ref = $re_matches[$site]; 

                # rest of restriction sites going right in vector
                foreach my $re_pos_r (@{$re_pos_r_ref}) {
    
                    if ($re_pos_r < $gene_stop) {  next  }
    
                    # right primers
                    my $more_modified_gene = $modified_gene;
                    substr($more_modified_gene, $re_pos_r, length $re_r, $re_r);                       
                    my @right = _create_right_primers($more_modified_gene, 
                                                      $re_pos_r, 
                                                      $re_r,
                                                      $clamp);
    
                    # generate primer pairs and have them checked
                    my $reply = _primer_pairs(\@left, 
                                              \@right, 
                                              $more_modified_gene, 
                                              $gene_start,
                                              $gene_stop,
                                              $search_shift);

                    # solution obtained if $reply is defined
                    if (defined $reply) { 

                        # see if solution is valid
                        # and check if one or all solutions are needed
                        if (_handle_reply($reply,
                                          $more_modified_gene,
                                          $re_l,
                                          $re_r,
                                          $solution_ar) == 1 and $all == 0) {
                        
                            # return first valid solution found
                            return $solution_ar;  
                        
                        }
                    } 
                }
                
                # next RE site for possible insertion
                $site = $site + 1;  
            }
        } 
    }

    # all solutions
    return $solution_ar;
}



#### handle reply from solution checker ####
sub _handle_reply {

    my ($reply,                 # solution checker reply
        $more_modified_gene,    # ATGC modified for left and right RE
        $re_l,                  # left RE
        $re_r,                  # right RE
        $solution_ar            # array reference to solution
       ) = @_;

    # check that RE sequences have only one match
    my $left_cnt  = 0;
    my $right_cnt = 0;
    while ($more_modified_gene =~ /$re_l/g) { 
        $left_cnt  += 1 
    }
    while ($more_modified_gene =~ /$re_r/g) { 
        $right_cnt += 1 
    }

    # one match for each?
    if ($left_cnt == 1 and $right_cnt == 1) { 

        # add anonymous hash to anonymous array
        push @{$solution_ar}, { primer3  => $reply,
                                left_re  => $re_l,
                                right_re => $re_r };
        return 1; # solution added
    }
    
    return 0; # solution not added
}



#### check primers with Primer3 ####
sub _primers_ok {

    my ($left,          # left primer string
        $right,         # right primer string
        $gene,          # ATGC string modified to match primers
        $gene_start,    # location of start codon
        $gene_stop,     # location of stop codon
        $search_shift   # move search zone into gene from left
       ) = @_;
 
    my $ok = 1;
    
    if (defined $right_bad{$right}) { 
        $ok = 0 
    }
    if (defined $left_bad{$left})   { 
        $ok = 0 
    }
    
    # outta here if one of the primers has been identified as bad
    if ($ok == 0) { 
        return undef 
    }
  
    # create Boulder file text for Primer3
    my $range = length $gene; 
    my $excluded_region_start = $gene_start + 30 + $search_shift;
    my $excluded_length =  $gene_stop - $excluded_region_start + 3;
    my @boulder = 
       ( "SEQUENCE=$gene",
         "PRIMER_PRODUCT_MAX_SIZE=$range",
         "PRIMER_PRODUCT_SIZE_RANGE=100-$range",
         "PRIMER_LEFT_INPUT=$left",
         "PRIMER_RIGHT_INPUT=$right",
         "EXCLUDED_REGION=$excluded_region_start,$excluded_length",
         "PRIMER_EXPLAIN_FLAG=1",
         "=" );

   # write Boulder file
   if ($^O =~ /^MSW/) {
    
        # MS Windows uses intermediate files
        open  $fh_write, ">p3in.txt";
        foreach (@boulder) { 
            print $fh_write "$_\n"; 
        }
    
 system(".\\$ENV{PRIMER3DIR}\\primer3.exe -format_output <p3in.txt >p3out.txt");
    
        open $fh_read, "<p3out.txt";

    } else {

       foreach (@boulder) {       
           print $fh_write "$_\n"; 
       } 
    }     
    
    my $p3_reply;
    
    # go through Primer3 output
    PRIMER3_READ: while (<$fh_read>) { 
    
        my $line = $_;
        
        $p3_reply .= $line;
    
        if ($line =~ /NO PRIMERS FOUND/) { 
        
            # solution fails primer3
            $ok = 0; 
        } 

        if ($line =~ /^primer3 release/) {

            # done with primer3 for this primer pair
            last PRIMER3_READ;
        }

        # check left and right
        if ($line =~ /^Left.*0$/)  { 
            $ok = 0; 
            $left_bad{$left} = 1 
        }
        if ($line =~ /^Right.*0$/) { 
            $ok = 0; 
            $right_bad{$right} = 1 
        }
    }
        
    # close intermediate files under windows
    if ($^O =~ /^MSW/) {
       close $fh_write;  
       close $fh_read;
    }
    
    # whew! a solution
    if ($ok == 1) { 
        return $p3_reply; 
    }

    # no solution
    return undef;
}



#### create the primers ####
sub _create_primers {

    my ($re,            # restriction enzyme sequence
        $gene,          # ATGC sequence modified to match RE
        $primers_ref,   # reference to primers array
        $clamp          # both left right
       ) = @_;

    my @qs = ( ['', ''], ['?', ''], ['', '?'], ['?', '?'] );
    
    # padding outside of RE
    for my $pad (3 .. 12) {
    
        # ? marks for different types of matching
        foreach my $q (@qs) {
        
            # left and right '?' for regular expression matches
            my $l = $q->[0];
            my $r = $q->[1];
            
            # establish proper pattern
            my $pattern;

            if ($clamp eq 'both')  {            
                $pattern = "[GC](.{3,$pad}$l)$re(.{3,$pad}$r)[GC]" 
            }
            if ($clamp eq 'left')  { 
                $pattern = ".(.{3,$pad}$l)$re(.{3,$pad}$r)[GC]" 
            }
            if ($clamp eq 'right') { 
                $pattern = "[GC](.{3,$pad}$l)$re(.{3,$pad}$r)." 
            }

            # pattern matches
            while ($gene =~ /($pattern)/g) { 
    
                my $location = $-[0]; 
                my $primer   = $1;
                
                # limit primer sizes
                my $l = length $1;
                if ($l >= 18 and $l <= 24) {               
                    $primers_ref->{$primer} = $location; 
                }
            }
        }
    }

    return undef;
}



#### sanity check gene ####
#### see if stop codon has been 'wacked' ####
sub _sanity_check_gene {

    my ($gene,         # ATGC as primers will make
        $gene_start,   # start codon
        $gene_stop     # stop codon
       ) = @_;

    my $stop = _stop_codon($gene, $gene_start);
    
    if (not defined $stop) {     
        return 0 
    }
    
    if ($gene_stop == $stop) {
        return 1;

    } else {
        return 0; 
    }
}



#### generate primer pairs, then process them one by one ####
sub _primer_pairs {

    my ($left_primers_ref,   # array reference
        $right_primers_ref,  # array reference
        $gene,               # ATGC
        $gene_start,         # start codon location
        $gene_stop,          # stop codon location
        $search_shift        # shift search in from right
       ) = @_;

    if (@{$left_primers_ref} == 0 or @{$right_primers_ref} == 0) {

        # need both left and right or go home
        return undef 
    }

    # lefties
    foreach my $left_pr (@{$left_primers_ref}) { 

        # righties
        foreach my $right_pr (@{$right_primers_ref}) { 

            # sequence to be made OK
            if (_sanity_check_gene($gene, $gene_start, $gene_stop) == 1) { 
                
                # primers OK
                my $reply = _primers_ok($left_pr, 
                                        $right_pr, 
                                        $gene, 
                                        $gene_start, 
                                        $gene_stop, 
                                        $search_shift);

                if (defined $reply) { 

                    # reply is OK here
                    return $reply 
                }
            }
        }
    }

    return undef;
}



#### how many '.' in pattern ####
sub _number_dots {

    my (@chars    # array of characters in pattern
       ) = @_;
       
    my $num = 0;

    foreach (@chars) { 
        if ($_ eq '.') { 
            $num += 1 
        } 
    }

    return $num;
}


#### see if there are too many substitutions in a row being requested ####
sub _too_many_substitutions_in_a_row {

    my ($max_dots,    # maximum nuber of dots
        @chars        # array of characters to check
       ) = @_;

    my $n_in_a_row = 0;

    # count '.' in a row, reset where needed
    foreach my $c (@chars) {

        if ($c eq '.') { 
            $n_in_a_row += 1 
        } else { 
            $n_in_a_row = 0 
        }

        if ($n_in_a_row == $max_dots) {   
            return 1 
        }
    }

    return 0;
}



#### recursively generate patterns ####
sub _generate_re_patterns {

    my ($max_dots,     # maximum number of dots also limits recursion
        $pattern_hr,   # pattern hash reference
        @r             # incoming pattern to modify
       ) = @_;  

    my @s;   # next pattern

    # add to hash, keep track of only one time it's found !!!!!!
    $pattern_hr->{ join '', @r } = '';
   
    # limit to annealing capability of primer
    if (_number_dots(@r) == $max_dots)  {   
        return undef 
    }

    if (_too_many_substitutions_in_a_row($max_dots, @r) == 1) {     
        return undef 
    }

    # successively generate next group of patterns
    for my $i (1 .. @r) { 
    
        # already have a '.' here
        if ($r[$i-1] eq '.') {      
            next 
        }
        
        # empty @s, generate clean pattern array
        while (@s) {        
            pop @s 
        }
    
        # build next pattern
        for my $j (1 .. @r) {
        
            if ($i == $j) {             
                push @s, '.'                
            } else {            
                push @s, $r[$j-1] 
            }
        }
            
        # keep going
        _generate_re_patterns( $max_dots, $pattern_hr, @s ); 
    }

    # all patterns stored in hash
    my @k=keys %{$pattern_hr};     
    
    return @k;
}



#### regular expression patterns with '.' generated ####
sub _get_re_patterns {

    my ($re,         # restriction enzyme
        $max_dots    # maximum number of '.'
       ) = @_;     
    
    my @re = split '', $re; # characters in RE
    my %l;                  # hash function with list of generated RE

    # generate patterns for requested enzyme
    my @pats = _generate_re_patterns($max_dots, \%l, @re);

    # sort patterns here
    my @sorted;

    foreach my $n (0 .. $max_dots) {
        foreach my $p (@pats) { 
            if (_number_dots(split '', $p) == $n) { 
                push @sorted, $p 
            }
        }
    }

    return @sorted;
}

    

#### get matches for re pattern in gene ####
sub _get_re_matches {

    my ($gene,       # ATGC sequence
        @patterns    # array of RE patterns
       ) = @_;

    my @positions;

    # loop through patterns
    my $p;
    foreach my $p (@patterns) {

        # loop across gene
        while ($gene =~ /($p)/g) {
            
            # 'magic' - no need to check for in-frame
            push @positions, $-[0]; 
        } 
    }   

    return @positions;
}



#### create left primers ####
sub _create_left_primers {
    
    my ($modified_gene,   # ATGC modified for left RE
        $re_pos_l,        # position of left RE
        $re_l,            # left RE
        $clamp            # type of GC clamp
       ) = @_;
    
    my $left_primers  = {};
    
    if ($clamp eq 'both') {            
        _create_primers($re_l, $modified_gene, $left_primers, 'both');                
    } else {            
        _create_primers($re_l, $modified_gene, $left_primers, 'left');
    }
    
    my @left;

    foreach my $l (keys %{$left_primers}) {    
    
        if (_check_3prime_gc($l) == 1) {
            push @left, $l
        }
    }
    
    return @left;
}



#### create right primers ####
sub _create_right_primers {
    
    my ($more_modified_gene,   # ATGC modified for left and right RE
        $re_pos_r,             # position of right RE
        $re_r,                 # right RE
        $clamp                 # type of GC clamp
       ) = @_;

    my $right_primers = {};
                   
    if ($clamp eq 'both') {                    
        _create_primers($re_r, $more_modified_gene, $right_primers, 'both');
    } else {                    
        _create_primers($re_r, $more_modified_gene, $right_primers, 'right');
    }
                    
    my @right;
     
    foreach my $r (keys %{$right_primers}) {  
        push @right, $r 
    }           
    
    # reverse complement right primer
    my @rev_comp;
    foreach my $r (@right) {
        $r =~ tr/ATGC/TACG/; 
        $r = reverse $r;
        if (_check_3prime_gc($r) == 1) {
            push @rev_comp, $r;
        }
    }
    
    return @rev_comp;
}



#### check 3' end for undesirable [GC] run ####
sub _check_3prime_gc {

    my ($primer   # 5' to 3' order for left or right
       ) =@_;

    my $num_at_end = 5;
    my $end = substr($primer, (length $primer) - $num_at_end, $num_at_end);

    if ($end =~ /[GC][GC][GC]/) { 
            
        # undesirable GC run found at 3' end   
        return 0 
        
    } else { 
  
        return 1 
        
    }
}



#### find in-frame start codon ####
sub _start_codon {

    my ($gene      # ATGC sequence
       ) = @_;

    my $gene_start = 0;
    
    if ($gene =~ /^((.{$CODON_SIZE})*?)(ATG)/) { 
        $gene_start = $-[3];
    }

    return $gene_start;
}



#### find in-frame stop codon location ####
sub _stop_codon {

    my ($gene,         # ATGC sequence
        $gene_start    # look for stop after start
       ) = @_;

    my $WAY_TOO_BIG = 100000000; # bigger than any anticipated pattern
    my $gene_stop = $WAY_TOO_BIG;

    # look for stop codon, keep track of first one in sequence after start codon
    foreach my $stop_codon (('TAA', 'TAG', 'TGA')) {

        if (substr($gene, $gene_start) =~ /^((.{$CODON_SIZE})*?)($stop_codon)/) {

            if ($-[3] < $gene_stop) { 
                $gene_stop = $-[3] 
            }
        }
    }

    # sanity check if stop codon was found
    if ($gene_stop == $WAY_TOO_BIG) { 
        return undef 
    } else { 
        return $gene_stop + $gene_start
    }
}


1;

__END__

=head1 NAME

Bio::MCPrimers and mcprimers
 
=head1 DESCRIPTION

Creates molecular cloning PCR primer pairs for a given gene so that the
gene can be directionally inserted into a vector. pET-32a is the only
vector currently supported. Solver is generic, restriction enzymes and 
their order in the vector are specified in the caller.
 
=head1 EXPORT SUBROUTINES

sub find_mc_primers

    $gene,      # ATGC string with extended regions left and right
    $flag_hr,   # anonymous hash reference to flags
    @re         # array of restriction enzyme strings
    
Not explicitily exported. Use Bio::MCPrimers::find_mc_primers

See mcprimers for an example of use.
 
=head1 INSTALLATION

MCPrimers.pm  - place into Perl Bio/MCPrimers.pm.
mcprimers     - place in a directory where it can be accessed by users.
mcprimers.bat - for MS Windows

PRIMER3DIR -   set this environment variable to point to Primer3 executable.
If PRIMER3DIR is not set, MCPrimers will set it to '.' by default.

MSWindows -    use primer3.exe
Other OS  -    use primer3_core, have IPC::Open3 in search path
 
=head1 DEPENDENCIES

Non-MSWindows use requires IPC::Open3

Primer3 used as primer3.exe on MSWindows and as primer3_core elsewhere.
Specify PRIMER3DIR for path to Primer3 executable

=head1 EXAMPLES

mcprimers -h

mcprimers ppib.fa 

mcprimers -clamp 3prime hemy.fa 

mcprimers -s 24 -all -clamp 3prime cyss.fa > cyss.pr3

See mcprimers for an example of the use of Bio::MCPrimers
 
=head1 LIMITATIONS

Limitation: Does not account for redundancy codes.

Note:       MS Windows runs use intermediate files, so don't use this in a
            manner that intermediate files will overwrite one another. CGI
            use under MSWindows is probably a bad idea.
            
There is no guarantee this code will find the best solution, or even any 
solution, or that the solutions it finds will be correct or useful to you.

Use at your risk. Check any solutions you obtain.
 
=head1 BUGS

Probably. Use at your own risk.

This software comes with no guarantee of usefulness.
Use at your own risk. Check any solutions you obtain.
Stephen G. Lenk assumes no responsibility for the use of this software.
 
=head1 SOFTWARE HEURISTIC EMPLOYED

- User of mcprimers gets FASTA for gene that is in-frame. The Kegg
site allows a FASTA with 21 NT upstream and 200 or so NT downstream to
be generated. This FASTA is specified in command line for mcprimers.

- mcprimers generates an array of ATGC representing extended gene.

- mcprimers generates an in-order array of restriction enzyme sequences.

- mcprimers generates a hash of flag values.

- MCPrimers gets extended gene sequence, RE array, and flag hash.

- An array of permutations of each RE with up to three '.', not all in a
row, is generated for regular expressions.

- An array of matches to the gene for each permutation is generated, using
regular expessions that generate variability and proper GC clamping. These
matches are checked to be sure that there are no extended GC runs (3 or more 
in a row) in the last 5 nucleotides at the 3' end.

- The restriction enzyme sequences and regular expression matches are 
processed so that the first site is processed against all remaining sites.
When that is done, the second site is processed against all remaining
sites. This continues until the last site is left, where processing
terminates.

- The current head against all the rest process generates successive pairs
of primers. These are tested against Primer3. Any primers identified by 
Primer3 as bad are excluded from future consideration, to reduce search
space size.

-The extended gene sequence is successively modified to match the restriction
enzyme sequences, allowing proper matching in Primer3.

- A solution primer pair validated by Primer3 is checked to ensure that it 
does not match anywhere else in the sequence.

- The modified gene sequence is checked to ensure that the STOP codon is 
still in the same place, which is needed when the protein is expressed in
the transfected vector.

- If only one solution is desired, the first solution is returned to the 
calling program as soon as it is found.

- If all solutions are desired, they are successively found and retained, 
being returned to the caller when the whole search is done.

- Individual solutions are held in anonymous hashes. The has references are
held in an array, which is returned to the caller for processing.
 
=head1 COPYRIGHT

Stephen G. Lenk (C) 2005. All rights reserved. 

This program is free software; you can redistribute it and/or  
modify it under the same terms as Perl itself.
 
Primer3 is used by this code to verify that the PCR primers are OK.
Primer3 is Copyright (c) 1996,1997,1998,1999,2000,2001,2004
Whitehead Institute for Biomedical Research. All rights reserved.

=head1 AUTHOR

Stephen G. Lenk, November 2005 

slenk@emich.edu
 
=head1 ACKNOWLEDGEMENTS

Primer3 is called by this code to verify that the PCR primers are OK.

Thanks to Tim Wiggin for algorithm suggestions and encouragement
    Use of direct string comparisons
    Modify gene for Primer3 check
    
Thanks to Dan Clemans for showing me molecular cloning in the first place
I am using Dr. Clemans's ideas about good MC primers in this code.
Any errors in interpretation or implementation are mine.

Ken Youens-Clark <kyclark@gmail.com> has provided guidance in the
proper behaviour and naming of this software, including a code review.

Other references:

http://www.premierbiosoft.com/tech_notes/PCR_Primer_Design.html
http://www.mcb.uct.ac.za/pcroptim.htm


=cut