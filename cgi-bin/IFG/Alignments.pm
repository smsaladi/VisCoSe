###################################################################################################
# IFG::Alignments.pm
# (c) 2001/2002/2003 Michael Spitzer, IFG Muenster (Integrated Funtional Genomics), Germany
#
# several subroutines that may be useful someday, somehow, for someone...
# update: ok, they are indeed useful, especially for me... ;-)
#
# 06.08.2003   added delete_seqs()
# 12.01.2003   added split_alignment()
# 04.12.2002   added write_phylip()

use strict;
use diagnostics;

package IFG::Alignments;
require Exporter;

our @EXPORT  = qw(read_aln                      write_aln                        write_phylip                     write_raw_aln                 find_lowercase
                 _lowercase2gap                 _find_gap_cols                   find_unaligned_columns           find_gap_columns              round sort_num
                 degap_alignment                uppercase_seqs                   mutate_alignment                 split_alignment               gap2x
                 delete_seqs                    add_text                         add_taxon                        cut_range_from_aln            despace
                 find_gap_only_columns          remove_columns                   merge_alignments                 run_mafft                     get_first_x_sequences
                 calc_matches                   calc_pairwise_distance           jukes_cantor_correction          kimura_protein_correction     pairwise_distance_matrix
                 pretty_print_seqpair           run_mview                        fasta2stockholm
                );
		  
our $VERSION = "06.08.2003";

use IFG::NutsnBolts;
use Data::Dumper;

use constant max_defline_len  => 60;

sub sort_num { $a <=> $b }

sub read_aln {
# reads a FASTA-formatted alignment from a file and returns it as a hash, together with an array holding the FASTA-codes in proper sequence
# Each defline has to be UNIQUE, empty deflines ">" are not allowed (and senseless in our context, anyway... ;-)
   my $aln              =   $_[0];                                               # filename
   #my $line             = "";
   my (%align, @codearray, %deflines)  = ((), (), ());                           # "%deflines" is used to check for double (triple, etc...) deflines
   
      
   open (_aln_, "<$aln") or die ">>> IFG::Alignments::read_aln -> [IOerr] Could not open '$aln'!\n";
   until (eof _aln_) {
      my $line = <_aln_>;
      #chomp($line);
      $line =~ s/\n|\r|\t//g;
      my $defline = '';
      if ($line =~ s/\A(?:\>|\%)//) {
         if ($deflines{$line}) {
            $defline = $line.'_copy'.($deflines{$line} + 1);                                                                                            # if the defline already exists, append a number to the current defline
            print "\n>>> [IFG::Alignments::read_aln()] WARNING: found multiple (n=".($deflines{$line} + 1).
                  ") sequences with defline '$line' => renaming sequence to '$defline'\n";
         }
         else { $defline = $line; }
         push(@codearray, $defline);                                                                                                         # if we hit a fasta defline we remove the ">" (or "%") and add it to the codearray
         if ($deflines{$line}) { $deflines{$line}++; }
         else { $deflines{$line} = 1; }
      }
      else { ($align{$codearray[-1]} .= $line) =~ s/\s//g; }                                                                                 # $codearray[-1] is always the last fasta-code that was 'push'ed! Read in all sequence lines and clean them from whitespaces
   }
   close(_aln_);
   return \%align, \@codearray;
}

sub write_aln {
# writes an alignment as FASTA format
   my $filename     = $_[0];                                                              # filename to write the alignment to
   my $alignment    = $_[1];                                                              # hash reference to the alignment (derived by IFG::read_aln)
   my $codearray    = $_[2];                                                              # array reference to the array holding the fasta-codes of the alignment (derived by IFG::read_aln)

   open(_aln_, ">$filename") or die ">>> IFG::Alignments::write_aln -> [IOerr] Could not open '$filename' for write-access!\n";
   foreach my $code (@{$codearray}) { print _aln_ ">$code\n".$alignment->{$code}."\n"; }
   close(_aln_);
}

sub write_phylip {
# writes an alignment as PHYLIP format
   my $filename               = $_[0];                                                    # filename to write the alignment to
   my $alignment              = $_[1];                                                    # hash reference to the alignment (derived by IFG::read_aln)
   my $codearray              = $_[2];                                                    # array reference to the array holding the fasta-codes of the alignment (derived by IFG::read_aln)
   my $len                    = length($alignment->{$codearray->[0]});
   my ($start, $end, $round)  = (0, 0, 0);
   my $line                   = "";

   open(_aln_, ">$filename") or die ">>> IFG::Alignments::write_phylip -> [IOerr] Could not open '$filename' for write-access!\n";
   print _aln_ scalar(@{$codearray})." $len\n";
   while ($end < $len) {                                                                  # go thru the complete sequence-length until the end is reached
      foreach my $code (@{$codearray}) {                                                  # write a subset of each sequence for the actual block
         if ($round == 0) {                                                               # in the zeroeth round: write the deflines!
            if (length($code) <= 10) { $line = ($code.(' ' x (11-length($code)))); }      # defline-length <= 10: eventually fill up with spaces
            else { $line = substr($code, 0, 10).' '; }                                    # defline-length >  10: use 1-10-substring
         }
         else { $line = ' ' x 11; }                                                       # not round zero: fill up with 10+1 spaces
         if ($len - $start >= 50) { $end = $start + 50; }                                 # print 50 chars of the sequence or up to it's end
         else { $end = $len; }
         print _aln_ $line.substr($alignment->{$code}, $start, ($end-$start))."\n";       # print to file
      }
      if ($end < $len) { print _aln_ "\n"; }                                              # insert a blank line if needed
      $round++;                                                                           # increase the round-counter
      $start += 50;                                                                       # increase the starting position by 50
   }
   close(_aln_);
}

sub write_raw_aln {
# writes an alignment as RAW (easy viewing with any text editor or browser)
   my $filename     = $_[0];                                                              # filename to write the alignment to
   my $alignment    = $_[1];                                                              # hash reference to the alignment (derived by RIPE_ALN::read_aln)
   my $codearray    = $_[2];                                                              # array reference to the array holding the fasta-codes of the alignment (derived by RIPE_ALN::read_aln)
   my $codelen      = 0;
   my $dummy        = "";

   for (my $i = 0; $i < scalar(@{$codearray}); $i++) {                                             # get length of the longest defline
      if (length($codearray->[$i]) > $codelen) { $codelen = length($codearray->[$i]); }
   }
   open(_aln_, ">$filename") or die ">>> IFG::Alignments::write_raw_aln -> [IOerr] Could not open '$filename' for write-access!\n";
   foreach my $code (@{$codearray}) {
      if (length($code) < $codelen) { $dummy = $code.(" " x ($codelen - length($code))); }         # eventually fill the defline up with whitespaces if it's to short
      else { $dummy = $code; }
      print _aln_ "$dummy:   ".$alignment->{$code}."\n";
   }
   close(_aln_);
}

sub _lowercase2gap {
# OBSOLETE???????
# internal subroutine, changes lowercase characters (= unaligned columns from dialign) into gap characters
   my $pep       = $_[0];
   my $pep_codes = $_[1];
   
   foreach my $code (@{$pep_codes}) { $pep->{$code} =~ s/[a-z]/-/g; }                     # change all lowercase characters to gap-characters ('-')
   # no returning necessary since we modify via *references*
}

sub _find_gap_cols {
# OBSOLETE???????
# internal subroutine
# searches for columns having *only* gap characters in them and pushes these positions into a handed array (array_ref '$cut_me')
   my $pep       =   $_[0];                                                               # hash reference to dialign alignment
   my @pep_codes = @{$_[1]};                                                              # array reference to fasta-codes of dialign alignment (dereferenced for performance reasons)
   my $cut_me    =   $_[2];                                                               # we store the 'to-be-deleted'-columns here
   my $len       = length($pep->{$pep_codes[0]});                                         # length of alignment
   my $no_gap    = 0;                                                                     # used as 'lowercase'-indicator
   
   for (my $i = 0; $i < $len; $i++) {                                                     # step thru the columns
      $no_gap = 0;
      foreach my $code (@pep_codes) {                                                     # step thru the sequences
         $no_gap = (substr($pep->{$code}, $i, 1) ne "-");                                 # here we look if any sequence has anything other than gaps in column '$i'
	 last if ($no_gap);                                                               # if we find anything other than a gap, we can immediately leave the inner 'foreach'-loop!
      }
      unless ($no_gap) { push(@{$cut_me}, $i); }                                          # if we found a gap_only-column, we add it to the array
   }
   # no returning necessary since we modify via *references*
}

sub find_unaligned_columns {
# OBSOLETE???????
# this subroutine searches an alignment for columns having lowercase characters or consisting completely of gap-symbols ('-').
# It returns a reference to an array with the found columns.
   my $pep        = $_[0];                                                                # hash reference to dialign peptide alignment
   my $pep_codes  = $_[1];                                                                # array reference which holds the fasta-codes of the peptide-alignment as they occured in the file
   my @cut_me     = ();
   
   &_lowercase2gap($pep, $pep_codes);
   &_find_gap_cols($pep, $pep_codes, \@cut_me);
   @cut_me = sort(sort_num @cut_me);
   return \@cut_me;
}

sub find_gap_columns {
# OBSOLETE???????
# this subroutine searches an alignment for columns having *one or more* gap characters
# It returns a reference to an array with the found columns.
   my $pep        =   $_[0];                                                              # hash reference to dialign peptide alignment
   my @pep_codes  = @{$_[1]};                                                             # array reference which holds the fasta-codes of the peptide-alignment as they occured in the file
   my $len        = length($pep->{$pep_codes[0]});                                        # length of alignment
   my @cut_me     = ();
   my $gap        = 0;

   for (my $i = 0; $i < $len; $i++) {                                                     # step thru the columns
      undef($gap);
      foreach my $code (@pep_codes) {                                                     # step thru the sequences
         $gap = (substr($pep->{$code}, $i, 1) eq "-");                                    # here we look if any sequence has a gap in column '$i'. One gap is enough to qualify for deletion!
        last if ($gap);                                                                   # if we find *one* gap, we can immediately leave the 'foreach'-loop!
      }
      if ($gap) { push(@cut_me, $i); }                                                    # if we found a gap-column, we add it to the array
   }
   @cut_me = sort(sort_num @cut_me);                                                      # sort the array position-wise
   return \@cut_me;
}

sub find_lowercase {
# searches for columns having *one or more* lowercase characters in them and pushes these positions into a handed array (array_ref '$cut_me')
   my $pep       =   $_[0];
   my @pep_codes = @{$_[1]};
   my $cut_me    =   $_[2];
   my $len       = length($pep->{$pep_codes[0]});                                         # length of alignment
   my $lower     = 0;                                                                     # used as 'lowercase'-indicator
   
   for (my $i = 0; $i < $len; $i++) {                                                     # step thru the columns
      $lower = 0;
      foreach my $code (@pep_codes) {                                                     # step thru the sequences
         $lower = (substr($pep->{$code}, $i, 1) =~ m/[a-z]/);                             # here we look if any sequence has lowercase chars in column '$i'
        last if ($lower);                                                                 # if we find one, we can immediately leave the 'foreach'-loop!
      }
      if ($lower) { push(@{$cut_me}, $i); }                                               # if we found a lowercase-column, we add it to the array
   }
}

sub degap_alignment {
# removes all gap characters from an alignment resulting in a virtually unaligned sequence set
   my $aln                    = $_[0];
   my $codes                  = $_[1];
   my $maxseqs                = $_[2];                                                    # will only keep $maxseqs number of sequences
   my (%degap, @degapcodes)   = ((), ());
   my $numseqs                = scalar(@{$codes});

   if (($maxseqs) && ($numseqs > $maxseqs)) { $numseqs = $maxseqs; }                      # only convert and return '$maxseqs' sequences, if specified, and drop all sequences >$maxseqs
   for (my $i = 0; $i < $numseqs; $i++) {
      ($degap{$codes->[$i]} = $aln->{$codes->[$i]}) =~ s/-//g;                            # copy the sequence data to %degap and remove all gap characters
      push(@degapcodes, $codes->[$i]);
   }
   return (\%degap, \@degapcodes);
}

sub uppercase_seqs {
# converts a complete alignment ($in, $in_codes) into uppercase-only sequences
   my $in         = $_[0];
   my $in_codes   = $_[1];

   foreach my $code (@{$in_codes}) { $in->{$code} = uc($in->{$code}); }
   return $in;                                                                            # returning $in is not really nessecary, but done here anyway
}

sub mutate_alignment {
# (very) simple mutation of protein sequences
   my $in               = $_[0];
   my $in_codes         = $_[1];

   foreach my $code (@{$in_codes}) { $in->{$code} =~ tr/ILVFMCAGPTSYWQNHEDKR/RKDEHNQWYSTPGACMFVLI/; }
   return $in;
}

sub split_alignment {                                    # divide alignments into portions (column-wise)
# splits an alignment at the given position and returns both subsets
   my $in                  = $_[0];
   my $in_codes            = $_[1];
   my $splitpos            = $_[2] + 1;                                                   # "convert" the split position to Perl-like counting (beginning with zero!)
   my (%split1, %split2)   = ((), ());
   my $aln_len             = length($in->{$in_codes->[0]});

   if ($splitpos < $aln_len) {                                                            # don't split if the split-position is bigger than the alignment-length
      foreach my $code (@{$in_codes}) {
         $split1{$code} = substr($in->{$code}, 0, $splitpos);
         $split2{$code} = substr($in->{$code}, ($splitpos + 1), $aln_len);
      }
   }
   return (\%split1, \%split2);
}

sub gap2x {
# converts gap characters to 'X'
   my $in                  = $_[0];                                                       # hash_ref to alignments
   my $in_codes            = $_[1];                                                       # array_ref to deflines

   foreach my $code (@{$in_codes}) { $in->{$code} =~ tr/-/X/; }
   return $in;
}

sub delete_seqs {
# deletes sequences (specified in @{$del}) from an alignment
   my $in                  = $_[0];                                                       # hash-ref to array
   my $in_codes            = $_[1];                                                       # array-ref to deflines of array
   my $del                 = $_[2];                                                       # array-ref to deflines of sequences which are to be deleted
   my $param               = $_[3];
   my @temp                = ();
   my $found               = 0;

   foreach my $code (@{$del}) { delete($in->{$code}); }                                   # delete key->value-pairs from alignment-hash
   foreach my $code (@{$in_codes}) {                                                      # generate new defline-array without deleted sequences
      $found = 0;
      foreach my $delseq (@{$del}) {
         $found = ($code eq $delseq);
         last if $found;
      }
      unless ($found) { push(@temp, $code); }                                             # only copy $code if it doesn't exist in @del
      else {
         if ($param->{'verbose'}) { print "   deleting '$code'\n"; }
      }
   }
   return \@temp;
}

sub add_text {
# adds a user-defined text to each defline of an alignment
   my $in                  = $_[0];
   my $in_codes            = $_[1];
   my $text                = $_[2];
   my $param               = $_[3];
   my (%new, @new_codes)   = ((), ());

   foreach my $code (@{$in_codes}) {
      push(@new_codes, $code.$text);
      $new{$new_codes[-1]} = $in->{$code};
   }
   return (\%new, \@new_codes);
}

sub add_taxon {
# adds a user-defined taxon to each defline of an alignment (after converting all square brackets to round brackets)
   my $in                  = $_[0];
   my $in_codes            = $_[1];
   my $taxon               = $_[2];
   my $param               = $_[3];
   my (%new, @new_codes)   = ((), ());

   for (my $i = 0; $i < scalar(@{$in_codes}); $i++) {
      ($new_codes[$i] = $in_codes->[$i]) =~ tr/[]/()/;
      $new_codes[$i] .= " [$taxon]";
      $new{$new_codes[$i]} = $in->{$in_codes->[$i]};
   }
   return (\%new, \@new_codes);
}

sub cut_range_from_aln {
# cuts a range from an alignment. The range must have the format "xx-yy"
   my $in            = $_[0];
   my $in_codes      = $_[1];
   my $range         = $_[2];
   my @range         = split(/\.\./, $range);
   my %cut           = ();

   foreach my $code (@{$in_codes}) { $cut{$code} = substr($in->{$code}, $range[0], ($range[1] - $range[0])); }
   return \%cut;
}

sub despace {
# changes some special characters (not liked by some alignment programs) into underscores
   my $in                  = $_[0];
   my $in_codes            = $_[1];
   my $param               = $_[2];
   my (%new, @new_codes)   = ((), ());
   my $newcode             = '';

   foreach my $code (@{$in_codes}) {
      ($newcode = $code) =~ s/\s|\:|\\|\//_/g;
      $new{$newcode} = $in->{$code};
      push(@new_codes, $newcode);
   }
   return (\%new, \@new_codes);
}

sub find_gap_only_columns {
# searches for columns having *only* gap characters in them and pushes these positions into a handed array (array_ref '$cut_me')
   my $pep                 =   $_[0];                                # hash reference to dialign alignment
   my @pep_codes           = @{$_[1]};                               # array reference to fasta-codes of dialign alignment (dereferenced for performance reasons)
   my %cut_me              = ();                                     # we store the 'to-be-deleted'-columns here
   my $len                 = length($pep->{$pep_codes[0]});          # length of alignment
   my $no_gap              = 0;                                      # used as 'lowercase'-indicator
   
   for (my $i = 0; $i < $len; $i++) {                                # step thru the columns
      $no_gap = 0;                                                   # init 'nogap' for each column
      foreach my $code (@pep_codes) {                                # step thru the sequences
         $no_gap = (substr($pep->{$code}, $i, 1) ne "-");            # here we look if any sequence has anything other than gaps in column '$i'
	 last if ($no_gap);                                          # if we find anything other than a gap, we can immediately leave the 'foreach'-loop!
      }
      unless ($no_gap) { $cut_me{$i} = 1; }                          # if we found a gap_only-column, we add it to the hash
   }
   return \%cut_me;
}

sub remove_columns {
   my $in                  = $_[0];
   my $in_codes            = $_[1];
   my $columns             = $_[2];
   my $param               = $_[3];             # not needed!!!!!!!!!!
   my $len                 = length($in->{$in_codes->[0]});          # length of alignment
   my %temp                = ();

   for (my $i = 0; $i < $len; $i++) {                                # step thru the columns
      unless ($columns->{$i}) {                                      # if the actual column is *NOT* a gap-only column...
         foreach my $code (@{$in_codes}) {                           # ...step thru the sequences and add the character at position '$i'
            $temp{$code} .= substr($in->{$code}, $i, 1);
         }
      }
   }
   return (\%temp, $in_codes);
}

sub merge_alignments {
   my $aln1                = $_[0];
   my $aln1codes           = $_[1];
   my $aln2                = $_[2];
   my $aln2codes           = $_[3];
   my (%temp, @tempcodes)  = ((), ());

   foreach my $code (@{$aln1codes}) {
      $temp{$code} = $aln1->{$code};
      push(@tempcodes, $code);
   }
   foreach my $code (@{$aln2codes}) {
      $temp{$code} = $aln2->{$code};
      push(@tempcodes, $code);
   }
   return (\%temp, \@tempcodes);
}

sub run_mafft {
   my $source              = $_[0];                # source filename
   my $destination         = $_[1];                # destination filename
   my $mafft_script        = $_[2];                # which mafft script is to be used for aligning (nwns, nwnsi, fftns, *fftnsi*, fftnsi_nj)
   my $iterations          = $_[3];                # only of use when using "fftns" or "nwns"
   my $shellstring         = '';
   
   unless ($mafft_script) { $mafft_script = 'fftnsi'; }                                                                          # if no mafft script was specified, use 'fftns' as default
   unless ($destination)  { $destination  = $mafft_script.'_'.$source; }                                                         # if no destination filename was specified, generate one from the source filename
   if (($mafft_script eq 'fftns') || ($mafft_script eq 'nwns')) {
      unless ($iterations) { $iterations = 2; }
      $shellstring = "$mafft_script $iterations $source >$destination 2>/dev/null";                                              # run mafft (generate shellstring) - NOTE: STDERR is redirected to "/dev/null" !!!
   }
   elsif ($mafft_script eq 'fftnsi_nj') { $shellstring = "fftnsi --nj $source >$destination 2>/dev/null"; }                      # customized fftnsi version, activating NJ instead of UPGMA (personal communication with MAFFT author)
   else { $shellstring = "$mafft_script $source >$destination 2>/dev/null"; }
   my $log = `$shellstring`;                                                                                                     # run mafft (execution)
}

sub divide_alignment {                             # divide alignments into portions (sequence-wise)
   my $codes               = $_[0];
   my $numseqs             = $_[1];
   my %tempcodes           = ();
   my $portion             = 0;
   my @dummy               = @{$codes};

   while (@dummy) { @{$tempcodes{$portion++}} = splice(@dummy, 0, $numseqs); }
   return (\%tempcodes);
}

sub get_first_x_sequences {      # returns the first "$numseqs" sequences from the alignment "$in, $in_codes"
   my $in                  = $_[0];
   my $in_codes            = $_[1];
   my $maxseqs             = $_[2];
   my (%temp, @tempcodes)  = ((), ());
   my $numseqs             = scalar(@{$in_codes});

   if (($maxseqs) && ($numseqs > $maxseqs)) { $numseqs = $maxseqs; }                                                 # only convert and return '$maxseqs' sequences, if specified, and drop all sequences >$maxseqs
   for (my $i = 0; $i < $numseqs; $i++) {
      $temp{$in_codes->[$i]} = $in->{$in_codes->[$i]};
      push(@tempcodes, $in_codes->[$i]);
   }
   return (\%temp, \@tempcodes);
}

sub calc_matches {
   my $seq1                      = $_[0];
   my $seq2                      = $_[1];
   my $len                       = $_[2] || length($$seq1);
   my ($matches, $gaps, $npos)   = (0, 0, 0);

   for (my $i = 0; $i < $len; $i++) {                                                                                # step thru the length of the sequence(s)
      my $char1 = uc(substr($$seq1, $i, 1));                                                                         # extract the characters at position 'i'
      my $char2 = uc(substr($$seq2, $i, 1));                                                                         # extract the characters at position 'i'
      if (($char1 ne '-') && ($char2 ne '-')) {                                                                      # positions featuring gaps (either seq1 or seq2) are not counted as matches, but are separately
         $npos++;                                                                                                    # remember all positions which have amino acids in *BOTH* sequences
         $matches += ($char1 eq $char2);                                                                             # full match for similar amino acids (1)!
         if ((($char1 eq 'B') && ($char2 eq 'D')) || (($char1 eq 'Z') && ($char2 eq 'E'))) { $matches += 0.5; }      # do we have ambiguous amino acids? (B->D / Z->E) -> half a match for ambiguous amino acids (0.5)!
      }
      else { $gaps++; }                                                                                              # count gaps (gap penalty is handled in the parent routine)
   }
   return ($matches, $gaps, $npos);
}

sub calc_pairwise_distance {
   my $matches             = $_[0];
   my $gaps                = $_[1];
   my $npos                = $_[2];
   my $gap_penalty         = $_[3];

   my $divisor = $npos + ($gaps * $gap_penalty);
   if ($divisor == 0) { return 1; }
   else { return (1 - ($matches / $divisor)); }
}

sub jukes_cantor_correction {
   my $pwdist              = $_[0];
   my $b                   = (19 / 20);

   return (-($b) * log(1 - ($pwdist / $b)));                                                                         # Perl's "log" is for "base e" !!!
}

sub kimura_protein_correction {
   my $pwdist              = $_[0];

   return -(log(1 - $pwdist - (0.2 * $pwdist**2)));                                                                  # Perl's "log" is for "base e" !!!
}

sub pairwise_distance_matrix {
   my $in                        = $_[0];
   my $codes                     = $_[1];
   my $param                     = $_[2];
   my $numseqs                   = scalar(@{$codes});
   my $seqlen                    = length($in->{$codes->[0]});
   my $numtrials                 = (($numseqs - 1) * $numseqs) / 2;              # number of comparisons the algorithm has to perform (for percentage display)
   my @matrix                    = ();
   my ($counter, $gap_penalty)   = (0 , 0);
   my $mod                       = 0;

   for (my $i = 0; $i < ($numseqs - 1); $i++) {
      for (my $j = ($i + 1); $j < $numseqs; $j++) {
         if ($param->{'verbose'}) {
            $counter++;
            my $percent = IFG::NutsnBolts::round((($counter / $numtrials) * 100), 0);
            if ($percent == $mod) {
               print '.';
               $mod += 10;
            }
         }
         my ($matches, $gaps, $npos) = calc_matches(\$in->{$codes->[$i]}, \$in->{$codes->[$j]}, $seqlen);            # calculate the matching rate / gap rate of two sequences
         my $pwdist = calc_pairwise_distance($matches, $gaps, $npos, $gap_penalty);                                  # calculate the *uncorrected* pairwise distance of 2 sequences
         #my $jccorr = jukes_cantor_correction($pwdist);                                                              # correct the pairwise distance following "Jukes-Cantor"
         #my $kpcorr = 0;
         #if ($gap_penalty == 0) { $kpcorr = kimura_protein_correction($pwdist); }                                    # correct the pairwise distance following "Kimura for proteins"
         #print "PWDIST: $pwdist\nJCCORR: $jccorr\nKPCORR: $kpcorr\n\n";
         $pwdist = IFG::NutsnBolts::round($pwdist, 10);
         $matrix[$i]->[$j] = $pwdist;                                                                                # seq_i is to seq_j    what    seq_j is to seq_i
         $matrix[$j]->[$i] = $pwdist;                                                                                # ...and vice-versa
         
      }
      $matrix[$i]->[$i] = IFG::NutsnBolts::round(0, 10);                                                             # the pairwise distance of a sequence to itself is ZERO!
   }
   return \@matrix;
}

sub pretty_print_seqpair {
   my $seq1                = $_[0];                   # seq1 (reference!)
   my $seq2                = $_[1];                   # seq2 (reference!)
   my $code1               = $_[2];                   # defline1 (reference!)
   my $code2               = $_[3];                   # defline2 (reference!)
   my $width               = $_[4];
   my $offset              = $_[5];
   my $max_defline_len     = $_[6];
   my $param               = $_[7];
   my ($code11, $code22)   = ($$code1, $$code2);      # used for making the deflines equal length

   if ($offset > 0) {                                                                                                                  # if an offset is wanted of '$offset' characters => insert '$offset' space chars
      $code11 = (' ' x $offset).$code11;
      $code22 = (' ' x $offset).$code22;
   }
   if (length($code11) != length($code22)) {                                                                                           # make the two defline equal length by appending space characters
      if (length($code11) > length($code22)) {
         my $diff = (length($code11) - length($code22));
         $code22 .= ' ' x $diff;
      }
      else {
         my $diff = (length($code22) - length($code11));
         $code11 .= ' ' x $diff;
      }
   }
   if (length($code11) > $max_defline_len) {                                                                                            # if defline length exceeds 'max_defline_len', cut it to 'max_defline_len'
      substr($code11, ($max_defline_len - 3), (length($code11) - $max_defline_len - 3), '...');
      substr($code22, ($max_defline_len - 3), (length($code22) - $max_defline_len - 3), '...');
   }
   my $len  = $width - length($code11);                                                                                                # both deflines have the same length at this point!
   my @seq1 = grep(/./, split(/(.{$len})/, $$seq1));                                                                                   # split sequence in parts of "$len" (or less)
   my @seq2 = grep(/./, split(/(.{$len})/, $$seq2));                                                                                   # split sequence in parts of "$len" (or less)
   print "\n\n\n";
   for (my $i = 0; $i < scalar(@seq1); $i++) {
      my $sims = '';
      for (my $x = 0; $x < length($seq1[$i]); $x++) {
         if ((substr($seq1[$i], $x, 1) eq substr($seq2[$i], $x , 1)) && (substr($seq1[$i], $x, 1) ne '-')) { $sims .= '|'; }           # only set marks for identical amino acids (not for positions with where both feature gaps!)
         elsif ((substr($seq1[$i], $x, 1) ne substr($seq2[$i], $x , 1)) &&
                (substr($seq1[$i], $x, 1) ne '-') && (substr($seq2[$i], $x, 1) ne '-')) { $sims .= '#'; }                              # only set marks for non-identical amino acids (not for positions with one bearing gaps!)
         else { $sims .= ':'; }
      }
      print $code11.' :   '.$seq1[$i]."\n".(' ' x (length($code11) + 5)).$sims."\n".$code22.' :   '.$seq2[$i]."\n\n";
   }
   print "\n\n\n";
}

sub run_mview {
   my $blastfile                 = $_[0];
   my $mviewfile                 = $_[1];
   my $param                     = $_[2];

   if (!(-s $mviewfile) || ($param->{'force_pipeline'} && ($param->{'force_pipeline'} == 2))) {                                                          # only run MVIEW if the MVIEW-alignment does not exist!
      if (-s $blastfile) {                                                                                                                               # if BLAST report exists and has a size greater than zero...
         my $shellstring = ($param->{'mviewdir'} ? $param->{'mviewdir'}.'/' : '').'mview -in blast -out pearson -hsp discrete '.$blastfile.' >'.$mviewfile;                               # generate MVIEW commandline
         my $log = `$shellstring`;                                                                                                                       # run MVIEW
      }
      else {
         if ($param->{'verbose'}) { print "[IFG::Alignments::run_mview()] BLAST report ('$blastfile')is empty!!!\n"; }                                   # BLAST report is empty, there was some error while executing BLAST
         exit 0;
      }
   }
   else {
      if ($param->{'verbose'}) { print "[IFG::Alignments::run_mview()] Skipping MVIEW, alignment '$mviewfile' already exists!!!\n"; }
   }
}

sub fasta2stockholm {
   my $infile  = $_[0];
   my $outfile = $_[1];
   my $param   = $_[2];

   if (-s $infile) {
      unless (-s $outfile) {
         my $shellstring = $param->{'sreformatdir'}."/sreformat --pfam selex $infile >$outfile";
         my $log = `$shellstring`;
      }
      return 1;
   }
   else { return -1; }
}
