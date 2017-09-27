###################################################################################################
# motif.pm
# (c) 2002 Michael Spitzer, IFG Muenster (Integrated Funtional Genomics), Germany
#

package VisCoSe::Motif;
require Exporter;

use strict;
use diagnostics;

our @EXPORT  = qw(search_motifs                       search_motifs_including_gaps                 _get_seeds                          _walk_sequence
                  _print_motif                        _refine_motifs                               _build_patterned_motifs             eval_motifs
                  write_motifs                        evaluate_conservation_rates                  keep_sequences_with_motifs          keep_sequences_with_motifs_for_conserved
                  print_conserved_positions           calc_mean_conservation_rate_per_file);

our $VERSION = "23.09.2002";

sub sort_num { $a <=> $b }

sub _print_raw_motif_results {
   my $conserve               = $_[0];
   my $max                    = $_[1];
   my $param                  = $_[2];

   if ($param->{"verbose"} && ($param->{"verbose"} == 2)) {                                                                      # only print when '-vv' was specified
      foreach my $column (sort(sort_num keys(%{$conserve}))) {
         print "[$column]\t";
         foreach my $symbol (sort(keys(%{$conserve->{$column}}))) { print "$symbol\t"; }
         print "\n\t";
         foreach my $symbol (sort(keys(%{$conserve->{$column}}))) { print $conserve->{$column}->{$symbol}."\t"; }
         print "\n\n";
      }
      print "Maximum rate found: $max\n";
   }
}

sub _print_rated_motif_results {
   my $rated                  = $_[0];
   my $param                  = $_[1];

   if ($param->{"verbose"} && ($param->{"verbose"} == 2)) {                                                                      # only print when '-vv' was specified
      foreach my $rate (sort(sort_num keys(%{$rated}))) {
         print "Symbols with conservation $rate:\n";
         foreach my $column (sort(sort_num keys(%{$rated->{$rate}}))) {
            if ($rated->{$rate}->{$column} ne "-") { print $rated->{$rate}->{$column}."$column\t"; }
         }
         print "\n\n";
      }
   }
}

sub search_motifs {
   my $clean                     =   $_[0];
   my @cleancodes                = @{$_[1]};
   my $param                     =   $_[2];
   my $len                       = length($clean->{$cleancodes[0]});
   my %conserve                  = ();
   my %rated                     = ();
   my %cont                      = ();
   my ($max, $rate)              = (0, 0);
   my ($symbol, $not_gap_only)   = ("", 0);

   $max = scalar(@cleancodes);                                                                                                   # maximum: number of sequences == number of times an AA may occur per column
   for (my $i = 0; $i < $len; $i++) {                                                                                            # step thru length of alignment                                      
      $not_gap_only = 0;
      foreach my $code (@cleancodes) {                                                                                           # step thru sequences of alignment                                   
         $symbol = substr($clean->{$code}, $i, 1);
         $not_gap_only += ($symbol ne '-');                                                                                      # will be TRUE if there occurs a non-gap character!                  
         if ($symbol ne '-') {                                                                                                   # don't count gaps                                                   
            $conserve{$i}->{$symbol}++;                                                                                          # per column: increase counter per aminoacid                         
            #if ($conserve{$i}->{$symbol} > $max) { $max = $conserve{$i}->{$symbol}; }                                            # maximum-determination (I think this way was wrong, dude!!!)
         }
      }
      unless ($not_gap_only) { $conserve{$i}->{'-'} = 1; }                                                                       # do we have a gap-only column?                                      
   }
   &_print_raw_motif_results(\%conserve, $param);                                                                                # print the results of the raw motif search method                   
   foreach my $column (sort(sort_num keys(%conserve))) {                                                                         # REFINEMENT: convert into percentages                               
      foreach my $symbol (sort(keys(%{$conserve{$column}}))) {                                                                   # for each column look at all found aminoacids                       
         $rate = &IFG::NutsnBolts::round((($conserve{$column}->{$symbol} / $max) * 100), 2);                                     # and calculate the percentage for number of sequences it occurs in  
         $conserve{$column}->{$symbol} = $rate;                                                                                  # save the percentage back                                           
      }
   }
   foreach my $column (sort(sort_num keys(%conserve))) {                                                                         # REFINEMENT: resort by "conservation rate" of an aminoacid          
      foreach my $symbol (sort(keys(%{$conserve{$column}}))) {                                                                   # for each column look at all found aminoacids                       
         $rate = $conserve{$column}->{$symbol};                                                                                  # get the rate (in % already!)                                       
         $rated{$rate}->{$column} = $symbol;                                                                                     # per rate, save the column number and it's aminoacid                
      }
   }
   foreach my $column (sort(keys(%conserve))) {                                                                                  # REFINEMENT: recalc for percentages                                 
      foreach my $symbol (sort(keys(%{$conserve{$column}}))) {
         $rate = $conserve{$column}->{$symbol};                                                                                  # get the rate (in % already!)                                       
         $cont{$column}->{$rate} = $symbol;
      }
   }
   &_print_rated_motif_results(\%rated, $param);                                                                                 # print the rated results                                            
   return (\%conserve, \%rated, \%cont);
}

sub search_motifs_including_gaps {
   my $clean                     =   $_[0];
   my @cleancodes                = @{$_[1]};
   my $param                     =   $_[2];
   my $len                       = length($clean->{$cleancodes[0]});
   my %conserve                  = ();
   my %rated                     = ();
   my %cont                      = ();
   my ($max, $rate)              = (0, 0);
   my ($symbol, $not_gap_only)   = ("", 0);

   $max = scalar(@cleancodes);                                                                                                   # maximum: number of sequences == number of times an AA may occur per column
   for (my $i = 0; $i < $len; $i++) {                                                                                            # step thru length of alignment
      $not_gap_only = 0;
      foreach my $code (@cleancodes) {                                                                                           # step thru sequences of alignment
         $symbol = substr($clean->{$code}, $i, 1);
         $conserve{$i}->{$symbol}++;                                                                                             # per column: increase counter per aminoacid
         #if ($conserve{$i}->{$symbol} > $max) { $max = $conserve{$i}->{$symbol}; }                                               # maximum-determination
      }
   }
   &_print_raw_motif_results(\%conserve, $param);                                                                                # print the results of the raw motif search method
   foreach my $column (sort(sort_num keys(%conserve))) {                                                                         # REFINEMENT: convert into percentages
      foreach my $symbol (sort(keys(%{$conserve{$column}}))) {                                                                   # for each column look at all found aminoacids
         $rate = &IFG::NutsnBolts::round((($conserve{$column}->{$symbol} / $max) * 100), 2);                                     # and calculate the percentage for number of sequences it occurs in
         $conserve{$column}->{$symbol} = $rate;                                                                                  # save the percentage back
      }
   }
   foreach my $column (sort(sort_num keys(%conserve))) {                                                                         # REFINEMENT: resort by "conservation rate" of an aminoacid
      foreach my $symbol (sort(keys(%{$conserve{$column}}))) {                                                                   # for each column look at all found aminoacids
         $rate = $conserve{$column}->{$symbol};                                                                                  # get the rate (in % already!)
         $rated{$rate}->{$column} = $symbol;                                                                                     # per rate, save the column number and it's aminoacid
      }
   }
   foreach my $column (sort(keys(%conserve))) {                                                                                  # REFINEMENT: recalc for percentages
      foreach my $symbol (sort(keys(%{$conserve{$column}}))) {
         $rate = $conserve{$column}->{$symbol};                                                                                  # get the rate (in % already!)
         $cont{$column}->{$rate} = $symbol;
      }
   }
   &_print_rated_motif_results(\%rated, $param);                                                                                 # print the rated results
   return (\%conserve, \%rated, \%cont);
}

sub _get_seeds {
   my %rated                  = %{$_[0]};
   my $param                  =   $_[1];
   my @seeds                  = ();

   foreach my $rate (reverse(sort(sort_num keys(%rated)))) {
      if ($rate >= $param->{"ti"}) { push(@seeds, sort(sort_num keys(%{$rated{$rate}}))); }                                      # if the conservation-rate if above the threshold, add all these columns...
   }                                                                                                                             # ...to the 'seeds'-array
   @seeds = sort(sort_num @seeds);                                                                                               # sort 'seeds'-array by ascending order
   return (\@seeds, scalar(@seeds));
}

sub _walk_sequence {
   my $rated                  = $_[0];                # rated conservation-rates, sorted by rate
   my $cont                   = $_[1];                # rated conservation-rates, sorted by column
   my $seed                   = $_[2];                # the seed for which we want to search for neighbouring elements
   my $len                    = $_[3];                # length of the alignment
   my $param                  = $_[4];
   my $break                  = 1;
   my %motif                  = ();

   for (my $i = $seed; $i > -1; $i--) {                                                                                          # walk backwards in sequence
      $break = 1;
      foreach my $rate (sort(sort_num keys(%{$cont->{$i}}))) {                                                                   # step thru all conservation-rates for that column
         if ($rate >= $param->{"te"}) {                                                                                          # if a cons.-rate is >= the motif-expansion-threshold...
            $motif{$i} = {"aa"   => $cont->{$i}->{$rate},                                                                        # ...save this position with the rate and the character in a hash.
                          "rate" => $rate,
                         };
            undef($break);                                                                                                       # if there's a cons.-rate above the threshold, we don't want to leave...
         }                                                                                                                       # ...the loop
      }
      last if $break;                                                                                                            # the loop is only left if there is no "big-enough" cons.-rate...
   }                                                                                                                             # ...for the actual column ($i) (== when '$break' remains defined)
   for (my $i = $seed; $i < $len; $i++) {                                                                                        # walk forward in sequence
      $break = 1;
      foreach my $rate (sort(sort_num keys(%{$cont->{$i}}))) {                                                                   # step thru all conservation-rates for that column
         if ($rate >= $param->{"te"}) {                                                                                          # if a cons.-rate is >= the motif-expansion-threshold...
            $motif{$i} = {"aa"   => $cont->{$i}->{$rate},                                                                        # ...save this position with the rate and the character in a hash.
                          "rate" => $rate,
                         };
            undef($break);                                                                                                       # if there's a cons.-rate above the threshold, we don't want to leave...
         }                                                                                                                       # ...the loop
      }
      last if $break;                                                                                                            # the loop is only left if there is no "big-enough" cons.-rate...
   }                                                                                                                             # ...for the actual column ($i) (== when '$break' remains defined)
   return \%motif;
}

sub _print_motif {
   my $motif                  = $_[0];
   my $param                  = $_[1];

   foreach my $column (sort(sort_num keys(%{$motif}))) { print "      [".$motif->{$column}->{"aa"}."-$column]\trate: ".$motif->{$column}->{"rate"}."\n"; }
}

sub _refine_motifs {
   my $motifs                       = $_[0];                   # array with references to hashes of motifs
   my $param                        = $_[1];
   my $num_motifs                   = scalar(@{$motifs});      # number of motifs in array '$motifs'
   my ($num_i, $num_j, $redundancy) = (0, 0, 0);
   my %redundant_motifs             = ();                      # used for indentification of redundant (= similar) motifs
   my @refined_motifs               = ();                      # used for storing the non-redundant motifs

   for (my $i = 0; $i < $num_motifs; $i++) {                                                                                     # step thru all motifs
      $num_i = scalar(keys(%{$motifs->[$i]}));                                                                                   # get the number of aminoacids of the actual motif
      for (my $j = ($i+1); $j < $num_motifs; $j++) {                                                                             # step thru all motifs starting at '$i+1'
         $num_j = scalar(keys(%{$motifs->[$j]}));                                                                                # get the number of aminoacids of the actual motif
         if ($num_j == $num_i) {                                                                                                 # if the number of aminoacids of two motifs are the same...
            if ($param->{"verbose"} && ($param->{"verbose"} == 2)) { print "   Motifs $i <-> $j\thave same number of aminoacids, checking for positional similarity: "; }
            $redundancy = 0;
            foreach my $column (sort(sort_num keys(%{$motifs->[$i]}))) {                                                         # ...check if the aminoacid-positions of both motifs are the same, too.
               $redundancy += exists($motifs->[$j]->{$column});
            }
            if ($redundancy == $num_i) {                                                                                         # if '$redundancy' equals number of aminoacids in the actual motif...
               if ($param->{"verbose"} && ($param->{"verbose"} == 2)) { print "100% redundancy, deleting motif #$j\n"; }         # ...the motif will be deleted.
               $redundant_motifs{$j} = 1;
            }
            else { if ($param->{"verbose"} && ($param->{"verbose"} == 2)) { print " no positional redundancy, keeping motif #$j\n"; } }
         }
      }
   }
   for (my $i = 0; $i < $num_motifs; $i++) {
      unless (exists($redundant_motifs{$i})) { push(@refined_motifs, $motifs->[$i]); }
   }
   return (\@refined_motifs, scalar(@refined_motifs));
}

sub _build_patterned_motifs {
   my $motif                  = $_[0];                # refined/combined motifs
   my $cont                   = $_[1];                # rated conservation-rates, sorted by column
   my $param                  = $_[2];
   my %patterned_motif        = ();

   foreach my $column (sort(sort_num keys(%{$motif}))) {                                                                         # look at each column of a motif for other aminoacids
      if ($param->{"verbose"}) { print "      $column\t- "; }
      foreach my $rate (reverse(sort(sort_num keys(%{$cont->{$column}})))) {                                                     # get *all* aminoacids with *all* rates
         if ($rate >= $param->{"ta"}) {                                                                                          # only include alternative AAs with a cons.-rate >= the "ta"-threshold
            if ($param->{"verbose"}) { print $cont->{$column}->{$rate}."($rate) | "; }
            $patterned_motif{$column}->{$rate} = $cont->{$column}->{$rate};                                                      # save the aminoacid-character with the columns and rates
         }
      }
      if ($param->{"verbose"}) { print "\n"; }
   }
   return \%patterned_motif;
}

sub eval_motifs {
# nothing in here right now
   my $rated                        = $_[0];                # rated conservation-rates, sorted by rate
   my $cont                         = $_[1];                # rated conservation-rates, sorted by column
   my $param                        = $_[2];
   my $len                          = $_[3];                # length of alignment
   my @patterned_motifs             = ();
   my ($motif, $patterned_motif)    = (0, 0);
   my @motifs                       = ();
   my @motif_range                  = ();
   my $motif_score                  = 0;

   my ($seeds, $numseeds) = &_get_seeds($rated, $param);                                                                         # search for seeds in the conservation-rates
   for (my $i = 0; $i < $numseeds; $i++) {
      if ($param->{"verbose"} && ($param->{"verbose"} == 2)) { print "   Examining seed at position ".$seeds->[$i]." for motif #$i:\n"; }
      $motif = &_walk_sequence($rated, $cont, $seeds->[$i], $len, $param);                                                       # "walk" back and forth the sequence beginning from the actual seed
      push(@motifs, $motif);                                                                                                     # save the found motif for the actual seed
      if ($param->{"verbose"} && ($param->{"verbose"} == 2)) { &_print_motif($motif, $param); }                                  # pretty printing of unrefined motifs
   }
   if ($param->{"verbose"} && ($param->{"verbose"} == 2)) { print "   Refining motifs:\n"; }
   my ($refined_motifs, $num_refined) = &_refine_motifs(\@motifs, $param);                                                       # refine/combine found motifs
   if ($param->{"verbose"} && ($param->{"verbose"} == 2)) {                                                                      # pretty printing of refined motifs
      for (my $i = 0; $i < scalar(@{$refined_motifs}); $i++) {
         print "   Refined motif #$i:\n";
         &_print_motif($refined_motifs->[$i], $param);
      }
   }
   for (my $i = 0; $i < $num_refined; $i++) {                                                                                    # "patterning" motifs:
      if ($param->{"verbose"}) { print "   Building patterned motif for motif #$i:\n"; }                                         # look at each column of the motif for alternative AAs
      $patterned_motif = &_build_patterned_motifs($refined_motifs->[$i], $cont, $param);
      push(@patterned_motifs, $patterned_motif);
   }
   for (my $i = 0; $i < $num_refined; $i++) {                                                                                    # pretty printing of patterned motifs
      if ($param->{"verbose"}) { print "   Printing final motif #$i:\n      "; }
      undef(@motif_range);
      $motif_score = 0;
      foreach my $column (sort(sort_num keys(%{$patterned_motifs[$i]}))) {
         unless ($motif_range[0]) {                                                                                              # save the beginning of the motif in the alignment
            $motif_range[0] = $column;                                                                                           # the starting position is saved two times, preventing uninitialized...
            $motif_range[1] = $column;                                                                                           # array-indices for motifs consisting of only one position.
         }
         else { $motif_range[1] = $column; }                                                                                     # successively save the end-position of the motif
         if (scalar(keys(%{$patterned_motifs[$i]->{$column}})) > 1) {
            if ($param->{"verbose"}) { print "["; }
            foreach my $rate (reverse(sort(sort_num keys(%{$patterned_motifs[$i]->{$column}})))) {
               $motif_score += $rate;
               if ($param->{"verbose"}) { print $patterned_motifs[$i]->{$column}->{$rate}; }
            }
            if ($param->{"verbose"}) { print "]"; }
         }
         else {
            foreach my $rate (reverse(sort(sort_num keys(%{$patterned_motifs[$i]->{$column}})))) {
               $motif_score += $rate;
               if ($param->{"verbose"}) { print $patterned_motifs[$i]->{$column}->{$rate}; }
            }
         }
      }
      $motif_score = $motif_score / (($motif_range[1] - $motif_range[0]) + 1);                                                   # divide the score by the length of the motif
      print "\n      [".$motif_range[0]."..".$motif_range[1]."]\n      Score: $motif_score\n";
   }
   return (\@patterned_motifs, $num_refined);
}

sub write_motifs {
   my %rated                  = %{$_[0]};    # rated conservation-rates, sorted by rate
   my @patterned_motifs       = @{$_[1]};
   my $param                  =   $_[2];
   my $paramline              =   $_[3];
   my $len                    =   $_[4];     # length of alignment
   my $num_refined            =   $_[5];
   my $motifseq               = "";
   my @motif_range            = ();
   my $motif_score            = 0;

   open(_motif_, ">$paramline.motifs") or die "[IOErr] Could not open motif_$paramline.txt for write access!\n";
   print _motif_ "Rate\tSequence\n";
   print _motif_ "----\t"."-" x $len."\n";
   foreach my $rate (reverse(sort(sort_num keys(%rated)))) {
      $motifseq = "." x $len;                                                                                                    # prepare the motifsequence-line, which will be printed out
      print _motif_ "$rate\t";
      foreach my $column (keys(%{$rated{$rate}})) {
         substr($motifseq, $column, 1, $rated{$rate}->{$column});
      }
      print _motif_ "$motifseq\n";
   }
   print _motif_ "\n---\n\n";
   print _motif_ "Final motifs:\n-------------\n";
   print _motif_ "Seed-threshold          : ".$param->{"ti"}."%\n";
   print _motif_ "Seed-sprouting-threshold: ".$param->{"te"}."%\n";
   print _motif_ "ALT_INCLUSION-threshold : ".$param->{"ta"}."%\n\n";
   for (my $i = 0; $i < $num_refined; $i++) {                                                                                    # printing of clean motifs
      print _motif_ "Final motif #$i:\n   ";
      undef(@motif_range);
      $motif_score = 0;
      foreach my $column (sort(sort_num keys(%{$patterned_motifs[$i]}))) {                                                       # step thru all columns of the actual motif
         unless ($motif_range[0]) {                                                                                              # save the beginning of the motif in the alignment
            $motif_range[0] = $column;                                                                                           # the starting position is saved two times, preventing uninitialized...
            $motif_range[1] = $column;                                                                                           # array-indices for motifs consisting of only one position.
         }
         else { $motif_range[1] = $column; }                                                                                     # successively save the end-position of the motif
         if (scalar(keys(%{$patterned_motifs[$i]->{$column}})) > 1) {                                                            # if there are more than one AA in this column...
            print _motif_ "[";                                                                                                   # ...put these alternative AAs in square brackets...
            foreach my $rate (reverse(sort(sort_num keys(%{$patterned_motifs[$i]->{$column}})))) {
               $motif_score += $rate;
               print _motif_ $patterned_motifs[$i]->{$column}->{$rate};
            }
            print _motif_ "]";
         }
         else {                                                                                                                  # ...else print out a single AA without any brackets
            foreach my $rate (reverse(sort(sort_num keys(%{$patterned_motifs[$i]->{$column}})))) {
               $motif_score += $rate;
               print _motif_ $patterned_motifs[$i]->{$column}->{$rate};
            }
         }
      }
      $motif_score = $motif_score / (($motif_range[1] - $motif_range[0]) + 1);                                                   # divide the score by the length of the motif
      print _motif_ "\n   [".$motif_range[0]."..".$motif_range[1]."]\n   Score: $motif_score\n";
   }
   close(_motif_);
}

sub keep_sequences_with_motifs_old {
# this subroutine sorts out those sequences in which we can't find all of the (previously) defined motifs
   my $clean                  =   $_[0];
   my @cleancodes             = @{$_[1]};
   my $blast_par              =   $_[2];
   my @compiled_motifs        = @{$_[3]};
   my $param                  =   $_[4];
   #my @motifs                 = @{&get_motifs()};
   my $score                  = 0;
   my (%temp, %del)           = ();
   my (@tempcodes, @delcodes) = ();
   my ($dummy, $delcount)     = (0, 0);
   
   foreach my $code (@cleancodes) {
      $score = 0;
      $dummy = 0;
      for (my $i = 0; $i < scalar(@compiled_motifs); $i++) {
         $compiled_motifs[$i] =~ m/(\d+)\.\.(\d+)/;                                                                              # extract the start and end of a motif ('x..y')
         if ($param->{"exact_motif"}) {                                                                                          # match motif exact? (means: each position within a motif mustn't be a gap-symbol!!!)
            for (my $i = ($1 - 1); $i < $2; $i++) {                                                                              # step thru the range which the actual motif covers ('$1'..'$2')
               $score += (substr($clean->{$code}, $i, 1) ne "-");                                                                # check: is a gap-symbol at the actual position?
               $dummy++;                                                                                                         # we check against the *number of aminoacids in the motifs*
            }
         }
         else {                                                                                                                  # "fuzzy" motif matching: just look if the subject sequence HSP starts before and ends after that motif
            $score += (($1 >= $blast_par->{$code}->{"query_start"}) && 
                       ($2 <= $blast_par->{$code}->{"query_end"}));                                                              # test: does motif[x] lie within the matched part of the query?
            $dummy = scalar(@compiled_motifs);                                                                                   # we check against the number of *motifs*
         }
      }
      if (($score / $dummy) != 1) {                                                                                              # if $score / number_of_motifs = 1, then it means we found all motifs in the query sequence matched
         if ($param->{"verbose"}) { print "   [motif] '$code'\tdoes not match motif requirements and will be deleted\n"; }
         $del{$code} = 1;                                                                                                        # Mark that sequence '$code' has to be deleted
         push(@delcodes, $code);                                                                                                 # dito
         $delcount++;
      }
   }
   if ($param->{"verbose"}) { print "   [motif] Number of deleted sequences: $delcount\n"; }
   foreach my $code (@cleancodes) {                                                                                              # copy sequences to a new hash and leave out those sequences which were
      unless (exists($del{$code})) {                                                                                             # marked as 'to-be-deleted' before
         $temp{$code} = $clean->{$code};
         push(@tempcodes, $code);
      }
   }
   return (\%temp, \@tempcodes);
}

sub keep_sequences_with_motifs_for_conserved {
# this subroutine sorts out those sequences in which we can't find all of the (previously) defined motifs
   my $clean                  =   $_[0];
   my @cleancodes             = @{$_[1]};
   my $blast_par              =   $_[2];
   my @compiled_motifs        = @{$_[3]};
   my $param                  =   $_[4];
   my $score                  = 0;
   my (%temp, %del)           = ();
   my (@tempcodes, @delcodes) = ();
   my ($delete, $delcount)    = (0, 0);
   my @motiflen               = ();
   my @score                  = ();
   my $d                      = 0;
   
   foreach my $code (@cleancodes) {
      undef($delete);
      undef(@score);
      undef(@motiflen);
      if ($d) { print ">>> $code:\t"; }
      for (my $i = 0; $i < scalar(@compiled_motifs); $i++) {
         $compiled_motifs[$i] =~ m/(\d+)\.\.(\d+)/;                                                                              # extract the start and end of a motif ('x..y')
         for (my $j = $1; $j < ($2 + 1); $j++) {                                                                                 # step thru the range which the actual motif covers ('$1'..'$2')
            $score[$i] += (substr($clean->{$code}, $j, 1) eq "-");                                                               # check: is a gap-symbol at the actual position?
            $motiflen[$i]++;                                                                                                     # we check against the *number of aminoacids in the motifs*
         }
         if ($d) { print "[gaps ".$score[$i]." / len ".$motiflen[$i]."]  -  "; }
         if ($score[$i] > 0) { $delete = 1; }                                                                                    # if there is a gap inside the motif range, delete the sequence
      }
      if ($d) { print "\n"; }
      if ($delete) {
         if ($param->{"verbose"}) { print "   [motif] '$code'\tdoes not match motif requirements and will be deleted\n"; }
         $del{$code} = 1;                                                                                                        # sequence '$code' has to be deleted
         $delcount++;
      }
   }
   if ($param->{"verbose"}) { print "   [motif] Number of deleted sequences: $delcount\n"; }
   foreach my $code (@cleancodes) {                                                                                              # copy sequences to a new hash and leave out those sequences which were
      unless (exists($del{$code})) {                                                                                             # marked as 'to-be-deleted' before
         $temp{$code} = $clean->{$code};
         push(@tempcodes, $code);
      }
   }
   return (\%temp, \@tempcodes);
}

sub evaluate_conservation_rates_1 {
   my $conserve                           = $_[0];
   my $rated                              = $_[1];          # "rated" is the way to go! => rated->{rate}->{column}->{symbol}
   my $cont                               = $_[2];
   my $param                              = $_[3];
   my ($mean, $num)                       = (0, 0);
   my %conserved_positions                = ();

   open(_file_, ">conssrates.txt") or die "[IOErr] Could not open 'conssrates.txt' for write access!\n";
   foreach my $rate (reverse(sort(sort_num keys(%{$rated})))) {                                                                  # '$key' holds the rate!
      if ($rate >= 80) {                                                                                                         # if the conservation rate is >= 80%...
         foreach my $column (sort(sort_num keys(%{$rated->{$rate}}))) {                                                          # ...get all the columns which are conserved by >= 80%...
            $conserved_positions{$column}->{'rate'}   = $rate;
            $conserved_positions{$column}->{'symbol'} = $rated->{$rate}->{$column};
         }
      }
      $mean += $rate;
      $num++;
      print _file_ $rate."\n";
   }
   close(_file_);
   $mean = $mean / $num;
   print "   Mean conservation rate: $mean\n";
   return \%conserved_positions;
}

sub evaluate_conservation_rates {
   my $conserve                           = $_[0];
   my $rated                              = $_[1];          # "rated" is the way to go! => rated->{rate}->{column}->{symbol}
   my $cont                               = $_[2];
   my $param                              = $_[3];
   my ($max, $mean, $num)                 = (0, 0, 0);
   my %conserved_positions                = ();
   my @dummy                              = ();

   foreach my $column (sort(sort_num keys(%{$cont}))) {
      ($max, $mean, $num) = (0, 0, 0);
      @dummy = ();
      foreach my $rate (reverse(sort(keys(%{$cont->{$column}})))) {
         $num++;
         $mean += $rate;
         if ($rate > $max) { $max = $rate; }
         push(@dummy, $cont->{$column}->{$rate}.', ');
      }
      $mean = ($mean / $num) * $max;                                                                                                # eval each column
      if ($mean >= 750) {
         print "[$column]\t($mean / $max)\t@dummy\n";
         $conserved_positions{$column}->{'rate'}   = $max;
         $conserved_positions{$column}->{'symbol'} = $cont->{$column}->{$max};
      }
   }
   return \%conserved_positions;
}

sub keep_sequences_with_motifs {
   print ">>>>>> VisCoSe::Motif::keep_sequences_with_motifs() was deprecated and moved to RiPE::Pratt::keep_sequences_with_motifs()!!! Exiting...\n";
   exit 0;
}

sub manual_motifs {
   my @motifs_abcg = ("117..122",      # Walker A     (ABCG1 as query)
                      "213..217",      # Signature    ("")
                      "236..244",      # Walker B     ("")
                     );
   return @motifs_abcg;
}

sub print_conserved_positions {
   my $conserved_positions    = $_[0];
   my $cont                   = $_[1];
   my $param                  = $_[2];
   my @conserved_ranges       = ();
   my ($begin, $end)          = ('', '');

   foreach my $column (sort(sort_num keys(%{$conserved_positions}))) {
      print "   [$column]\t".$conserved_positions->{$column}->{'rate'}."\t- >".$conserved_positions->{$column}->{'symbol'}."<\n";
      foreach my $rate (reverse(sort(sort_num keys(%{$cont->{$column+1}})))) {                                                      # print (for information only!) the first NEXT following position
         print "   (".($column+1).")\t$rate\t-  ".$cont->{$column+1}->{$rate}."\t";
         last;
      }
      if ($begin && $end) {
         if ($column != ($end + 1)) {                                                                                               # consecutive positions? If not, compile a range "x..y" and start over with actual markers
            push(@conserved_ranges, $begin.'..'.$end);
            $begin = $column;
            $end   = $column;
            print "|".$conserved_ranges[-1]."|\n";
         }
         else {
            $end = $column;                                                                                                         # if we have a consecutive number, just adjust the end-marker
            print "\n";
         }
      }
      else {                                                                                                                        # initially set markers (only executed in the first roll of the loop)
         $begin = $column;
         $end   = $column;
         print "\n";
      }
   }
   push(@conserved_ranges, $begin.'..'.$end);                                                                                       # add the last range using the last actual markers
   print "\t\t\t\t|".$conserved_ranges[-1]."|\n";
   return \@conserved_ranges;
}

sub _calc_mean_conservation_rate {
   my $conserve         = $_[0];
   my $cont             = $_[1];
   my $param            = $_[2];
   my %consensus        = ();

   foreach my $column (sort(sort_num keys(%{$conserve}))) {                                                                         # step thru each column of the conserve-hash
      my $num_symbols  = scalar(keys(%{$conserve->{$column}}));
      my $sum_consrate = 0;
      my @rates = reverse(sort(sort_num keys(%{$cont->{$column}})));                                                                # $rates[0] will hold the *highest* rate
      foreach my $symbol (keys(%{$conserve->{$column}})) {                                                                          # for each column, step thru each symbol existing in that particular column
         if ($symbol ne '-') {
            $sum_consrate += $conserve->{$column}->{$symbol};
         }
      }
      my $mean = $sum_consrate / $num_symbols;
      $consensus{$column} = {'symbol'  => $cont->{$column}->{$rates[0]},                                                            # save rate (= noise) and symbol per column
                             'rate'    => $mean,
                            };
   }
   return \%consensus;
}

sub calc_mean_conservation_rate_per_file {
   my $conserve         = $_[0];
   my $cont             = $_[1];
   my $param            = $_[2];
   my %consensus        = ();

   foreach my $file (keys(%{$cont})) {
      if ($param->{'verbose'}) { print "   $file..."; }
      ($consensus{$file}) = &_calc_mean_conservation_rate($conserve->{$file}, $cont->{$file}, $param);                              # extracts the consensus sequence from conservation rates per column/AA
      if ($param->{'verbose'}) { print "done\n"; }
   }
   return \%consensus;
}

sub _calc_noise {
   my $conserve         = $_[0];
   my $cont             = $_[1];
   my $param            = $_[2];
   my %consensus        = ();

   foreach my $column (sort(sort_num keys(%{$conserve}))) {                                                                         # step thru each column of the conserve-hash
      my $num_symbols  = scalar(keys(%{$conserve->{$column}}));
      my $sum_consrate = 0;
      my @rates = reverse(sort(sort_num keys(%{$cont->{$column}})));                                                                # $rates[0] will hold the *highest* rate
      foreach my $symbol (keys(%{$conserve->{$column}})) {                                                                          # for each column, step thru each symbol existing in that particular column
         if ($symbol ne '-') {
            $sum_consrate += $conserve->{$column}->{$symbol};
         }
      }
      my $noise = ($sum_consrate * $num_symbols) / 20;
      $consensus{$column} = {'symbol'  => $cont->{$column}->{$rates[0]},                                                            # save rate (= noise) and symbol per column
                             'rate'    => $noise,
                            };
   }
   return \%consensus;
}

sub calc_noise_per_file {
   my $conserve         = $_[0];
   my $cont             = $_[1];
   my $param            = $_[2];
   my %consensus        = ();

   foreach my $file (keys(%{$cont})) {
      if ($param->{'verbose'}) { print "   $file..."; }
      ($consensus{$file}) = &_calc_noise($conserve->{$file}, $cont->{$file}, $param);                                               # extracts the consensus sequence from conservation rates per column/AA
      if ($param->{'verbose'}) { print "done\n"; }
   }
   return \%consensus;
}
