###################################################################################################
# ConsensusTools.pm
# (c) 2002 Michael Spitzer, IFG Muenster (Integrated Funtional Genomics), Germany
#
# 18.03.2003   added subroutines

package VisCoSe::ConsensusTools;
require Exporter;

use strict;
use diagnostics;

our @EXPORT  = qw(get_consensus                       get_mean_gap_column_conservation_rate     get_consensus_including_gaps        prepare_sequence_consensus
                  write_html_consensus                read_alignments                           get_AA_frequency_for_alignments     get_consensus_from_alignments
                  format_consensus_for_alignments     generate_HTML_consensus_filename          write_consensus_for_alignments      negate_consensus);

our $VERSION = "18.03.2003";



sub sort_num { $a <=> $b }

sub get_consensus {
   my $cont             = $_[0];
   my $param            = $_[1];
   my @rates            = ();
   my %consensus        = ();

   foreach my $column (sort(sort_num keys(%{$cont}))) {                             # step thru all columns...
      @rates = reverse(sort(sort_num keys(%{$cont->{$column}})));                   # $rates[0] will hold the *highest* rate
      $consensus{$column} = {'symbol'  => $cont->{$column}->{$rates[0]},            # save rate and symbol per column
                             'rate'    => $rates[0],
                            };
   }
   return \%consensus;
}

sub get_mean_gap_conservation_rate {
   my $cont                         = $_[0];
   my @rates                        = ();
   my ($mean_gap_conss, $gap_num)   = (0, 0);

   foreach my $column (sort(sort_num keys(%{$cont}))) {                                # step thru all columns...
      foreach my $rate (keys(%{$cont->{$column}})) {
         if ($cont->{$column}->{$rate} eq '-') {
            $mean_gap_conss += $rate;
            $gap_num++;
         }
      }
   }
   if ($gap_num == 0) { return (0, 0, 0); }                                            # no gaps? => return all zero!
   else { return (($mean_gap_conss / $gap_num), $mean_gap_conss, $gap_num); }          # to get the mean conservation rate for gaps here: divide the summed conservation rate by the number of columns
}

sub get_mean_gap_column_conservation_rate {
   my $cont                         = $_[0];
   my @rates                        = ();
   my ($mean_gap_conss, $gap_num)   = 0;

   foreach my $column (sort(sort_num keys(%{$cont}))) {                                # step thru all columns...
      @rates = reverse(sort(sort_num keys(%{$cont->{$column}})));                      # $rates[0] will hold the *highest* rate
      if ($cont->{$column}->{$rates[0]} eq '-') {                                      # if $rates[0] is a gap, measure it's conservation rate
         $mean_gap_conss += $rates[0];                                                 # add the conservation rate
         $gap_num++;                                                                   # gap_num counts the number of columns with gaps as the highest conserved character
      }
   }
   return (($mean_gap_conss / $gap_num), $mean_gap_conss, $gap_num);                   # to get the mean conservation rate for gaps here: divide the summed conservation rate by the number of columns
}

sub get_consensus_including_gaps {
   my $cont                         = $_[0];
   my $threshold4gaps               = $_[1];
   my $param                        = $_[2];
   my @rates                        = ();
   my %consensus                    = ();

   foreach my $column (sort(sort_num keys(%{$cont}))) {                                # step thru all columns...
      @rates = reverse(sort(sort_num keys(%{$cont->{$column}})));                      # $rates[0] will hold the *highest* rate
      if ($cont->{$column}->{$rates[0]} ne '-') {                                      # if the first entry (= highest conserved character) in @rates is *NOT* a gap, get the amino acid for the consensus sequence
         $consensus{$column} = {'symbol'  => $cont->{$column}->{$rates[0]},
                                'rate'    => $rates[0],
                               };
      }
      elsif ($rates[0] > $threshold4gaps) {                                            # if $rates[0] is a gap *AND* it's above the threshold, use it for the consensus
         $consensus{$column} = {'symbol'  => $cont->{$column}->{$rates[0]},
                                'rate'    => $rates[0],
                               };
      }
      else {                                                                           # else ($rates[0] is a gap but below threshold): use next highest conserved character for consensus
         $consensus{$column} = {'symbol'  => $cont->{$column}->{$rates[1]},
                                'rate'    => $rates[1],
                               };
      }
   }
   return \%consensus;
}



sub prepare_sequence_consensus {
   my $consensus                       = $_[0];
   my $filename                        = $_[1];
   my $param                           = $_[2];
   my (%conss, @consscode)             = ((), ());
   my $seq                             = '';

   foreach my $column (sort(sort_num keys(%{$consensus}))) { $seq .= $consensus->{$column}->{'symbol'}; }
   $conss{$filename} = $seq;
   push(@consscode, $filename);
   return (\%conss, \@consscode);
}



sub write_html_consensus {
   my $file                = $_[0];
   my $filename            = $_[1];
   my $param               = $_[2];
   my ($header, $footer)   = &VisCoSe::HTMLstuff::prep_header($filename, $param);

   open(_file_, ">$filename") or die "[IOErr] Could not open 'conss_".$param->{'in'}.".html' for write access!\n";
   for (my $i = 0; $i < scalar(@{$header}); $i++) { print _file_ $header->[$i],"\n"; }
   for (my $i = 0; $i < scalar(@{$file}); $i++)   { print _file_ $file->[$i],"\n";   }
   for (my $i = 0; $i < scalar(@{$footer}); $i++) { print _file_ $footer->[$i],"\n"; }
   close(_file_);
}



sub read_alignments {
   my $param            = $_[0];
   my (%in, %incodes)  = ((), (), (), ());              # we'll store alignments here ($in->{filename}->{deflines}->{sequence_data})

   foreach my $file (@{$param->{'in'}}) {
      if ($param->{'verbose'}) { print "   $file..."; }
      ($in{$file}, $incodes{$file}) = &IFG::Alignments::read_aln($file);
      $in{$file} = &IFG::Alignments::uppercase_seqs($in{$file}, $incodes{$file});
      &IFG::Alignments::write_aln($file, $in{$file}, $incodes{$file});                       # rewrite the alignment back to disk in case there were doubled deflines
      if ($param->{'verbose'}) { print "done\n"; }
   }
   return (\%in, \%incodes);
}



sub get_AA_frequency_for_alignments {
   my $in                        = $_[0];
   my $incodes                   = $_[1];
   my $param                     = $_[2];
   my (%conserve, %rated, %cont) = ((), (), ());

   foreach my $file (keys(%{$in})) {
      if ($param->{'verbose'}) { print "   $file..."; }
      ($conserve{$file}, $rated{$file}, $cont{$file}) = &VisCoSe::Motif::search_motifs($in->{$file}, $incodes->{$file}, $param);    # search for motifs / aminoacid frequencies per alignment-position (PER FILE!!!)
      if ($param->{'verbose'}) { print "done\n"; }
   }
   return (\%conserve, \%rated, \%cont);
}



sub get_consensus_from_alignments {
   my $cont                      = $_[0];
   my $param                     = $_[1];
   my %consensus                 = ();

   foreach my $file (keys(%{$cont})) {
      if ($param->{'verbose'}) { print "   $file..."; }
      ($consensus{$file}) = &get_consensus($cont->{$file}, $param);                 # extracts the consensus sequence from conservation rates per column/AA
      if ($param->{'verbose'}) { print "done\n"; }
   }
   return \%consensus;
}



sub format_consensus_for_alignments {
   my $in                        = $_[0];
   my $incodes                   = $_[1];
   my $consensus                 = $_[2];
   my $conserve                  = $_[3];
   my $prevconss                 = $_[4]; # eventually holds previously computed and formatted consensus sequences (using the default aa alphabet)
   my $alphabet_name             = $_[5]; # if $prevconss is defined, $alphabet_name has to be defined to!
   my $param                     = $_[6];
   my (%htmlfile, %seq)          = ();
   my (%conss, %consscode)       = ((), ());

   foreach my $file (keys(%{$consensus})) {
      if ($param->{'verbose'}) { print "   $file..."; }
      ($conss{$file}, $consscode{$file}) = &prepare_sequence_consensus($consensus->{$file}, $file, $param);
      if ($prevconss) {
         ($htmlfile{$file}, $seq{$file}) = &VisCoSe::HTMLstuff::prepare_single_html_consensus($consensus->{$file}, $in->{$file}, $incodes->{$file}, $file,
                                                                                              $conserve->{$file}, $prevconss->{$file}, $alphabet_name, $param);
      }
      else {
         ($htmlfile{$file}, $seq{$file}) = &VisCoSe::HTMLstuff::prepare_single_html_consensus($consensus->{$file}, $in->{$file}, $incodes->{$file}, $file,
                                                                                              $conserve->{$file}, '', '', $param);
      }
      if ($param->{'linebreak'}) {
         my $maxcodelen = &VisCoSe::HTMLstuff::get_max_defline_length($incodes->{$file}, $param);
         $htmlfile{$file} = &VisCoSe::HTMLstuff::split_html_file($htmlfile{$file}, $maxcodelen, $param);
      }
      if ($param->{'verbose'}) { print "done\n"; }
   }
   return (\%conss, \%consscode, \%htmlfile, \%seq);
}

sub generate_HTML_consensus_filename {
   my $file                      = $_[0];
   my $filename                  = "";

   ($filename = $file) =~ s/(?:\.[A-Za-z0-9_]+)\Z//;
   $filename = 'consensus_'.$filename;
   return $filename;
}

sub write_consensus_for_alignments {
   my $conss                     = $_[0];
   my $consscode                 = $_[1];
   my $htmlfile                  = $_[2];
   my $suffix                    = $_[3];
   my $param                     = $_[4];
   my $filename                  = "";

   foreach my $file (keys(%{$conss})) {
      if ($param->{'verbose'}) { print "   $file..."; }
      $filename = &generate_HTML_consensus_filename($file);
      if ($suffix) { $filename .= '_'.$suffix; }
      &write_html_consensus($htmlfile->{$file}, ($filename.'.html'), $param);
      &IFG::Alignments::write_aln(($filename.'.fasta'), $conss->{$file}, $consscode->{$file});
      if ($param->{'verbose'}) { print "done\n"; }
   }
}

sub negate_consensus {
   my $full_conss    = $_[0];                # the gapped consensus (plain string)
   my $gap_conss     = $_[1];                # the full consensus (plain string)
   my $param         = $_[2];
   my $neg_conss     = '';

   for (my $i = 0; $i < length($full_conss); $i++) {                                                                                      # step thru the original consensus char by char
      if ((substr($gap_conss, $i, 1) eq '-') || (substr($gap_conss, $i, 1) eq 'X')) { $neg_conss .= substr($full_conss, $i, 1); }         # if the actual position is a gap append the corresponding amino acid character...
      else { $neg_conss .= 'X'; }                                                                                                         # ...else (the actual position is an amino acid character) append an 'X'
   }
   return $neg_conss;
}
