###################################################################################################
# simplified.pm
# (c) 2003 Michael Spitzer, IFG Muenster (Integrated Funtional Genomics), Germany
#

package VisCoSe::Simplified;
require Exporter;

use strict;
use diagnostics;
use Data::Dumper;

our @EXPORT  = qw(read_alphabets          verify_alphabet_groups_for_ambiguity         duplicate_alignments       translate_alignments);

our $VERSION = "25.07.2003";

sub sort_num { $a <=> $b }

sub read_alphabets {
   my $param                                              = $_[0];
   my ($line, $alphabet_name, $group_name, $group_symbol) = ('', '', '', '');
   my (%alphabets, @group_symbols)                        = ((), ());
   my $counter                                            = 0;

   open(_file_, '<'.$param->{'aa'}) or die "[IOErr] Could not open file with simplified amino acid alphabets '".$param->{'aa'}."'\n";
   until (eof(_file_)) {
      $line = <_file_>;
      chomp($line);
      $line =~ s/\#.*//;                                                                           # remove comments from line if there's any
      if ($line =~ m/\[(.+)\]/) {                                                                  # "[xyz]" ? -> New alphabet!
         $alphabet_name = $1;                                                                      # remember alphabet name
         $counter = 0;
         if ($param->{'debug'}) { print ">>> '$alphabet_name'\n"; }
      }
      elsif ($line =~ m/([A-Za-z]+)(?: |\t)+([A-Za-z])(?: |\t)*(?:\((.+)\))?/) {                       # "alphabet_group alphabet_symbol group_name" ?
         $counter++;
         @group_symbols = split(//, $1);                                                           # split the group symbols and save them in an array
         $group_symbol  = $2;                                                                      # remember the common symbol for the group members
         $group_name    = $3 || $counter;                                                          # remember the group name (if specified, else use $counter)
         if ($param->{'debug'}) { print "    '$1'   -   '$2'   -   '$3'   -   @group_symbols\n"; }
         foreach my $symbol (@group_symbols) {                                                     # save each symbol of a group member with its common...
            $alphabets{$alphabet_name}->{$group_name}->{$symbol} = $group_symbol;                  # ...group symbol in a hash using the alphabet name as key
         }
      }
   }
   close(_file_);
   return \%alphabets;
}

sub _compare_alphabet_groups {
   my $group1           = $_[0];
   my $group2           = $_[1];
   my $param            = $_[2];
   my @ambiguities      = ();

   foreach my $symbol1 (@{$group1}) {
      foreach my $symbol2 (@{$group2}) {
         if ($symbol1 eq $symbol2) { push(@ambiguities, $symbol1); }
      }
   }
   if (@ambiguities) { return \@ambiguities; }
   else { return undef; }
}

sub verify_alphabet_groups_for_ambiguity {
   my $alphabets           = $_[0];
   my $param               = $_[1];
   my @alphabet_groups     = ();
   my $ambiguities         = 0;
   my (@group1, @group2)   = ((), ());

   foreach my $alphabet (keys(%{$alphabets})) {                                                    # step thru all defined alphabets
      @alphabet_groups = keys(%{$alphabets->{$alphabet}});                                         # get groups of actual alphabet
      if ($param->{'debug'}) { print "[$alphabet]\n"; }
      for (my $i = 0; $i < (scalar(@alphabet_groups) - 1); $i++) {                                 # step thru array from 0..n-1
         for (my $j = ($i + 1); $j < scalar(@alphabet_groups); $j++) {                             # step thru array from 1..n
            if ($param->{'debug'}) { print "   ".$alphabet_groups[$i]." vs. ".$alphabet_groups[$j]."\n"; }
            @group1 = keys(%{$alphabets->{$alphabet}->{$alphabet_groups[$i]}});
            @group2 = keys(%{$alphabets->{$alphabet}->{$alphabet_groups[$j]}});
            $ambiguities = _compare_alphabet_groups(\@group1, \@group2, $param);
            if ($ambiguities) { print "Found ambiguities in [$alphabet]:\n".
                                      "   Group ".$alphabet_groups[$i]." vs. ".$alphabet_groups[$j].": @{$ambiguities}\n"; }
         }
      }
   }
}

sub _translate_single_alignment {
   my $aln           = $_[0];
   my $alncodes      = $_[1];
   my $alphabet      = $_[2];
   my $param         = $_[3];

   foreach my $code (@{$alncodes}) {                                                                                 # step thru all sequences
      for (my $i = 0; $i < length($aln->{$code}); $i++) {                                                            # step thru the actual sequence character by character
         my $aa = substr($aln->{$code}, $i, 1);
         foreach my $agroup (keys(%{$alphabet})) {                                                                   # step thru each *group* of the chosen alphabet
            if ($alphabet->{$agroup}->{$aa}) {                                                                       # if the current amino acid at position $i is existant in the current alphabet group "$agroup" (a hash!)...
               substr($aln->{$code}, $i, 1, $alphabet->{$agroup}->{$aa});                                            # ...then translate!
            }
         }
      }
   }
   return $aln;
}

sub translate_alignments {
   my $in            = $_[0];
   my $incodes       = $_[1];
   my $alphabet      = $_[2];
   my $param         = $_[3];

   foreach my $file (keys(%{$incodes})) {                                                    # step thru all the files...
      if ($param->{'verbose'}) { print "   $file..."; }
      $in->{$file} = _translate_single_alignment($in->{$file}, $incodes->{$file}, $alphabet, $param);   # ...and translate each single alignment
      if ($param->{'verbose'}) { print "done\n"; }
   }
   return $in;
}
