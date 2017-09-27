###################################################################################################
# AlignTools.pm
# (c) 2003 Michael Spitzer, IFG Muenster (Integrated Funtional Genomics), Germany
#
# 20.03.2003   added _extract_group_sequences(), split_dataset_into_groups(), write_group_alignments()
# 18.03.2003   initial version

package VisCoSe::AlignTools;
require Exporter;

use strict;
use diagnostics;

our @EXPORT  = qw(test_for_aligned_data      align_unaligned_datasets      generate_consensi_dataset     align_consensus_dataset
                  split_dataset_into_groups     write_group_alignments);

our $VERSION = "20.03.2003";



sub sort_num { $a <=> $b }



sub test_for_aligned_data {
   my $in            = $_[0];
   my $incodes       = $_[1];
   my $param         = $_[2];
   my ($isaln, $len) = (0, 0);

   foreach my $file (keys(%{$in})) {
      if ($param->{'verbose'}) { print "   $file..."; }
      $len   = length($in->{$file}->{$incodes->{$file}->[0]});                                        # get length of the first sequence in the alignment '$file'
      $isaln = 1;
      foreach my $code (@{$incodes->{$file}}) {
         if (length($in->{$file}->{$code}) != $len) { $isaln = 0; }
      }
      if ($isaln && (scalar(@{$incodes->{$file}}) > 1)) {                                             # a dataset with only one sequence is always considered "unaligned"!
         $in->{$file}->{'aligned'} = 'yes';                                                           # store a flag that alignment '$file' is an aligned dataset!
         $param->{'aligned_source'} = 1;                                                              # this is not completely sane...I should work that out in a better way, seperated per dataset...
      }
      else { $in->{$file}->{'aligned'} = 'no'; }
      if ($param->{'verbose'}) {
         if ($isaln) { print "is aligned already (all sequences have the same lengths)\n"; }
         else        { print "unaligned\n"; }
      }
   }
   return ($in, $incodes);
}

sub align_unaligned_datasets {            # 18.03.2003
   my $in                  = $_[0];
   my $incodes             = $_[1];
   my $param               = $_[2];
   my ($log, $shellstring) = ("", "");

   foreach my $file (keys(%{$in})) {
      if ($param->{'verbose'}) { print "   $file..."; }
      if (($in->{$file}->{'aligned'} eq 'no') && (scalar(@{$incodes->{$file}}) > 1)) {                # if we have unaligned data:
         unless (-s $file.'_aligned') {                                                               # if we noticed earlier that a dataset contains no gaps...
            $shellstring = 'nice -19 fftnsi '.$file.' >'.$file.'_aligned';                            # ...we align it using MAFFT/FFTNS (closing STDERR = default MAFFT output)
            $log = `$shellstring`;                                                                    # RUN FORREST, RUN!!!
            if ($param->{'verbose'}) { print "aligning done\n"; }
         }
         else { print "aligned data already exists\n"; }
         ($in->{$file}, $incodes->{$file}) = &IFG::Alignments::read_aln($file.'_aligned');            # read back the aligned dataset for the corresponding filename
      }
      else { if ($param->{'verbose'}) { print "no need to align\n"; } }
   }
   return ($in, $incodes);
}



sub generate_consensi_dataset {              # 18.03.2003
   my $conss                  = $_[0];       # the consensi for all files passed by the user ($conss->{file}->{sequence_name} - sequence_name == filename)
   my $consscode              = $_[1];       # the corresponding codes ($consscode->{file}->[array_of_sequence_names] - array for compatibility reasons)
   my $param                  = $_[2];
   my (%allconss, @allcodes)  = ((), ());

   foreach my $file (keys(%{$conss})) {
      if ($param->{'verbose'}) { print "   $file..."; }
      $allconss{$file} = $conss->{$file}->{$file};                   # consensus sequence data is store under $conss->{filename}->{filename}
      push(@allcodes, $file);                                        # the first 'filename' points to the file from which the consensus was derived
      if ($param->{'verbose'}) { print "done\n"; }                   # the second 'filename' points to the actual consensus sequence
   }
   return (\%allconss, \@allcodes);
}



sub align_consensus_dataset {             # 18.03.2003
   my $allconss            = $_[0];
   my $allcodes            = $_[1];
   my $param               = $_[2];
   my ($log, $shellstring) = ("", "");

   unless (-s 'consensi_aligned.fasta') {                                                    # only align if aligned dataset does not exist already
      &IFG::Alignments::write_aln('consensi_unaligned.fasta', $allconss, $allcodes);         # write the unaligned consensus sequences
      $shellstring = 'nice -19 fftnsi consensi_unaligned.fasta >consensi_aligned.fasta';
      $log = `$shellstring`;
   }
   else { if ($param->{'verbose'}) { print "  consensi_aligned.fasta already exists\n"; } }
   ($allconss, $allcodes) = &IFG::Alignments::read_aln('consensi_aligned.fasta');            # read back in the aligned consensus sequences
   return ($allconss, $allcodes);
}



sub _extract_group_sequences {
   my $in                     = $_[0];          # reference to a *SINGLE* alignment
   my $incodes                = $_[1];          # reference to a *SINGLE* defline-array
   my $seqarray               = $_[2];          # reference to the array containing the sequence numbers of the actual group from 'split_dataset_into_groups()'
   my $param                  = $_[3];
   my (%group, @groupcodes)   = ((), ());

   foreach my $seqnum (@{$seqarray}) {
      if ($param->{'verbose'} && ($param->{'verbose'} == 2)) { print "$seqnum: ".$incodes->[$seqnum-1]."\n"; }
      $group{$incodes->[$seqnum-1]} = $in->{$incodes->[$seqnum-1]};                       # 'seqnum-1' is necessary to get from human-wise counting (beginning with '1')...
      push(@groupcodes, $incodes->[$seqnum-1]);                                           # ...to perlwise-counting (beginning with '0')
   }
   return (\%group, \@groupcodes);                                                        # return the 'alignment' of the grouped sequences
}



sub split_dataset_into_groups {
   my $in               = $_[0];
   my $incodes          = $_[1];
   my $param            = $_[2];
   my @groups           = split(/\;/, $param->{'groups'});     # split the 'groups'-string at ";"
   my $seqvector        = 0;
   my @seqarray         = ();

   foreach my $file (keys(%{$in})) {
      $seqvector = Bit::Vector->new(scalar(@{$incodes->{$file}}) + 100);                  # create new Bitvector with length (number of sequences + 100 security range ;-)
      if ($param->{'fra'}) {                                                              # do not degap if '-fra' was not defined on commandline
         ($in->{$file}, $incodes->{$file}) = &IFG::Alignments::degap_alignment($in->{$file}, $incodes->{$file});  # degap the alignment we want to extract the groups from
      }
      foreach my $group (@groups) {
         if ($param->{'verbose'}) { print "   group '$group'..."; }
         $seqvector->from_Enum($group);                                                   # generating Bitvector from user-defined range
         @seqarray = $seqvector->Index_List_Read();                                       # converting Bitvector into array with fully interpreted ranges
         ($in->{'grp_'.$group.'_'.$file}, $incodes->{'grp_'.$group.'_'.$file}) = &_extract_group_sequences($in->{$file}, $incodes->{$file}, \@seqarray, $param);
         unless ($param->{'fra'}) { $in->{'grp_'.$group.'_'.$file}->{'aligned'} = 1; }
         if ($param->{'verbose'}) { print "done\n"; }
      }
      delete($in->{$file});                                                               # delete the data of the original file from the hash with alignment data
      delete($incodes->{$file});                                                          # delete the data of the original file from the hash with defline data
   }
   $param->{'numfiles'} = scalar(@groups);
   return ($in, $incodes, $param);
}



sub write_group_alignments {
   my $in               = $_[0];
   my $incodes          = $_[1];
   my $param            = $_[2];

   foreach my $file (keys(%{$in})) {
      unless (-s $file) {
         &IFG::Alignments::write_aln($file, $in->{$file}, $incodes->{$file});             # write the unaligned consensus sequences
      }
   }
}

