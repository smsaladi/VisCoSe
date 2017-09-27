#! /usr/bin/perl -w
###################################################################################################
# conss_aln.pl
# (c) 2003 Michael Spitzer, IFG-IZKF
#
# 17.02.2004   fixed a bug with range checking in "test_group_params()"
# 18.03.2003   - added checking of input dataset(s) for gaps (aligned vs. unaligned)
#              - realigns input dataset(s) if necessary using MAFFT/FFTNSI
#              - added routines for generation and alignment of multiple consensus seqs
# 17.03.2003   preparing for handling of multiple alignment *OR* group-defining with an alignment
# 14.03.2003   updated version, split into single modules, see seperate history for each module

#require v5.6.1;
use diagnostics;                                            # verbose warning- and error-messages.
use strict;                                                 # This is a "must-have".

use lib '/var/www/vhosts/viscose/cgi-bin';                  # add path where my custom Perl-Modules reside (~/bin/perl/)

use IFG::Alignments;                                        # collected useful stuff done (mostly) by me, myself, and I ;-)
use IFG::NutsnBolts;                                        # useful stuff again
use VisCoSe::Motif;                                         # subroutines dealing with motif-finding stuff from RiPE
use VisCoSe::HTMLstuff;                                     # subroutines for preparing the HTML-consensus pages
use VisCoSe::ConsensusTools;                                # subroutines for Consensus-IO and -preparation
use VisCoSe::AlignTools;                                    # subroutines for multiple Alignment-Handling and -Preparation
use VisCoSe::Simplified;                                    # subroutines for Consensus-IO and -preparation
use CGI::Carp qw(fatalsToBrowser);
use Data::Dumper;
if (eval('require Bit::Vector::Overload')) { require Bit::Vector::Overload; }
else { require Bit::Vector; }

sub sort_num { $a <=> $b }

sub parse_args {
# parses the arguments given via command line and builds a hash which will be returned.
# Does a little error checking and dies if something's wrong with the parameters (eg. some switch is unknown or so).
   my @args          = @{$_[0]};        # holds @ARGV
   my %param         = ();              # will be filled with the parameters from @ARGV

   for (my $i = 0; $i < scalar @args; $i++) {
      CASE: { 
         if ($args[$i] eq '-in')             { $param{'in'}             = $args[++$i]; last CASE; }
         if ($args[$i] eq '-groups')         { $param{'groups'}         = $args[++$i]; last CASE; }
         if ($args[$i] eq '-fra')            { $param{'fra'}            = 1;           last CASE; }
         if ($args[$i] eq '-gh')             { $param{'graph_height'}   = $args[++$i]; last CASE; }
         if ($args[$i] eq '-cg')             { $param{'colorgraph'}     = 1;           last CASE; }
         if ($args[$i] eq '-sa')             { $param{'sa'}             = 1;           last CASE; }
         if ($args[$i] eq '-ca')             { $param{'ca'}             = 1;           last CASE; }
         if ($args[$i] eq '-fa')             { $param{'fa'}             = 1;           last CASE; }
         if ($args[$i] eq '-v')              { $param{'verbose'}        = 1;           last CASE; }
         if ($args[$i] eq '-vv')             { $param{'verbose'}        = 2;           last CASE; }
         if ($args[$i] eq '-cgi')            { $param{'cgi'}            = 1;           last CASE; }
         if ($args[$i] eq '-aa')             { $param{'aa'}             = $args[++$i]; last CASE; }
         if ($args[$i] eq '-aaname')         { $param{'aaname'}         = $args[++$i]; last CASE; }
         if ($args[$i] eq '-bg')             { $param{'bg'}             = $args[++$i]; last CASE; }
         if ($args[$i] eq '-fg')             { $param{'fg'}             = $args[++$i]; last CASE; }
         if ($args[$i] eq '-lb')             { $param{'linebreak'}      = $args[++$i]; last CASE; }
         if ($args[$i] eq '-fs')             { $param{'fs'}             = $args[++$i]; last CASE; }
         if ($args[$i] eq '-gc')             { $param{'graph_char'}     = $args[++$i]; last CASE; }
         if ($args[$i] eq '-cm')             { $param{'cmethod'}        = $args[++$i]; last CASE; }
         die "[PARAM] Bad parameters! Argument '".$args[$i]."' unknown!\n";
      }
   }
   unless ($param{'graph_height'})  { $param{'graph_height'} = 33; }                   # default height of the conservation-rate graph
   unless ($param{'aln_groups'})    { $param{'aln_groups'}   = 1;  }                   # split groups are realigned by default
   unless ($param{'fs'})            { $param{'fs'}           = 8;  }                   # default font size for <PRE>-text
   if ($param{'bg'}) {                                                                 # if a background color was defined by the user, test for aliases!
      if ($param{'bg'} eq 'white')     {
         $param{'bg'} = '#FFFFFF';                                                     # aliases for background colors
         unless ($param{'fg'}) { $param{'fg'} = '#000000'; }                           # if no foreground was specified, use a default
      }
      if ($param{'bg'} eq 'black')     {
         $param{'bg'} = '#000000';
         unless ($param{'fg'}) { $param{'fg'} = '#FFFFFF'; }                           # if no foreground was specified, use a default
      }
      if ($param{'bg'} eq 'lightgray') {
         $param{'bg'} = '#DADADA';
         unless ($param{'fg'}) { $param{'fg'} = '#000000'; }                           # if no foreground was specified, use a default
      }
      if ($param{'bg'} eq 'darkgray')  {
         $param{'bg'} = '#505050';
         unless ($param{'fg'}) { $param{'fg'} = '#FFFFFF'; }                           # if no foreground was specified, use a default
      }
      if ($param{'bg'} eq 'midgray')   {
         $param{'bg'} = '#A0A0A0';
         unless ($param{'fg'}) { $param{'fg'} = '#000000'; }                           # if no foreground was specified, use a default
      }
   }
   unless ($param{'bg'})            { $param{'bg'} = '#FFFFFF';   }                    # default background color
   unless ($param{'fg'})            { $param{'fg'} = '#000000';   }                    # default foreground color
   unless ($param{'graph_char'})    { $param{'graph_char'} = '|'; }                    # default character used to build the graph bars with
   $param{'defline_len'} = 4;                                                          # default indention for graph and stuff (DO NOT CHANGE, or indention WILL fail!!!)
   return \%param;
}

sub test_inputfile_parameters {
   my $param   = $_[0];
   my @files   = split(/ /, $param->{'in'});

   if ((scalar(@files) > 1) && $param->{'groups'}) {
      print "\n[Commandline] Multiple Alignments *AND* Group-Defining is not allowed! Exiting...\n\n";
      exit 0;
   }
   if ($param->{'graph_char'} && (length($param->{'graph_char'}) > 1)) {
      print "\n[Commandline] Length of string for parameter '-gc' too long (must not be longer than 1 character)! Exiting...\n\n";
      exit 0;
   }
   if ($param->{'graph_height'} && ($param->{'graph_height'} > 100)) {
      print "\n[Commandline] A graph-height > 100 lines is not supported! Exiting...\n\n";
      exit 0;
   }
   $param->{'in'} = \@files;
   $param->{'numfiles'} = scalar(@files);
   return $param;                                                                      # just to be "clean" ;-)
}

sub test_group_params {
   my $incodes       = $_[0];
   my $param         = $_[1];
   my @nums          = split(/\;|\-|\,/, $param->{'groups'});
   my $nums_hash     = &IFG::NutsnBolts::array2hash(\@nums);
   my $error         = 0;

   @nums = sort(sort_num keys(%{$nums_hash}));
   foreach my $filename (keys(%{$incodes})) {
      if ($nums[-1] > scalar(@{$incodes->{$filename}})) {
         print "\n[Commandline] Alignment '$filename' has only ".scalar(@{$incodes->{$filename}})." sequences! Faulty group definition '".$param->{'groups'}."'...\n\n";
         $error = 1;
      }
   }
   return $error;
}

sub consensus_generation_and_formatting {
   my $in            = $_[0];
   my $incodes       = $_[1];
   my $suffix        = $_[2];  # will be appended to the filename if specified
   my $prevconss     = $_[3];  # eventually holds previously computed and formatted consensus sequences
   my $simplified    = $_[4];  # if specified (eq name_of_alphabet) VisCoSe will prepend consensus-names to the consensus sequence itself!
   my $param         = $_[5];


   if ($param->{'verbose'}) { print "Analysing aminoacid frequency in alignment(s):\n"; }
   my ($conserve, $rated, $cont) = &VisCoSe::ConsensusTools::get_AA_frequency_for_alignments($in, $incodes, $param);          # search for motifs / aminoacid frequencies per alignment-position
   if ($param->{'verbose'}) { print "done\n\n"; }

   if ($param->{'verbose'}) { print "Getting consensus sequence(s) "; }
   my $consensus = 0;
   CASE: {
      if ($param->{'cmethod'} && ($param->{'cmethod'} eq 'mean')) {
         if ($param->{'verbose'}) { print "- using using mean conservation rate:\n"; }
         $consensus = &VisCoSe::Motif::calc_mean_conservation_rate_per_file($conserve, $cont, $param);
         last CASE;
      }
      if ($param->{'cmethod'} && ($param->{'cmethod'} eq 'noise')) {
         if ($param->{'verbose'}) { print "- using noise measure:\n"; }
         $consensus = &VisCoSe::Motif::calc_noise_per_file($conserve, $cont, $param);
         last CASE;
      }
      if ($param->{'verbose'}) { print "- using relative majority:\n"; }
      ($consensus) = &VisCoSe::ConsensusTools::get_consensus_from_alignments($cont, $param);                                  # extracts the consensus sequence from conservation rates per column/AA
   }
   if ($param->{'verbose'}) { print "done\n\n"; }


   if ($param->{'verbose'}) { print "Formatting single consensus:\n"; }
   my ($conss, $consscode, $htmlfile, $seq) = VisCoSe::ConsensusTools::format_consensus_for_alignments($in, $incodes, $consensus, $conserve, $prevconss, $simplified, $param);
   if ($param->{'verbose'}) { print "done\n\n"; }

   if ($param->{'verbose'}) { print "Writing single consensus in HTML and FASTA format:\n"; }
   &VisCoSe::ConsensusTools::write_consensus_for_alignments($conss, $consscode, $htmlfile, $suffix, $param);                  # writes FASTA- and HTML-consensus for each *single* input dataset
   if ($param->{'verbose'}) { print "done\n\n"; }
   
   return ($conss, $consscode, $consensus, $seq);
}

sub consensus_comparison {
   my $conss         = $_[0];
   my $consscode     = $_[1];
   my $consensus     = $_[2];
   my $suffix        = $_[3];
   my $param         = $_[4];
   my $filename      = 'consensi_aligned';

   if ($param->{'verbose'}) { print "Generating consensus-dataset of ".$param->{'numfiles'}." files/groups:\n"; }
   my ($allconss, $allcodes) = &VisCoSe::AlignTools::generate_consensi_dataset($conss, $consscode, $param);                # $allconss, $allcodes: unaligned consensus!
   if ($param->{'verbose'}) { print "done\n\n"; }

   if ($param->{'verbose'}) { print "Aligning consensus-dataset:\n"; }
   if ((!($param->{'fra'})) && ($param->{'aligned_source'})) {
      print "   Not *RE*aligning on user-request ('-fra' not specified)!\n";
   }
   else {
      ($allconss, $allcodes) = &VisCoSe::AlignTools::align_consensus_dataset($allconss, $allcodes, $param);                # $allconss, $allcodes: realigned consensus!
   }
   if ($param->{'verbose'}) { print "done\n\n"; }

   if ($param->{'verbose'}) { print "Formatting Consensus-Comparison:\n"; }
   my $htmlfile = &VisCoSe::HTMLstuff::format_consensus_comparison($allconss, $allcodes, $consensus, $param);
   if ($param->{'linebreak'}) { $htmlfile = &VisCoSe::HTMLstuff::split_html_file($htmlfile, $param); }
   if ($suffix) { $filename .= '_'.$suffix; }
   &VisCoSe::ConsensusTools::write_html_consensus($htmlfile, $filename.'.html', $param);
   if ($param->{'verbose'}) { print "done\n\n"; }
}

### main ##########################################################################################

$| = 1;      # flushes buffer for prints etc. immediately!

if ((!(@ARGV)) || ($ARGV[0] eq "-h") || ($ARGV[0] eq "?") || ($ARGV[0] eq "help")) {
print <<EOF;
\nVisCoSe V1.0-10092003

SYNOPSIS:
viscose.pl [-h/-v] -in alignment.fasta [options...]

-h             : display this short description
-v             : be really noisy (-vv is even more noisy!)

Alignment options:
-------------------------------------------------------------------------------------------------------------------
-in <file>     : alignment(s) to analyse for motifs (FASTA/PEARSON/GDE format).
                 More than one alignment must be enclosed in "
-groups <def>  : defines groups of sequences (within a single alignment '-in'!) which will
                 be treated as single datasets (format: -groups "1-11,14;12,13,15-26;27-32")
                 (resulting 'group-datasets' will be realigned per default!)

Coloring options:
-------------------------------------------------------------------------------------------------------------------
-fra           : "force realignment"; if specified VisCoSe realigns the single groups after     [default:     off]
                 splitting *even if the source dataset was aligned*
-cg            : "colorgraph" - draw the graph in the same color as the consensus sequence      [default:     off]
-gh <num>      : "graph height" - height of the graph in characters                             [default:      33]
-sa/-ca/-fa    : "color alignment"/"simple alignment"/"full alignment"                          [default:     off]
                 - prints the alignment on top of the consensus, "-ca" in color, "-sa" in
                   black&white, "-fa" prints all characters in it's conservation rate color

Viewing options:
-------------------------------------------------------------------------------------------------------------------
-bg <color>    : background color (accepts "white", "black", "lightgray", "midgray",            [default:   white]
                 "darkgray" or color definitions in *HTML code*)
-fg <color>    : foreground color (accepts "white", "black", "lightgray", "midgray",            [default:   black]
                 "darkgray" or color definitions in *HTML code*) - a specification of a
                 foreground color is not nessecary when using the standard aliases since
                 it's determined automatically by VisCoSe basing on the background color
-fs <num>      : size of the fixed width font (-fs 6 / -lb 150 for printing!)                   [default:       8]
-gc <char>     : "graphchar" - specifies which character is to be used for building the         [default:     "|"]
                 conservation graph. CAREFUL: length = 1 !!! Longer strings will mess up
                 the graph!
-lb <num>      : introduces linebreaks for the alignment and the consensus sequence at          [default:     off]
                 position "num"


Alphabet/Consensus options:
-------------------------------------------------------------------------------------------------------------------
-aa <file>     : filename with the configuration data for simplified amino acid alphabets
                 (if none is provided, VisCoSe uses the standard alphabet only)
-aaname <name> : name of the alphabet which should be used for translating
-cm <name>     : "conservation method" - chooses the method used to calculate the rate          [default:  relmaj]
                 for the histogram:
                  relmaj: "relative majority", the highest conservation rate is displayed
                  mean  : calculates mean conservation rate per column

(c) 2002/2003 Michael Spitzer (Integrated Functional Genomics, Muenster, Germany)\n
EOF
exit 0;
}

my ($in, $incodes, $trans, $alphabets) = (0, 0, 0, 0);
my $param = &parse_args(\@ARGV);
$param    = &test_inputfile_parameters($param);

### reading and configuring simplified amino acid alphabets
if ($param->{'aa'}) {
   if ($param->{'verbose'}) { print "\nReading alphabet configuration file '".$param->{'aa'}."'..."; }
   $alphabets = VisCoSe::Simplified::read_alphabets($param);
   VisCoSe::Simplified::verify_alphabet_groups_for_ambiguity($alphabets, $param);
   if ($param->{'verbose'}) { print "done\n"; }
}
### reading and configuring simplified amino acid alphabets done

### Bit::Vector Configuration
Bit::Vector->Configuration("in=enum");
Bit::Vector->Configuration("semantic=arithmetic");
Bit::Vector->Configuration("out=enum");
### Bit::Vector Configuration done

if ($param->{'verbose'}) { print "\nReading datasets:\n"; }
($in, $incodes) = &VisCoSe::ConsensusTools::read_alignments($param);                                              # read *all* datasets
if (($param->{'numfiles'} == 1) && $param->{'groups'}) {
   if (test_group_params($incodes, $param)) { exit 0; }
}
if ($param->{'verbose'}) { print "done\n\n"; }

if ($param->{'groups'}) {
   if ($param->{'verbose'}) { print "Splitting dataset into groups:\n"; }
   ($in, $incodes, $param) = &VisCoSe::AlignTools::split_dataset_into_groups($in, $incodes, $param);              # get corresponding sequences for all defined groups
   &VisCoSe::AlignTools::write_group_alignments($in, $incodes, $param);                                           # write the resulting (degapped!) "group-datasets"
   if ($param->{'verbose'}) { print "done\n\n"; }
}

if ($param->{'verbose'}) { print "Testing dataset(s) for status:\n"; }
($in, $incodes) = &VisCoSe::AlignTools::test_for_aligned_data($in, $incodes, $param);                             # Test dataset(s) for gaps to determine aligned/unaligned state of dataset(s)
if ($param->{'verbose'}) { print "done\n\n"; }

if ($param->{'verbose'}) { print "Aligning dataset(s):\n"; }
($in, $incodes) = &VisCoSe::AlignTools::align_unaligned_datasets($in, $incodes, $param);                          # align any unaligned datasets (if groups were defined or unaligned datasets were detected by VisCoSe::AlignTools::test_alignments_for_gaps())
if ($param->{'verbose'}) { print "done\n\n"; }

my ($conss, $consscode, $consensus, $prevconss) = consensus_generation_and_formatting($in, $incodes, '', '', 0, $param);

my ($conss_smpl, $consscode_smpl, $consensus_smpl, $prevconss_smpl) = (0, 0, 0);
if ($param->{'aa'}) {
   if ($param->{'verbose'}) { print "Translating using alphabet '".$param->{'aaname'}."':\n"; }
   $in = VisCoSe::Simplified::translate_alignments($in, $incodes, $alphabets->{$param->{'aaname'}}, $param);
   if ($param->{'verbose'}) { print "done\n\n"; }
   ($conss_smpl, $consscode_smpl, $consensus_smpl, $prevconss_smpl) = consensus_generation_and_formatting($in, $incodes, 'simplified_'.$param->{'aaname'},
                                                                                                          $prevconss, $param->{'aaname'}, $param);
}

if (($param->{'numfiles'} > 1) || ($param->{'groups'})) {
   consensus_comparison($conss, $consscode, $consensus, '', $param);
}

if ($param->{'cgi'}) {                                                                                            # if the script was started by my CGI-Skript, write a file indicating that we finished!
   open(_file_, ">viscose_finished.txt");
   print _file_ "blah";
   close(_file_);
}

exit 0;
