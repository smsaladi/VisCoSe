###################################################################################################
# IFG::NutsnBolts.pm
# (c) 2001/2002 Michael Spitzer, IFG Muenster (Integrated Funtional Genomics), Germany
#
# 15.10.2003   read_config() updated / bugfixed
#
# several subroutines that may be useful someday, somehow, for someone...
# update: ok, they are indeed useful, especially for me... ;-)

package IFG::NutsnBolts;
require Exporter;

our @EXPORT  = qw(cnv_seqnum_to_sequential_fmt        cnv_seqnum_to_range_fmt          get_accession_from_description      round
                  write_array                         write_matrix                     derivate                            smooth
                  log10                               log_n                            a_mean                              q_mean
                  g_mean                              h_mean                           sd                                  sqr
                  array2hash                          read_config                      cut_file_extension                  get_directory_content
                  bzip2                               bunzip2
                 );
		  
our $VERSION = "06.09.2002";

sub sort_num { $a <=> $b }

sub cnv_seqnum_to_sequential_fmt {
# takes a sequence number from a range of ranges and calculates the corresponding sequential number
   my $seqnum     = $_[0];                    # sequence number in range format
   my @range      = split(",", $_[1]);        # Range information, e.g. "0-33,35-55,..."
   my $i          = 0;
   my $found      = 0;
   
   foreach my $range (@range) {
      $range =~ m/(\d+)\-(\d+)/;
      for (my $j = ($1-1); $j < ($2); $j++) {
         if ($j == $seqnum) { $found = 1; }
	 last if ($found == 1);
         $i++;
      }
      last if ($found == 1);
   }
   return $i;
}

sub cnv_seqnum_to_range_fmt {
# takes a sequential sequence number and calculates the corresponding sequence number in a range of ranges
   my $seqnum     = $_[0];                    # sequence number in sequential format
   my @range      = split(",", $_[1]);        # Range information, e.g. "0-33,35-55,..."
   my $i          = 0;
   my $found      = 0;
   
   foreach my $range (@range) {
      $range =~ m/(\d+)\-(\d+)/;
      for (my $j = ($1-1); $j < ($2); $j++) {
         if ($i == $seqnum) { $found = $j; }
	 last if ($found != 0);
         $i++;
      }
      last if ($found != 0);
   }
   return $found;
}

sub get_accession_from_description {
# parses a description-line (e.g. from a blast report) and extracts the accession code (e.g. NM_xxxxxx)
   my $description      = $_[0];
   my $accession        = "";

   $description =~ m/\(([\w\d]+)\)/;      # normally, accession-codes are in between '(...)' and consist of A-z, 0-9, _
   return $1 || "NA";                     # if an accession was found ($1 == defined) return it, else "NA"
}

sub round {
# from perlmonks.org: by 'btrott' on Apr 25, 2000 (http://perlmonks.org/index.pl?node_id=8786)
   my $number        = $_[0];
   my $precision     = $_[1];

   sprintf "%.".$precision."f", $number;
}

sub write_array {
   my $filename               = $_[0];
   my $array                  = $_[1];

   open(_file_, ">$filename") or die "[IOErr] Could not open '$filename' for write access!\n";
   foreach my $value (@{$array}) { print _file_ $value."\n"; }
   close(_file_);
}

sub write_matrix {
   my $filename               = $_[0];
   my $array                  = $_[1];

   open(_file_, ">$filename") or die "[IOErr] Could not open '$filename' for write access!\n";
   for (my $i = 0; $i < scalar(@{$array}); $i++) {
      for (my $j = 0; $j < scalar(@{$array->[$i]}); $j++) {
         if ($j < (scalar(@{$array->[$i]}) - 1)) { print _file_ $array->[$i]->[$j].'   '; }
         else { print _file_ $array->[$i]->[$j]; }
      }
      print _file_ "\n";
   }
   close(_file_);
}

sub derivate {
   my $array                  = $_[0];
   my $stepwidth              = $_[1];
   my $num                    = scalar(@{$array});
   my @delta                  = ();
   my $delta                  = 0;

   for (my $i = 0; $i < $num; $i += $stepwidth) {
      if (($i + $stepwidth) >= $num) { $delta = $array->[$i] - $array->[$num-1]; }
      else { $delta = $array->[$i] - $array->[$i + $stepwidth]; }
      push(@delta, $delta);
      last if (($i + $stepwidth) >= $num);
   }
   return \@delta;
}

sub smooth {
   my $array                  = $_[0];
   my $stepwidth              = $_[1];
   my @smoothed               = ();
   my ($upper, $sum, $num)    = (0, 0, 0);

   for (my $i = 0; $i < scalar(@{$array}); $i++) {
      if (($i + $stepwidth) > scalar(@{$array})) { $upper = scalar(@{$array}); }
      else { $upper = $i + $stepwidth; }
      ($sum, $num) = (0, 0); 
      for (my $j = $i; $j < $upper; $j++) {
         $sum += $array->[$j];                                                            # summarize all values in the range
         $num++;                                                                          # count the number of values in the range
      }
      push(@smoothed, ($sum / $num));                                                     # store the mean value
   }
   return \@smoothed;
}

sub log10 {
   return (log($_[0]) / log(10));
}
                   
sub log_n {
   return (log($_[0]) / log($_[1]));                                                      # "$_[1]" is the base!
}
                   
sub a_mean { # calculates the arithmetical means 
   my $array         = $_[0];
   my $num           = scalar(@{$array});
   my $sum           = 0;

   for (my $i = 0; $i < $num; $i++) { $sum += $array->[$i]; }
   return ($sum / $num);
}

sub q_mean { # calculates the quadratic means
   my $array         = $_[0];
   my $num           = scalar(@{$array});
   my $sum           = 0;

   for (my $i = 0; $i < $num; $i++) { $sum += sqr($array->[$i]); }
   return sqrt($sum / $num);
}

sub g_mean { # calculates the geometric means >>> ZERO if one of the elements in the array is zero!
   my $array         = $_[0];
   my $num           = scalar(@{$array});
   my $mult          = $array->[0];

   for (my $i = 1; $i < $num; $i++) { $mult *= $array->[$i]; }
   #print "MULT: $mult   -> ";
   return ($mult ** (-($num)));
}

sub h_mean { # calculates the harmonic means >>> UNDEF if one of the elements in the array is zero!
   my $array         = $_[0];
   my $num           = scalar(@{$array});
   my ($mult, $sum)  = ($array->[0], 0);

   for (my $i = 1; $i < $num; $i++) { $mult *= $array->[$i]; }
   for (my $i = 0; $i < $num; $i++) {
      for (my $j = ($i + 1); $j < $num; $j++) {
         $sum += $array->[$i] * $array->[$j];
      }
   }
   #print "MULT: $mult / SUM: $sum   -> ";
   return (($num * $mult) / $sum);
}

sub sd { # calculates the standard deviation
   my $array         = $_[0];
   my $a_mean        = a_mean($array);
   my $num           = scalar(@{$array});
   my $sum           = 0;

   for (my $i = 0; $i < $num; $i++) { $sum += sqr($array->[$i] - $a_mean); }
   return sqrt($sum / $num);
}

sub sqr { #calculates the square of a value
   return ($_[0] * $_[0]);
}

sub array2hash {
   my %hash    = ();
   
   @hash{@{$_[0]}} = (1) x @{$_[0]};
   return \%hash;
}

sub read_config {
   my $filename                                          = $_[0];
   my $param                                             = $_[1];
   my %config                                            = ();
   my ($line, $section, $parameter, $value, $ioerr)      = (0, 0, 0, 0, 0);

   open(_file_, "<$filename") or $ioerr = 1;
   unless ($ioerr) {                                                                                  # don't continue if the file couldn't be opened
      until (eof(_file_)) {
         $line = <_file_>;
         chomp($line);
         $line =~ s/\#.*//;                                                                           # remove comments from line if there's any
         if ($line =~ m/\[(.+)\]/) {                                                                  # "[xyz]" ? -> New config file section!
            $section = $1;                                                                            # remember name of the config file section
            if ($param->{'debug'}) { print ">>> [CONFIG] '[$section] section'\n"; }
         }
         elsif ($line =~ m/([A-Za-z0-9\_\.\-\" ]+)(?:\s)*\=(?:\s)*(.+)/) {                            # read the parameter specification {parameter_name} {=} {parameter_value}
            $parameter = $1;
            $value     = $2;
            $parameter =~ s/\s*\Z//;                                                                  # remove all trailing whitespaces from the parameter name
            $value     =~ s/\s*\Z//;                                                                  # remove all trailing whitespaces from the parameter value
            if ($param->{'debug'}) { print ">>>    [CONFIG] '$parameter'\t- '$value'\n"; }
            $config{$section}->{$parameter} = $value;
         }
      }
      close(_file_);
   }
   return \%config;                                                                                   # on file-open-error, an empty hash is returned!
}

sub cut_file_extension {
   my $filename   = $_[0];

   $filename =~ s/(?:\.[^.]+)\Z//;
   return $filename;
}

sub get_directory_content {
   my $dir        = $_[0];
   my @cleandir   = ();

   opendir(_dir_, $dir);
   my @files = readdir(_dir_);
   closedir(_dir_);
   foreach my $file (@files) {
      if (($file ne '.') && ($file ne '..')) { push(@cleandir, $file); }
   }
   return \@cleandir;
}

sub bzip2 {
   my $file    = $_[0];

   my $log = `bzip2 -v9 $file 2>/dev/null`;
}

sub bunzip2 {
   my $file    = $_[0];

   my $log = `bzip2 -vd $file 2>/dev/null`;
}
