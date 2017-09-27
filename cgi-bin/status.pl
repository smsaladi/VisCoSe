#! /usr/bin/perl -wT
####################################################################################################
## status.pl
#
use diagnostics;    # verbose warning- and error-messages.
use strict;         # This is a "must-have".
#use CGI::Carp qw(fatalsToBrowser);
use CGI qw/:standard/;

sub _remove_grep_lines {
    my @log = split(/\n/, $_[0]);
    my $new = '';

    foreach my $line (@log) {
        if ($line !~ m/grep/) { $new .= "$line\n"; }
    }
    return ($new || 'none');
}
####################################################################################################

$ENV{'PATH'} = '/bin:/usr/bin';

print header(),
      start_html('Status of VisCoSe server'),
      h1('Status of the VisCoSe server');


my $shellstring  = "ps -Af |grep 'viscose.pl'";
my $viscose_jobs = `$shellstring`;
   $viscose_jobs = _remove_grep_lines($viscose_jobs);
   $shellstring  = "ps -Af |grep 'dvtditr'";
my $aln_jobs1    = `$shellstring`;
   $aln_jobs1    = _remove_grep_lines($aln_jobs1);
   $shellstring  = "ps -Af |grep 'fftns'";
my $aln_jobs2    = `$shellstring`;
   $aln_jobs2    = _remove_grep_lines($aln_jobs2);

print br(), 'Currently active VisCoSe jobs:',   br(), pre($viscose_jobs), br(),
      br(), 'Currently active alignment jobs:', br(), pre($aln_jobs1), pre($aln_jobs2);

print end_html();

exit 0;

