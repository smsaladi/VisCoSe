#! /usr/bin/perl -w

# interface.pl
# (c) 2003 Michael Spitzer, IFG-IZKF
# (c) 2017 Shyam Saladi, Caltech

use diagnostics;                # verbose warning- and error-messages.
use strict;                     # This is a "must-have".
#use CGI::Carp qw(fatalsToBrowser);
use CGI qw/:standard/;
use POSIX;
use File::Basename;
use Data::Dumper;
use FindBin qw($Bin);

use constant debug => 1;        # set to '1' for debug messages

use constant temp_dir    => $Bin . '/../temp';
use constant cgi_dir     => $Bin;
use constant mafft_bin   => $Bin . '/mafft';
use constant link_url    => 'https://viscose.herokuapp.com/temp/';

$ENV{'MAFFT_BINARIES'}  = mafft_bin;
sub print_header {
    print 'Content-type: text/html', "\n\n", '<HTML>', "\n", '<HEAD>', "\n", '<TITLE>Status of VisCoSe: interface.pl</TITLE>', "\n", '</HEAD>', "\n", '<PRE>', "\n", 'Status of VisCoSe: interface.pl (interface script for consensus visualization):', "\n",
      '--', "\n\n";
}


sub print_footer {
    print '</PRE>', "\n", '</HTML>', "\n";
}


sub extract_data {
    my $err      = 0;
    my %param    = ();
    my @gmt      = gmtime();
    my $date     = ( 1900 + $gmt[5] ) . '_';            # render date-string
    my $boundary = 0;

    ( length( $gmt[4] ) == 1 ) ? ( $date .= '0' . ( $gmt[4] + 1 ) . '_' ) : ( $date .= ( $gmt[4] + 1 ) . '_' );
    ( length( $gmt[3] ) == 1 ) ? ( $date .= '0' . ( $gmt[3] ) . '_' )     : ( $date .= ( $gmt[3] ) . '_' );
    ( length( $gmt[2] ) == 1 ) ? ( $date .= '0' . ( $gmt[2] + 2 ) . '_' ) : ( $date .= ( $gmt[2] + 2 ) . '_' );
    ( length( $gmt[1] ) == 1 ) ? ( $date .= '0' . ( $gmt[1] ) . '_' )     : ( $date .= ( $gmt[1] ) . '_' );
    $param{'dir'} = $date . basename(POSIX::tmpnam);                                            # prepare job-ID (job-ID is saved in $param{'dir'} !!!)

    if (debug) { print Dump(); }

    if (param('alignment')) {
        $param{'alignment'} = param('alignment');
        $err++;                                                                                 # used to check a doubled alignment via upload and copy&paste
        if (debug) { print ">>> Alignment found (textbox):\n".$param{'alignment'}."'\n"; }
    }

    if (param('localfile')) {
        $param{'source_filename'} = param('localfile');
        ($param{'source_filename'}) = $param{'source_filename'} =~ m#^(.*)$#;                   # "naive" untainting of the filename
        my $fh                    = upload('localfile');
        $param{'alignment'}       = '';
        while (my $line = <$fh>) { $param{'alignment'} .= $line; }
        $err++;                                                                                 # used to check a doubled alignment via upload and copy&paste
        if (debug) { print ">>> Alignment found (upload '".$param{'source_filename'}."'):\n".$param{'alignment'}."\n"; }
    }

    if (param('groups')) {
        $param{'groups'} = param('groups');
        if (debug) { print ">>> Groups found: '".$param{'groups'}."'\n"; }
    }

    if (param('gh' )) {
        $param{'gh'} = param('gh' );
        if (debug) { print ">>> Graph-Height found: '".$param{'gh'}."'\n"; }
    }

    if (param('fra')) {
        $param{'fra'} = 1;
        if (debug) { print ">>> FRA (Force re-alignment) activated\n"; }
    }

    if (param('cg')) {
        $param{'cg'} = 1;
        if (debug) { print ">>> CG (Colored graph) activated\n"; }
    }

    if (param('radiogroup')) {
        $param{param('radiogroup')} = 1;
        if (debug) { print ">>> Alignment coloring: ".$param{param('radiogroup')}."\n"; }
    }

    if (param('bg')) {
        $param{'bg'} = param('bg');
        if (debug) { print ">>> BG (Background) found: '".$param{'bg'}."'\n"; }
    }

    if (param('fg')) {
        $param{'fg'} = param('fg');
        if (debug) { print ">>> FG (Foreground) found: '".$param{'fg'}."'\n"; }
    }

    if (param('fs' )) {
        $param{'fs'} = param('fs' );
        if (debug) { print ">>> FS (Fontsize) found: '".$param{'fs'}."'\n"; }
    }

    if (param('lb')) {
        $param{'lb'} = 1;
        if (debug) { print ">>> LB (Linebreak) found\n"; }
    }

    if (param('lbpos')) {
        $param{'lbpos'} = param('lbpos');
        if (debug) { print ">>> LBPOS (Linebreak position) found: '".$param{'lbpos'}."'\n"; }
    }

    if (param('gc')) {
        $param{'gc'} = param('gc');
        if (debug) { print ">>> GC (Graph character) found: '".$param{'gc'}."'\n"; }
    }

    if (param('alphabet') && (param('alphabet') ne 'none')) {
        $param{'aa'}     = cgi_dir . '/alphabets.dat';
        $param{'aaname'} = param('alphabet');
        if (debug) { print ">>> ALPHABET found: '".$param{'aaname'}."'\n"; }
    }

    if ( $err >= 2 ) { print "Textbox alignment *AND* upload alignment found. Ignoring textbox alignment!\n"; }
    unless ( $param{'source_filename'} ) { $param{'source_filename'} = 'src_alignment.fasta'; }

    return \%param;
}


sub get_number_of_sequences {
    my $param  = $_[0];
    my @aln    = split( /\n/, $param->{'alignment'} );
    my $numseq = 0;

    foreach my $line (@aln) {
        if ( $line =~ m/\A(?:\>|\%)/ ) { $numseq++; }
    }
    return $numseq;
}


sub test_parameters {
    my $param   = $_[0];
    my $numseq  = &get_number_of_sequences($param);
    my $dummy   = 0;
    my @numbers = ();
    my $error   = 0;

    foreach my $key ( keys( %{$param} ) ) {
      CASE: {
            if ( $key eq 'alignment' ) {

                # no sanity check here 'coz alignment data and deflines can be *very* variable
                last CASE;
            }
            if ( $key eq 'groups' ) {
                if ( $param->{$key} =~ m/[^0-9\,\-\;]/ ) {    # check for illegal characters in group-definition (allowed: [0-9,-;] )
                    print "\n>>> Illegal characters in groups-definition '" . $param->{$key} . "'! Only 0-9 and , and - and ; are allowed!\n";
                    $error = 1;
                }
                else {
                    $dummy = $param->{'groups'};              # save the group-definition to a temporary variable
                    $dummy =~ tr/-,;/:/;                      # translate all valid special characters to ':'
                    @numbers = split( /\:/, $dummy );         # split at ':', so we get all the numbers and we can easily look at each of them
                    foreach my $number (@numbers) {
                        if ( $number > $numseq ) {
                            print "\n>>> Group-range-number '$number' is too high! Alignment has only $numseq sequences!\n";
                            $error = 1;
                        }
                    }
                }
                last CASE;
            }
            if ( $key eq 'gh' ) {
                if ( $param->{$key} =~ m/[^0-9]/ ) {
                    print "\n>>> Illegal characters in graph-height '" . $param->{$key} . "'! Only 0-9 are allowed!\n";
                    $error = 1;
                }
                last CASE;
            }
            if ( $key eq 'lbpos' ) {
                if ( $param->{$key} =~ m/[^0-9]/ ) {
                    print "\n>>> Illegal characters in linebreak position '" . $param->{$key} . "'! Only 0-9 are allowed!\n";
                    $error = 1;
                }
                last CASE;
            }
            if ( $key eq 'gc' ) {
                if ( length( $param->{$key} ) > 1 ) {
                    print "\n>>> Illegal length of graph character/string '" . $param->{$key} . "'! Length must not exceed '1'!\n";
                    $error = 1;
                }
                last CASE;
            }
        }
    }
    unless ( $param->{'alignment'} ) {
        print "\n>>> You must at least provide an alignment!\n";
        $error = 1;
    }
    return $error;
}


sub print_parameters {
    my $param = $_[0];

    print "\nParameter status information:\n=============================\n";
    print "Job ID          : " . $param->{'dir'} . "\n";
    foreach my $key ( keys( %{$param} ) ) {
      CASE: {
            if ( $key eq 'alignment' ) {
                if ( $param->{'source_filename'} ) { print "Alignment       : present (uploaded: '" . $param->{'source_filename'} . "')\n"; last CASE; }
                else { print "Alignment     : present / copy&paste\n"; last CASE; }
            }
            if ( $key eq 'groups' ) { print "Groups          : " . $param->{$key} . "\n";            last CASE; }
            if ( $key eq 'fra' )    { print "Force ReAlign   : enabled, if pre-aligned!\n";          last CASE; }
            if ( $key eq 'gh' )     { print "Graph Height    : " . $param->{$key} . " characters\n"; last CASE; }
            if ( $key eq 'sa' )     {
                unless ( $param->{'ca'} || $param->{'fa'} ) { print "Show Alignment  : enabled / simple\n"; }
                last CASE;
            }
            if ( $key eq 'ca' ) {
                unless ( $param->{'fa'} ) { print "Show Alignment  : enabled / colored\n"; }
                last CASE;
            }
            if ( $key eq 'fa' ) { print "Show Alignment  : enabled / full coloring\n";              last CASE; }
            if ( $key eq 'cg' ) { print "Graph           : colored\n";                              last CASE; }
            if ( $key eq 'bg' ) { print "Background Color: " . $param->{$key} . "\n";               last CASE; }
            if ( $key eq 'fg' ) { print "Foreground Color: " . $param->{$key} . "\n";               last CASE; }
            if ( $key eq 'fs' ) { print "Font Size       : " . $param->{$key} . " pts\n";           last CASE; }
            if ( $key eq 'lb' ) { print "Linebreak After : " . $param->{'lbpos'} . " characters\n"; last CASE; }
            if ( $key eq 'gc' ) { print "Graph Character : '" . $param->{'gc'} . "'\n";             last CASE; }
            if ( $key eq 'aa' ) { print "Alphabet        : '" . $param->{'aaname'} . "'\n";         last CASE; }
        }
    }
    print "\n";
    mkdir( temp_dir . '/' . $param->{'dir'} ) or die ">>> [IOErr] Could not create temporary directory for alignment data (job-ID '" . $param->{'dir'} . "')! Contact spitzem\@uni-muenster.de!\n";
    chdir( temp_dir . '/' . $param->{'dir'} );
    write_temp_resultpage($param);
    print 'Temporary result page (may be used for bookmarking prior to job completion):' . "\n" . '<A HREF="' . link_url . '/' . $param->{'dir'} . '/viscose_results.html" target="_blank">viscose_results.html</A>' . "\n\n\n\n";
    return $param;
}


sub write_temp_resultpage {
    my $param = $_[0];

    my @html = (
        '<HTML>',
        '<HEAD>',
        '<TITLE>VisCoSe-Results for job-ID ' . $param->{'dir'} . '</TITLE>',
        '</HEAD>',
        '<font face="Arial, Helvetica, sans-serif">',
        '<TABLE BORDER=0 CELLSPACING="2" CELLPADDING="4" WIDTH="100%">',
        '   <TR bgcolor="#EEEEFF">',
        '      <TD>',
        '         <H3>VisCoSe-Results for job-ID "' . $param->{'dir'} . '":</H3>',
        '      </TD>',
        '   </TR>',
        '</TABLE>',
        '<P>',
        '<TABLE BORDER=0 CELLSPACING="2" CELLPADDING="4" WIDTH="100%">',
        '   <TR bgcolor="#DDDDFF">',
        '      <TD><font face="Arial, Helvetica, sans-serif" size="2">',
        '         Your VisCoSe-job is not ready yet...<BR><BR>Please wait until this page shows your results automatically<BR>or manually reload the page in your browser.',
        '      </TD>',
        '   </TR>',
        '</TABLE>',
        '</P>',
        '<P>',
        '<TABLE BORDER=0 CELLSPACING="0" CELLPADDING="4" WIDTH="100%">',
        '   <TR bgcolor="#EEEEFF">',
        '      <TD><font face="Arial, Helvetica, sans-serif">',
        '         This page refreshes every 30 seconds',
        '      </TD>',
        '   </TR>',
        '</TABLE>',
        '</P>',
        '<meta http-equiv="refresh" content="30; URL=' . link_url . '/' . $param->{'dir'} . '/viscose_results.html">',
        '</HTML>',
    );
    open( _file_, '>' . temp_dir . '/' . $param->{'dir'} . '/viscose_results.html' ) or die ">>> [IOErr] Could not create temporary result page (job-ID '" . $param->{'dir'} . "')! Contact spitzem\@uni-muenster.de!\n";
    foreach my $line (@html) { print _file_ $line . "\n"; }
    close(_file_);
}


sub write_alignment_to_disk {
    my $param = $_[0];

    open( _aln_, '>'.temp_dir.'/'.$param->{'dir'}.'/'.$param->{'source_filename'}) or die ">>> [IOErr] Could not write alignment data (job-ID '".$param->{'dir'}."')! Contact spitzem\@uni-muenster.de!\n";
    print _aln_ $param->{'alignment'};
    close(_aln_);
}


sub start_viscose_script {
    my $param = $_[0];
    my ( $shellstring, $log ) = ( '', '' );
    my $seconds = 0;

    $shellstring = cgi_dir . '/viscose.pl -v -cgi -in ' . $param->{'source_filename'} . ' ';
    foreach my $key ( keys( %{$param} ) ) {
      CASE: {
            if ( $key eq 'groups' ) { $shellstring .= '-groups ' . "'" . $param->{'groups'} . "'" . ' '; last CASE; }
            if ( $key eq 'fra' )    { $shellstring .= '-fra ';                                           last CASE; }
            if ( $key eq 'gh' )     { $shellstring .= '-gh ' . $param->{'gh'} . ' ';                     last CASE; }
            if ( $key eq 'sa' )     {
                unless ( $param->{'ca'} || $param->{'fa'} ) { $shellstring .= '-sa '; }    # '-ca' or '-fa' overrides '-sa'
                last CASE;
            }
            if ( $key eq 'ca' ) {
                unless ( $param->{'fa'} ) { $shellstring .= '-ca '; }                      # '-fa' overrides '-ca'
                last CASE;
            }
            if ( $key eq 'fa' )     { $shellstring .= '-fa ';                                  last CASE; }
            if ( $key eq 'cg' )     { $shellstring .= '-cg ';                                  last CASE; }
            if ( $key eq 'bg' )     { $shellstring .= '-bg ' . $param->{'bg'} . ' ';           last CASE; }
            if ( $key eq 'fg' )     { $shellstring .= '-fg ' . $param->{'fg'} . ' ';           last CASE; }
            if ( $key eq 'fs' )     { $shellstring .= '-fs ' . $param->{'fs'} . ' ';           last CASE; }
            if ( $key eq 'lb' )     { $shellstring .= '-lb ' . $param->{'lbpos'} . ' ';        last CASE; }
            if ( $key eq 'gc' )     { $shellstring .= '-gc "' . $param->{'gc'} . '" ';         last CASE; }
            if ( $key eq 'aa' )     { $shellstring .= '-aa "' . $param->{'aa'} . '" ';         last CASE; }
            if ( $key eq 'aaname' ) { $shellstring .= '-aaname "' . $param->{'aaname'} . '" '; last CASE; }
        }
    }
    if (debug) { print "\n\ncommandline: '$shellstring'\n"; }

    ($shellstring) = $shellstring =~ m#^(.*)$#;    # untaint '$shellstring' the *INSECURE* way coz we know what's inside!!! Otherwise correct untainting would be advised!
    system($shellstring . ' | tee viscose.log');
    system(cgi_dir . '/make_result_page.pl -dir ' . $param->{'dir'} . ' -src ' . $param->{'source_filename'});

    $log         = `$shellstring`;
    if (debug) { chomp($log); print "Shellmessage: '$log'\nShellstring: '$shellstring'\n"; }

    until ( -s 'viscose_finished.txt' ) {                                           # test if the main result file has been written
        sleep(10);                                                                  # wait 10 seconds between each test
        $seconds += 10;                                                             # increase our timer by 10
        if ( ( $seconds % 60 ) == 0 ) { print '.'; }                                # At each minute of waiting time print a dot to the browser to avoid a timeout
    }

    $log = `rm viscose_finished.txt`;
    chdir(cgi_dir);

    return $seconds;
}


sub make_links {
    my $param = $_[0];

    print 'Results can be found at:' . "\n";
    print '<A HREF="' . link_url . '/' . $param->{'dir'} . '/viscose_results.html">viscose_results.html</A>' . "\n\n";
    print 'You will be forwarded to the results page in 5 seconds...' . "\n";
    print '<meta http-equiv="refresh" content="-1; URL=' . link_url . '/' . $param->{'dir'} . '/viscose_results.html">' . "\n";
}

### main ##########################################################################################

if (debug) {    # redirect STERR output so it is visible in the browser
    open( STDERR, ">&STDOUT" ) or print "Could not redirect STDERR!\n";
}

$| = 1;
umask(0000);    # set 'rw(x)' for all user, group, other

&print_header();    # print correct HTML header information

if (debug) { print "Extracting parameters from FORM data..."; }
my $param = &extract_data();             # get the single data from all FORM fields
if (debug) { print "done\n"; }

if (debug) { print "Testing parameters for integrity..."; }
my $error = &test_parameters($param);    # test parameters for sanity
if (debug) { print "done\n"; }

unless ($error) {                        # only continue if the parameter-testing was finished w/o any errors!
    &print_parameters($param);           # print status information about used parameters

    if (debug) { print "Writing source alignment to disk..."; }
    &write_alignment_to_disk($param);
    if (debug) { print "done\n"; }

    print "Running VisCoSe (runtime depends on dataset size and the necessity of (re-)alignment)\n-> printing one dot per minute: ";
    my $seconds = &start_viscose_script($param);
    print "done in $seconds seconds (rounded to a 10s)\n\n";

    &make_links($param);
}
else {
    print "\n>>> STOPPING PIPELINE EXECUTION BECAUSE OF INPUT ERRORS! PLEASE ADJUST YOUR PARAMETERS! <<<\n";
}

&print_footer();

exit(0);
