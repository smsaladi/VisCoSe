#! /usr/bin/perl -w

# make_result_page.pl
# (c) 2003 Michael Spitzer, IFG-IZKF
# (c) 2017 Shyam Saladi, Caltech

use diagnostics;                          # verbose warning- and error-messages.
use strict;                               # This is a "must-have".
use CGI::Carp qw(fatalsToBrowser);
use Data::Dumper;
use FindBin qw($Bin);

use constant debug => 0;                  # set to '1' for debug messages
use constant temp_dir    => $Bin . '/../temp';
use constant cgi_dir     => $Bin;
use constant mafft_bin   => $Bin . '/mafft';
use constant link_url    => 'https://viscose.herokuapp.com/temp/';

sub parse_args {
# parses the arguments given via command line and builds a hash which will be returned.
# Does a little error checking and dies if something's wrong with the parameters (eg. some switch is unknown or so).
   my @args          = @{$_[0]};        # holds @ARGV
   my %param         = ();              # will be filled with the parameters from @ARGV

   for (my $i = 0; $i < scalar @args; $i++) {
      CASE: {
         if ($args[$i] eq '-dir')            { $param{'dir'}             = $args[++$i]; last CASE; }
         if ($args[$i] eq '-src')            { $param{'source_filename'} = $args[++$i]; last CASE; }
         die "[PARAM] Bad parameters! Argument '".$args[$i]."' unknown!\n";
      }
   }
   return \%param;
}

sub make_links {
   my $param                         = $_[0];
   my ($color1, $color2, $htmlcolor) = ('#DDDDFF', '#CCCCFF', '#DDDDFF');
   my $dummy                         = '';
   my @html                          = ();

   opendir(_dir_, temp_dir.'/'.$param->{'dir'}) or die "Could not open result directory! Weird...please contact spitzem\@uni-muenster.de!\n";
   my @dir = readdir(_dir_);                                                  # read directory content from result directory
   closedir(_dir_);
   @html = ('<HTML>',
            '<HEAD>',
            '<TITLE>VisCoSe-Results for job-ID '.$param->{'dir'}.'</TITLE>',
            '</HEAD>',
            '<font face="Arial, Helvetica, sans-serif">',
            '<TABLE BORDER=0 CELLSPACING="2" CELLPADDING="4" WIDTH="100%">',
            '   <TR bgcolor="#EEEEFF">',
            '      <TD>',
            '         <H3>VisCoSe-Results for job-ID "'.$param->{'dir'}.'":</H3>',
            '      </TD>',
            '   </TR>',
            '</TABLE>',
            '<P>',
            '<TABLE BORDER=0 CELLSPACING="2" CELLPADDING="4">',
           );
   ($dummy = $param->{'source_filename'}) =~ s/(?:\.[A-Za-z0-9_]+)\Z//;                    # get rid of the file extension ".xxxxx"
   foreach my $file (sort(@dir)) {
      CASE: {                                                                             # in this CASE-construct explanations are added to each generated file
         if (($file eq $param->{'source_filename'})) {
            push(@html, '   <TR bgcolor="'.$htmlcolor.'">',
                        '      <TD><font face="Arial, Helvetica, sans-serif" size="2">',
                        '         <A HREF="'.$file.'">'.$file.'</A>',
                        '      </TD>',
                        '      <TD><font face="Arial, Helvetica, sans-serif" size="2">',
                        '         Source alignment as it was provided by the user',
                        '      </TD>',
                        '   </TR>',
                );
            last CASE;
         }
         if (($file eq $param->{'source_filename'}.'_aligned')) {
            push(@html, '   <TR bgcolor="'.$htmlcolor.'">',
                        '      <TD><font face="Arial, Helvetica, sans-serif" size="2">',
                        '         <A HREF="'.$file.'">'.$file.'</A>',
                        '      </TD>',
                        '      <TD><font face="Arial, Helvetica, sans-serif" size="2">',
                        '         Source alignment as it was provided by the user, aligned',
                        '      </TD>',
                        '   </TR>',
                );
            last CASE;
         }
         if (($file =~ m/\.tar\.gz\Z/) || ($file =~ m/\.tgz\Z/)) {
            push(@html, '   <TR bgcolor="'.$htmlcolor.'">',
                        '      <TD><font face="Arial, Helvetica, sans-serif" size="2">',
                        '         <A HREF="'.$file.'">'.$file.'</A>',
                        '      </TD>',
                        '      <TD><font face="Arial, Helvetica, sans-serif" size="2">',
                        '         Your results and all intermediate data (like aligned data subsets etc.) as a GZIPed TAR-archive',
                        '      </TD>',
                        '   </TR>',
                );
            last CASE;
         }
         if ($file eq 'consensi_aligned.html') {
            push(@html, '   <TR bgcolor="'.$htmlcolor.'">',
                        '      <TD><font face="Arial, Helvetica, sans-serif" size="2">',
                        '         <A HREF="'.$file.'">'.$file.'</A>',
                        '      </TD>',
                        '      <TD><font face="Arial, Helvetica, sans-serif" size="2">',
                        '         <B>The aligned consensi, colorized according to their column-wise conservation rate</B>',
                        '      </TD>',
                        '   </TR>',
                );
            last CASE;
         }
         if ($file eq 'consensi_aligned.fasta') {
            push(@html, '   <TR bgcolor="'.$htmlcolor.'">',
                        '      <TD><font face="Arial, Helvetica, sans-serif" size="2">',
                        '         <A HREF="'.$file.'">'.$file.'</A>',
                        '      </TD>',
                        '      <TD><font face="Arial, Helvetica, sans-serif" size="2">',
                        '         The aligned consensi in FASTA format',
                        '      </TD>',
                        '   </TR>',
                );
            last CASE;
         }
         if ($file eq 'consensi_unaligned.fasta') {
            push(@html, '   <TR bgcolor="'.$htmlcolor.'">',
                        '      <TD><font face="Arial, Helvetica, sans-serif" size="2">',
                        '         <A HREF="'.$file.'">'.$file.'</A>',
                        '      </TD>',
                        '      <TD><font face="Arial, Helvetica, sans-serif" size="2">',
                        '         The <I>un</I>aligned consensi in FASTA format',
                        '      </TD>',
                        '   </TR>',
                );
            last CASE;
         }
         if (($file =~ m/\Agrp/) && ($file =~ m/\.fasta\Z/)) {
            push(@html, '   <TR bgcolor="'.$htmlcolor.'">',
                        '      <TD><font face="Arial, Helvetica, sans-serif" size="2">',
                        '         <A HREF="'.$file.'">'.$file.'</A>',
                        '      </TD>',
                        '      <TD><font face="Arial, Helvetica, sans-serif" size="2">',
                        '         A subset of "<I>'.$param->{'source_filename'}.'</I>", unaligned <B>or</B> not re-aligned',
                        '      </TD>',
                        '   </TR>',
                );
            last CASE;
         }
         if (($file =~ m/\Agrp/) && ($file =~ m/\.fasta_aligned\Z/)) {
            push(@html, '   <TR bgcolor="'.$htmlcolor.'">',
                        '      <TD><font face="Arial, Helvetica, sans-serif" size="2">',
                        '         <A HREF="'.$file.'">'.$file.'</A>',
                        '      </TD>',
                        '      <TD><font face="Arial, Helvetica, sans-serif" size="2">',
                        '         A subset of "<I>'.$param->{'source_filename'}.'</I>", aligned/re-aligned',
                        '      </TD>',
                        '   </TR>',
                );
            last CASE;
         }
         if (($file =~ m/\Aconsensus_grp/) && ($file =~ m/\.fasta\Z/)) {
            push(@html, '   <TR bgcolor="'.$htmlcolor.'">',
                        '      <TD><font face="Arial, Helvetica, sans-serif" size="2">',
                        '         <A HREF="'.$file.'">'.$file.'</A>',
                        '      </TD>',
                        '      <TD><font face="Arial, Helvetica, sans-serif" size="2">',
                        '         The consensus of a subset of "<I>'.$param->{'source_filename'}.'</I>" in FASTA format',
                        '      </TD>',
                        '   </TR>',
                );
            last CASE;
         }
         if (($file =~ m/\Aconsensus_grp/) && ($file =~ m/\.html\Z/)) {
            push(@html, '   <TR bgcolor="'.$htmlcolor.'">',
                        '      <TD><font face="Arial, Helvetica, sans-serif" size="2">',
                        '         <A HREF="'.$file.'">'.$file.'</A>',
                        '      </TD>',
                        '      <TD><font face="Arial, Helvetica, sans-serif" size="2">',
                        '         <B>The consensus of a subset of "<I>'.$param->{'source_filename'}.'</I>" in colorized HTML format</B>',
                        '      </TD>',
                        '   </TR>',
                );
            last CASE;
         }
         if ($file eq 'consensus_'.$dummy.'.fasta') {
            push(@html, '   <TR bgcolor="'.$htmlcolor.'">',
                        '      <TD><font face="Arial, Helvetica, sans-serif" size="2">',
                        '         <A HREF="'.$file.'">'.$file.'</A>',
                        '      </TD>',
                        '      <TD><font face="Arial, Helvetica, sans-serif" size="2">',
                        '         <B>The consensus of the source dataset "<I>'.$param->{'source_filename'}.'</I>" in FASTA format</B>',
                        '      </TD>',
                        '   </TR>',
                );
            last CASE;
         }
         if ($file eq 'consensus_'.$dummy.'.html') {
            push(@html, '   <TR bgcolor="'.$htmlcolor.'">',
                        '      <TD><font face="Arial, Helvetica, sans-serif" size="2">',
                        '         <A HREF="'.$file.'">'.$file.'</A>',
                        '      </TD>',
                        '      <TD><font face="Arial, Helvetica, sans-serif" size="2">',
                        '         <B>The consensus of the source dataset "<I>'.$param->{'source_filename'}.'</I>" in colorized HTML format</B>',
                        '      </TD>',
                        '   </TR>',
                );
            last CASE;
         }
         if (($file =~ m/\Aconsensus_$dummy/) && ($file =~ m/_simplified_[A-Za-z0-9]+\.fasta/)) {
            push(@html, '   <TR bgcolor="'.$htmlcolor.'">',
                        '      <TD><font face="Arial, Helvetica, sans-serif" size="2">',
                        '         <A HREF="'.$file.'">'.$file.'</A>',
                        '      </TD>',
                        '      <TD><font face="Arial, Helvetica, sans-serif" size="2">',
                        '         <B>The consensus of the simplified source dataset "<I>'.$param->{'source_filename'}.'</I>" in FASTA format</B>',
                        '      </TD>',
                        '   </TR>',
                );
            last CASE;
         }
         if (($file =~ m/\Aconsensus_$dummy/) && ($file =~ m/_simplified_[A-Za-z0-9]+\.html\Z/)) {
            push(@html, '   <TR bgcolor="'.$htmlcolor.'">',
                        '      <TD><font face="Arial, Helvetica, sans-serif" size="2">',
                        '         <A HREF="'.$file.'">'.$file.'</A>',
                        '      </TD>',
                        '      <TD><font face="Arial, Helvetica, sans-serif" size="2">',
                        '         <B>The consensus of the simplified source dataset "<I>'.$param->{'source_filename'}.'</I>" in colorized HTML format</B>',
                        '      </TD>',
                        '   </TR>',
                );
            last CASE;
         }
         if ($file eq 'viscose.log') {
            push(@html, '   <TR bgcolor="'.$htmlcolor.'">',
                        '      <TD><font face="Arial, Helvetica, sans-serif" size="2">',
                        '         <A HREF="'.$file.'">'.$file.'</A>',
                        '      </TD>',
                        '      <TD><font face="Arial, Helvetica, sans-serif" size="2">',
                        '         <B><I>This is the logfile of VisCoSe. If you experience any errors with your output or data, please look here first for error messages or warnings concerning dataset integrity or violations of the FASTA format etc.</B></I>',
                        '      </TD>',
                        '   </TR>',
                );
            last CASE;
         }
      }
      if ($htmlcolor eq $color1) { $htmlcolor = $color2; }
      else { $htmlcolor = $color1; }
   }
   push(@html, '</TABLE>',
               '</P>',
               '<P>',
               '<TABLE BORDER=0 CELLSPACING="0" CELLPADDING="4" WIDTH="100%">',
               '   <TR bgcolor="#EEEEFF">',
               '      <TD><font face="Arial, Helvetica, sans-serif">',
               '         <B>Job results are ephemeral. They will probably be gone within the hour!</B>',
               '      </TD>',
               '   </TR>',
               '</TABLE>',
               '</P>',
               '</HTML>'
       );

   open(_html_, '>'.temp_dir.'/'.$param->{'dir'}.'/viscose_results.html') or die "Could not open 'viscose_results.html' for write access...please contact spitzem\@uni-muenster.de!\n";
   foreach my $line (@html) { print _html_ $line."\n"; }
   close(_html_);
}

sub tar_gzip_results {
   my $param               = $_[0];
   my ($log, $shellstring) = (0, 0);

   chdir(temp_dir);
   $shellstring = 'tar -czvf '.$param->{'dir'}.'.tgz '.$param->{'dir'};
   $log = `$shellstring`;
   $shellstring = 'mv '.$param->{'dir'}.'.tgz '.$param->{'dir'}.'/';
   $log = `$shellstring`;
   chdir(cgi_dir);
}

sub delete_files {
   my $param      = $_[0];
   my $log        = '';

   chdir(temp_dir.'/'.$param->{'dir'});
}

### main ##########################################################################################

if (debug) {                                                                                       # redirect STDERR output so it is visible in the browser
   open (STDERR, ">&STDOUT") or print "Could not redirect STDERR!\n";
}

$| = 1;
umask(0000);                                                                                       # set 'rw(x)' for all user, group, other

my $param = &parse_args(\@ARGV);
&delete_files($param);
&tar_gzip_results($param);
&make_links($param);

exit(0);
