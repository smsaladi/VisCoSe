VisCoSe Installation instructions
=================================

Contents of the archive:
------------------------
- install_instructions.txt
  This text
- viscose.tbz2
  Archive of the VisCoSe CGI scripts and required binaries (MAFFT 3.89)
- mafft-3.89-040205.tgz
  The source files of MAFFT 3.89 for custom compilation
  (Note: AFAIK VisCoSe does not work with current MAFFT versions due to
  changes in the command line interface)
- vhost_viscose
  Example VHOST configuration for Apache2.x

The archive already contains compiled MAFFT binaries of version 3.89
(compiled for i686 generic architecture, no specific optimizations apart
from -O2 and other basic GCC options). I you would like to use your own
compiled binaries, either use the MAFFT source within this distribution
or search for MAFFT 3.89 archives on the Internet. Good luck, though,
this version is hard to find.

VisCoSe was developed under Linux and never tested under MS Windows. No
support whatsoever is given for installation under MS Windows. Other
required software is:
- Any *nix system will do AFAIK (Solaris, Linux, *BSD, etc.)
- Apache 2.x
- mod_cgi (probably mod_perl, yet never tested)
- MAFFT V3.89
- Perl modules:
  - Bit::Vector
    (http://guest.engelschall.com/~sb/download/ or search CPAN)
  - CGI.pm
  - POSIX
  - File::Basename
  - Data::Dumper (not necessarily, but left within the script for
    debugging reasons)
- probably Sun Grid Engine, if you want to use job distribution features
  for clustered environments. Basic job queue submission is implemented
  in the VisCoSe CGI script.

I hope I didn't forget anything.



Installation
------------

1. Configure and activate a VHOST for VisCoSe, if needed and/or wanted.
You need to have at least basic knowledge about Apache configuration.
Be sure to correctly configure the ScriptAlias in the VHOST config.
The VHOST base dir must point to the "../viscose/" path.

Also, install all of the above mentioned required third party software
(like Bit::Vector etc.). I recommend using CPAN for installation of
Perl modules.


2. Copy "viscose.tbz2" into a suitable directory and untar the VisCoSe
archive:
shell> tar xjvf viscose.tbz2

- Make sure that all files are readable by the webserver.
- Make sure that the TEMP directory where VisCoSe is supposed to store
  it's temporary files and results (see next topic "3.") is WRITEABLE for
  the webserver!
- Make sure that the MAFFT binaries can be executed by the webserver.
- Make sure that the following Perl scripts are executable for the webserver:
  - interface.pl
  - make_result_page.pl
  - viscose.pl


3. Configure ./viscose/cgi-bin/interface.pl (sorry, parameters are hard-
coded in the Perl script):
-> Server root where the VisCoSe directory is located:
use constant server_root => '/var/www/vhosts/viscose';

-> Path to TEMP dir that VisCoSe should use:
use constant temp_dir    => '/var/www/vhosts/viscose/temp';

-> Path to the CGI directory (for locating config files. alphabets, etc.)
use constant cgi_dir     => '/var/www/vhosts/viscose/cgi-bin';

-> Path to the MAFFT binary directory:
use constant mafft_bin   => '/var/www/vhosts/viscose/cgi-bin/mafft';

-> URL that is used to link (via HTML!) the VisCoSe results to. This is
essential, otherwise users won't be able to view and download their results.
use constant link_url    => 'http://viscose.host.de/temp';


4. Configure "make_result_page.pl". Same parameters as for "interface.pl",
see above topic "3."


5. Configure "viscose.pl":
-> Add path of the VisCoSe modules to the Perl module path INC (set this path
to the correct location of YOUR installation):
use lib '/var/www/vhosts/viscose/cgi-bin';

NOTE1: the script "status.pl" is COMPLETELY unsupported! I used this myself as
a basic, browser-accessible, monitor for the Sun GridEngine environment to
watch running VisCoSe jobs online. USE AT YOUR OWN RISK!

NOTE2: this software is used at your own risk! I do not take responsibility
in loss of data or any other harmful effects that may be caused by the usage
of this software!


===================================================
Jan 2008, Michael Spitzer (spitzem@uni-muenster.de)