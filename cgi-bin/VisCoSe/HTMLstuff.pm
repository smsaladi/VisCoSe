###################################################################################################
# HTMLstuff.pm
# (c) 2002 Michael Spitzer, IFG Muenster (Integrated Funtional Genomics), Germany
#
# 20.03.2003   added format_consensus_sequence(), format_consensus_comparison()
# 18.03.2003   fixed and extended some things for use with multiple datasets
# 14.03.2003   initial version, bugs in _prepare_html_color() fixed, y-axis now next to graph

package VisCoSe::HTMLstuff;
require Exporter;

use strict;
use diagnostics;
use Data::Dumper;

use constant range0           => 25;                       
use constant range1           => 50;                       
use constant range2           => 75;                       
use constant range3           => 100;                      
use constant color0           => '<font color="#00xxFF">'; 
use constant color1           => '<font color="#00FFxx">'; 
use constant color2           => '<font color="#xxFF00">'; 
use constant color3           => '<font color="#FFxx00">'; 
use constant offset           => 8;
use constant minscaleoffset   => 20;
use constant def_alphabet     => '20 amino acids';

our @EXPORT  = qw(prep_header         prepare_single_html_consensus       format_consensus_comparison    get_max_defline_length
                  split_html_file);

our $VERSION = "18.03.2003";

sub sort_num { $a <=> $b }

sub prep_header {
   my $filename      = $_[0];
   my $param         = $_[1];
   my @header = ('<html>',
                 "<head><title>Consensus sequence for alignment '".$filename."'</title></head>",
                 '<basefont size="5" color="'.$param->{'fg'}.'" face="Courier">',
                 '<body bgcolor="'.$param->{'bg'}.'" text="'.$param->{'fg'}.'">',
                 '<H3>Consensus sequence for alignment <a href="'.$filename.'">"'.$filename.'"</a>:</H3>',
                 '<pre font-family:Fixedsys,Courier; padding:10px; style="font-size:'.$param->{'fs'}.'pt">',
                );
   my @footer = ('</pre>',
                 '</body>',
                 '</html>',
                );
   return (\@header, \@footer);
}

sub _prepare_html_color {
   my $value                                 = $_[0];
   my $param                                 = $_[1];
   my ($color, $hexcolor, $htmlcolor)        = (0, 0, 0);

   CASE: { 
      if ($value <= range0) {                                                                   # $value is 0-25
         $color     = $value * (255 / range0);                                                  # 25(%) * 10,2 = max. 255 as colorcode / increase green
         $hexcolor  = uc(sprintf("%lx", $color));
         if (length($hexcolor) == 1) { $hexcolor = "0".$hexcolor; }                             # if the HEX number has only length of 1, put a zero in front of it
         ($htmlcolor = color0) =~ s/xx/$hexcolor/;
         last CASE;
      }
      if ($value <= range1) {                                                                   # $value is 26-50
         $color     = (((range1- $value) / 25) * 255);                                          # 75(%) * 3,4 = max. 255 as colorcode / increase red
         $hexcolor  = uc(sprintf("%lx", $color));
         if (length($hexcolor) == 1) { $hexcolor = "0".$hexcolor; }                             # if the HEX number has only length of 1, put a zero in front of it
         ($htmlcolor = color1) =~ s/xx/$hexcolor/;
         last CASE;
      }
      if ($value <= range2) {                                                                   # $value is 51-75
         $color     = 255 - (((range2 - $value) / 25) * 255);                                   # 75(%) * 3,4 = max. 255 as colorcode / increase red
         $hexcolor  = uc(sprintf("%lx", $color));
         if (length($hexcolor) == 1) { $hexcolor = "0".$hexcolor; }                             # if the HEX number has only length of 1, put a zero in front of it
         ($htmlcolor = color2) =~ s/xx/$hexcolor/;
         last CASE;
      }
      if ($value <= range3) {                                                                   # $value is 76-100
         $color     = (((range3 - $value) / 25) * 255);                                         # 75(%) * 3,4 = max. 255 as colorcode / increase red
         $hexcolor  = uc(sprintf("%lx", $color));
         if (length($hexcolor) == 1) { $hexcolor = "0".$hexcolor; }                             # if the HEX number has only length of 1, put a zero in front of it
         ($htmlcolor = color3) =~ s/xx/$hexcolor/;
         last CASE;
      }
   }
   return $htmlcolor;
}

sub _prepare_single_graph_column {
   my $graph                           = $_[0];                                                                   # reference to array '@graph'
   my $column                          = $_[1];                                                                   # the actual column we prepare a graph for
   my $rate                            = $_[2];                                                                   # conservation rate
   my $aa                              = $_[3];                                                                   # aminoacid-symbol at actual position
   my $param                           = $_[4];
   my $htmlcolor                       = $_[5];                                                                   # if user wants a colored graph, here is the color
   my $height                          = IFG::NutsnBolts::round(($rate / 100) * $param->{'graph_height'}, 0);     # calculate height of the actual graph column

   for (my $y = ($param->{'graph_height'} - 1); $y > -1; $y--) {
      if ($aa ne '-') {
         ($y <= $height) ? ($graph->[$column]->[$y] = 1) : ($graph->[$column]->[$y] = -1);                        # mark the positions occupied by the graph column with '1', the empty space with '-1'
      }
      else {
         $graph->[$column]->[$y] = -1;
      }
   }
   if ($htmlcolor) { $graph->[$column]->[$param->{'graph_height'}] = $htmlcolor; }                                # save the right color for that graph column, we need it later in 'prepare_html_graph()'
}

sub generate_y_axis {
   my $htmlgraph                       = $_[0];          # *reference* to @htmlgraph from '_prepare_html_graph()'
   my $param                           = $_[1];
   my ($step, $htmlcolor, $text)       = (0, 0, 0);

   for (my $y = 0; $y < $param->{'graph_height'}; $y++) {
      $step = &IFG::NutsnBolts::round((($y+1) * (100 / $param->{'graph_height'})), 2);
      $htmlcolor = &_prepare_html_color($step, $param);
      if ($step < 10)     { $text = '  '.$step; }
      elsif ($step < 100) { $text = ' '.$step; }
      else {$text = $step; }
      $htmlgraph->[$y] .= $htmlcolor.$text.'%</font> ';
   }
}

sub _prepare_html_graph {
   my $graph                           = $_[0];
   my $consensus                       = $_[1];
   my $maxcodelen                         = $_[2];
   my $param                           = $_[3];
   my (@htmlgraph, @axis)              = ((), ());
   my $break                           = 0;

   if ($maxcodelen < minscaleoffset) { $maxcodelen = minscaleoffset; }
   &generate_y_axis(\@htmlgraph, $param);                                                                                     # Careful! '@htmlgraph' will be modified in 'generate_y_axis()'
   foreach my $line (@htmlgraph) {                                                                                            # extract the axis information (my, that code got really crappy...
      $line =~ m/(<font.+<\/font>)/;                                                                                          # ...if I had the time, I'd rewrite it from scratch...)
      push(@axis, $1);
   }
   for (my $y = 0; $y < $param->{'graph_height'}; $y++) {                                                                     # step down the graph-height
      foreach my $column (sort(sort_num keys(%{$consensus}))) {                                                               # step thru each graph-column
         $break = 0;
         if ($param->{'colorgraph'}) {
            if ($param->{'linebreak'} && ($column > 0) && (($column % $param->{'linebreak'}) == 0)) {                         # mark the position where a linebreak shall be introduced later
               $htmlgraph[$y] .= '</font>'.'###break###'.(' ' x ($maxcodelen - 4)).$axis[$y].' ';
               $break = 1;
            }
            if ($column == 0) {                                                                                               # first column? Apply color if ON-state!
               if ($graph->[$column]->[$y] == -1) {                                                                           # both graph-columns have 'OFF' state
                  $htmlgraph[$y] .= ' ';                                                                                      # don't change color!
               }
               else {                                                                                                         # both graph-columns have different colors and 'ON' state
                  $htmlgraph[$y] .= $graph->[$column]->[$param->{'graph_height'}].$param->{'graph_char'};                     # change color!
               }
            }
            elsif ($graph->[$column]->[$y] == $graph->[$column-1]->[$y]) {                                                    # both graph-columns have the same STATE at height 'y'
               if ($graph->[$column]->[$y] == -1) {                                                                           # both graph-columns have 'OFF' state
                  $htmlgraph[$y] .= ' ';                                                                                      # don't change color!
               }
               elsif ($graph->[$column]->[$param->{'graph_height'}] eq $graph->[$column-1]->[$param->{'graph_height'}]) {     # both graph-columns have the same COLOR and 'ON' state
                  if ($break) { $htmlgraph[$y] .= $graph->[$column]->[$param->{'graph_height'}].$param->{'graph_char'}; }     # change color at linebreak!!!
                  else { $htmlgraph[$y] .= $param->{'graph_char'}; }                                                          # don't change color!
               }
               else {                                                                                                         # both graph-columns have different colors and 'ON' state
                  if ($break) { $htmlgraph[$y] .= $graph->[$column]->[$param->{'graph_height'}].$param->{'graph_char'}; }     # we don't need a </font>-tag at linebreaks!
                  else { $htmlgraph[$y] .= '</font>'.$graph->[$column]->[$param->{'graph_height'}].$param->{'graph_char'}; }  # change color!
               }
            }
            else {                                                                                                            # both graph-columns have different states at height 'y'
               if ($graph->[$column]->[$y] == -1) {                                                                           # 'y' has 'OFF' state, so 'y-1' must have 'ON' state!
                  $htmlgraph[$y] .= ' ';                                                                                      # don't change color! Not needed when changing from 'ON' to 'OFF'!!!
               }
               else {                                                                                                         # 'y' has 'ON' state, so 'y'-1' must have 'OFF' state!
                  if ($break) { $htmlgraph[$y] .= $graph->[$column]->[$param->{'graph_height'}].$param->{'graph_char'}; }     # we don't need a </font>-tag at linebreaks!
                  else { $htmlgraph[$y] .= '</font>'.$graph->[$column]->[$param->{'graph_height'}].$param->{'graph_char'}; }  # change color! Needed when changing from 'OFF' to 'ON'!!!
               }
            }
         }
         else {                                                                                                               # here we don't need to apply any color changes at all (since the graph is single-colored)
            if ($param->{'linebreak'} && ($column > 0) && (($column % $param->{'linebreak'}) == 0)) {                         # mark the position where a linebreak shall be introduced later
               $htmlgraph[$y] .= '###break###'.(' ' x ($maxcodelen - 4)).$axis[$y].' ';
               #if (($param->{'sa'} || $param->{'ca'} || $param->{'fa'}) && ($maxcodelen > offset)) {
               #   $htmlgraph[$y] .= (' ' x ($maxcodelen - 4));                                                               # if defline length is >offset: fill up with spaces. '5' is the minimal offset due to the y-axis!
               #}
               #$htmlgraph[$y] .= $axis[$y].' ';
               $break = 1;
            }
            if ($graph->[$column]->[$y] == -1) { $htmlgraph[$y] .= ' '; }                                                     # 'y' has 'OFF' state
            else { $htmlgraph[$y] .= $param->{'graph_char'}; }                                                                # 'y' has 'ON' state
         }
      }
      $htmlgraph[$y] .= '</font>';
   }
   return \@htmlgraph;
}

sub _prepare_scale_axis {     # prepares the x-axis with the numbering of the alignment/consensus
   my $len                             = $_[0];
   my $maxcodelen                      = $_[1];
   my $consensi                        = $_[2];    # used to determine if we're doing axis for the consensus comparison or the alignment visualization
   my $param                           = $_[3];
   my (@axis, @splitaxis, @splitpos)   = ((), (), ());
   my ($dummy, $start, $times)         = ('', 0, 1);

   if ($maxcodelen < minscaleoffset) { $maxcodelen = minscaleoffset; }
   push(@splitpos, 0);                                                                                               # push initial (starting) position into array with split-positions
   for (my $i = 0; $i < $len; $i++) {                                                                                # step thru the complete length (this is the actual sequence data here in "$len"!!!!)
      if ($param->{'linebreak'} && ($i > 0) && (($i % $param->{'linebreak'}) == 0)) {                                # mark the position where a linebreak shall be introduced later
         push(@splitpos, $i);                                                                                        # remember positions where to split the scale axis
      }
      if (($i % 20) == 0) {                                                                                          # are we at a 20th*i position?
         $axis[0] .= '|';
         $axis[1] .= $i;
      }
      elsif (($i % 10) == 0) {                                                                                       # or are we at a 10th*i position?
         $axis[0] .= ':';
         if (length($axis[1]) <= $i) { $axis[1] .= ' '; }                                                            # only add spaces if the linelengths are equal!
      }
      else {                                                                                                         # any other position, not 20th*i nor 10th*i
         $axis[0] .= '.';
         if (length($axis[1]) <= $i) { $axis[1] .= ' '; }                                                            # only add spaces if the linelengths are equal!
      }
   }
   push(@splitpos, $len);                                                                                            # push last (ending) position into array with split-positions
   $dummy = (' ' x ($maxcodelen + 4));                                                                               # '+4' needed for correct indention of consensus sequence
   for (my $i = 1; $i < scalar(@splitpos); $i++) {
      if ($i == 1) {
         $splitaxis[0] .= $dummy.substr($axis[0], $splitpos[$i-1], ($splitpos[$i] - $splitpos[$i-1]));
         $splitaxis[1] .= $dummy.substr($axis[1], $splitpos[$i-1], ($splitpos[$i] - $splitpos[$i-1]));
      }
      else {
         $splitaxis[0] .= '###break###'.$dummy.substr($axis[0], $splitpos[$i-1], ($splitpos[$i] - $splitpos[$i-1]));
         $splitaxis[1] .= '###break###'.$dummy.substr($axis[1], $splitpos[$i-1], ($splitpos[$i] - $splitpos[$i-1]));
      }
   }
   return \@splitaxis;
}

sub get_max_defline_length {
   my $codes      = $_[0];
   my $param      = $_[1];
   my $codelen    = 0;

   for (my $i = 0; $i < scalar(@{$codes}); $i++) {                                                                   # get length of the longest defline
      if (length($codes->[$i]) > $codelen) { $codelen = length($codes->[$i]); }
   }
   return $codelen;
}

sub _prepare_alignment {
   my $aln                                = $_[0];                # alignment-hash
   my $codes                              = $_[1];                # deflines-array
   my $consensus                          = $_[2];                # the consensus-hash
   my $graph                              = $_[3];                # the graph (we need only the color "$graph->[$column]->[$param->{'graph_height'}]" here!)
   my $conserve                           = $_[4];                # conservation rates: {column} -> {character} -> {rate}
   my $param                              = $_[5];
   my ($dummy, $aa, $seq)                 = (0, 0, 0);
   my $maxcodelen                         = get_max_defline_length($codes, $param);
   my $conss_before                       = 0;                    # if defined: indicates that the *previous* position was consensus-colored
   my ($htmlcolor, $precolor, $break)     = (0, 0, 0);
   my @alignment                          = ();

   if ($maxcodelen < minscaleoffset) { $maxcodelen = minscaleoffset; }
   foreach my $code (@{$codes}) {                                                                                    # step thru all sequences of the alignment
      if (length($code) < $maxcodelen) {
         $dummy = $code.(' ' x ($maxcodelen - length($code)));                                                       # eventually fill the defline up with whitespaces if it's to short
      }
      else { $dummy = $code; }
      $seq = $dummy.':   ';
      if ($param->{'ca'} || $param->{'fa'}) {                                                                        # prepare a colorcoded alignment
         $conss_before = 0;
         foreach my $column (sort(sort_num keys(%{$consensus}))) {                                                   # step thru all columns of the consensus (=columns of alignment!)
            $break = 0;
            if ($param->{'linebreak'} && ($column > 0) && (($column % $param->{'linebreak'}) == 0)) {                # mark the position where a linebreak shall be introduced later
               $seq .= '</font>'.'###break###'.$dummy.':   ';
               $break = 1;
            }
            $aa = substr($aln->{$code}, $column, 1);                                                                 # get the amino acid character at position '$column'
            if ($param->{'fa'}) {                                                                                    # here we handle a "full colored" alignment
               if ($aa ne '-') {
                  CASE_FA:{
                     $htmlcolor = &_prepare_html_color($conserve->{$column}->{$aa}, $param);
                     if ($column == 0) {                                                                             # at position 0, we always need a initial '<font...>'-tag
                        $seq .= $htmlcolor.$aa;
                        $precolor = $htmlcolor;
                        last CASE_FA;
                     }
                     if (($htmlcolor eq $precolor) && (!($break))) {                                                 # if the colors of the actual and the previous position...
                        ($conss_before) ? ($seq .= $aa) : ($seq .= $htmlcolor.$aa);                                  # ...are the same we don't need to set a new '<font>'-tag! (saves some tags!)
                        $precolor = $htmlcolor;
                        last CASE_FA;
                     }
                     if (($htmlcolor ne $precolor) || $break) {                                                      # if the colors of the actual and the previous position...
                        ($conss_before) ? ($seq .= '</font>'.$htmlcolor.$aa) : ($seq .= $htmlcolor.$aa);             # ...are *NOT* the same we need to set a new '<font>'-tag, but...
                        $precolor = $htmlcolor;                                                                      # ...*NOT* if we had a gap-character *before*! (saves some tags!)
                        last CASE_FA;
                     }
                     print ">>> [$column] $aa\t<- !ERROR!\n";
                  }
                  $conss_before = 1;                                                                                 # define $conss_before to remember we had an AA here
               }
               else {
                  if ($conss_before) {
                     $seq .= '</font>'.$aa;                                                                          # only close the '<font>'-tag if we tagged the previous position!
                     $conss_before = 0;                                                                              # undefine $conss_before, because we tagged a gap now!
                  }
                  else { $seq .= $aa; }                                                                              # if $conss_before is not defined, we tagged a gap before! Therefore no need to close the FONT-tag here!
               }
            }
            else {                                                                                                   # here we handle a "consensus colored" alignment
               if (($aa eq $consensus->{$column}->{'symbol'}) && ($aa ne '-')) {                                     # is the actual aminoacid the same as in the consensus and not a gap?
                  $htmlcolor = $graph->[$column]->[$param->{'graph_height'}];                                        # HTML-color at the actual position
                  $precolor  = $graph->[$column-1]->[$param->{'graph_height'}];                                      # HTML-color at the previous position
                  CASE_CA:{
                     if ($column == 0) {                                                                             # at position 0, we always need a initial '<font...>'-tag
                        $seq .= $htmlcolor.$aa;
                        last CASE_CA;
                     }
                     if ($htmlcolor eq $precolor) {                                                                  # if the colors of the actual and the previous position...
                        ($conss_before) ? ($seq .= $aa) : ($seq .= $htmlcolor.$aa);                                  # ...are the same we don't need to set a new '<font>'-tag! ...
                        last CASE_CA;
                     }
                     if ($htmlcolor ne $precolor) {                                                                  # if the colors of the actual and the previous position...
                        ($conss_before) ? ($seq .= '</font>'.$htmlcolor.$aa) : ($seq .= $htmlcolor.$aa);             # are *NOT* the same, we need to set a new color, but we need to close the...
                        last CASE_CA;
                     }
                  }
                  $conss_before  = 1;
               }
               else {
                  if ($conss_before) {                                                                               # if we had an AA before we need to close the FONT-tag...
                     $seq .= '</font>'.$aa;
                     $conss_before = 0;
                  }
                  else { $seq .= $aa; }                                                                              # ...otherwise (gap before) we can continue to append gaps
               }
            }
         }
         push(@alignment, $seq.'</font>');
      }
      else {                                                                                                         # prepare a plain, non-colored alignment
         foreach my $column (sort(sort_num keys(%{$consensus}))) {                                                   # step thru all columns of the consensus (=columns of alignment!)
            $break = 0;
            if ($param->{'linebreak'} && ($column > 0) && (($column % $param->{'linebreak'}) == 0)) {                # mark the position where a linebreak shall be introduced later
               $seq .= '###break###'.$dummy.':   ';
               $break = 1;
            }
            $seq .= substr($aln->{$code}, $column, 1);                                                               # get the amino acid character at position '$column'
         }
         push(@alignment, $seq);
      }
   }
   push(@alignment, ' ');
   return \@alignment;
}

sub prepare_single_html_consensus {
   my $consensus                       = $_[0];
   my $aln                             = $_[1];
   my $codes                           = $_[2];
   my $filename                        = $_[3];
   my $conserve                        = $_[4];                         # conservation rates: {column} -> {character} -> {rate}
   my $prevconss                       = $_[5];                         # previously calculated 20AA consensus (if defined)
   my $alphabet_name                   = $_[6];                         # alphabet name (defined if $prevconss is defined!!!)
   my $param                           = $_[7];
   my (@file, @graph)                  = ((), ());                      # '@graph' is modified in prepare_single_graph_column()
   my ($seq, $axis, $alignment)        = ("", "");                      # '$seq' will hold the HTML-formatted consensus sequence
   my ($htmlcolor, $precolor, $len)    = (0, 0, 0);
   my $maxcodelen                      = get_max_defline_length($codes, $param);
   my ($break, $breaknum)              = (0, 0);
   my $def_alphabet_var                = def_alphabet;

   if ($maxcodelen < minscaleoffset) { $maxcodelen = minscaleoffset; }
   if (($param->{'sa'} || $param->{'ca'} || $param->{'fa'}) && (!($prevconss))) {
      if (length(def_alphabet) <= $maxcodelen) {
         $def_alphabet_var .= (' ' x ($maxcodelen - length(def_alphabet)));                                          # generate consensus name for default consensus
      }
      $seq = $def_alphabet_var.':   ';
   }
   elsif (($param->{'sa'} || $param->{'ca'} || $param->{'fa'}) && $prevconss) {
      if (length($alphabet_name) <= $maxcodelen) {
         $alphabet_name .= (' ' x ($maxcodelen - length($alphabet_name)));                                           # generate consensus name for alphabet consensus
      }
      if (length($alphabet_name) < length(def_alphabet)) {
         $alphabet_name .= (' ' x (length(def_alphabet) - length($alphabet_name)));
      }
      $seq = $alphabet_name.':   ';
   }
   elsif ($prevconss) {
      if (length($alphabet_name) <= length(def_alphabet)) {
         $alphabet_name .= (' ' x (length(def_alphabet) - length($alphabet_name)));                                  # generate consensus name for alphabet consensus
      }
      else {
         $alphabet_name = '...'.substr($alphabet_name, (length($alphabet_name) - (length(def_alphabet) - 3)), (length(def_alphabet) - 3));
      }
      $seq = $alphabet_name.':   ';
   }
   else { $seq = def_alphabet.':   '; }
   foreach my $column (sort(sort_num keys(%{$consensus}))) {                                                         # prepare colorcoded consensus and store in '$seq'
      $break = 0;
      if ($param->{'linebreak'} && ($column > 0) && (($column % $param->{'linebreak'}) == 0)) {                      # mark the position where a linebreak shall be introduced later
         $seq .= '</font>'.'###break###';
         if (($param->{'sa'} || $param->{'ca'} || $param->{'fa'}) && (!($prevconss))) {
            $seq .= $def_alphabet_var.':   ';
         }
         elsif (($param->{'sa'} || $param->{'ca'} || $param->{'fa'}) && (!($prevconss))) {
            $seq .= $alphabet_name.':   ';
         }
         elsif ($prevconss) {
            $seq .= $alphabet_name.':   ';
         }
         else {
            $seq .= $def_alphabet_var.':   ';
         }
         $breaknum++;
         $break = 1;
      }
      $htmlcolor = &_prepare_html_color($consensus->{$column}->{'rate'}, $param);
      if ($consensus->{$column}->{'symbol'} eq '-') { $htmlcolor = '<font color="'.$param->{'fg'}.'">'; }            # if a gap occurs in the consensus sequence: paint it black!!!
      if (($htmlcolor eq $precolor) && (!($break))) {                                                                # if the actual position is the same color as the preceding...
         $seq .= $consensus->{$column}->{'symbol'};                                                                  # ...we don't need to change the color!
         &_prepare_single_graph_column(\@graph, $column, $consensus->{$column}->{'rate'}, $consensus->{$column}->{'symbol'}, $param, $htmlcolor); # prepare a single column at '$column' of the graph
      }
      else {                                                                                                         # if the actual position has another color than the preceding...
         $seq .= '</font>'.$htmlcolor.$consensus->{$column}->{'symbol'};                                             # ...close the <font...>-tag and add a new one!
         $precolor = $htmlcolor;
         &_prepare_single_graph_column(\@graph, $column, $consensus->{$column}->{'rate'}, $consensus->{$column}->{'symbol'}, $param, $htmlcolor); # prepare a single column at '$column' of the graph
      }
      $len++;                                                                                                        # count sequence length, needed later
   }
   if ($param->{'sa'} || $param->{'ca'} || $param->{'fa'}) {                                                         # add alignment to HTML file?
      if ($maxcodelen > length(def_alphabet)) { $param->{'defline_len'} = $maxcodelen; }                             # !!! deprecated !!! SUPERFLUOUS !!! NEEDS A REWRITE !!!!!!! geeeeez...........
      else { $param->{'defline_len'} = length(def_alphabet); }                                                       # no alignment shown: set defline_length to a default value
      $axis = &_prepare_scale_axis($len, $maxcodelen, 0, $param);                                                    # generate x-axis
      for (my $i = (scalar(@{$axis}) - 1); $i > -1; $i--) { push(@file, $axis->[$i]); }                              # add the scale-text to the file-array
      $alignment = &_prepare_alignment($aln, $codes, $consensus, \@graph, $conserve, $param);
      push(@file, @{$alignment});                                                                                    # add prepared alignment-HTML-code to '@file'
      for (my $i = (scalar(@{$axis}) - 1); $i > -1; $i--) { push(@file, $axis->[$i]); }                              # add the scale-text to the file-array
      if ($prevconss) {
         push(@file, $prevconss.'</font>');
         push(@file, $axis->[0]);
      }
      push(@file, $seq.'</font>');
      my $htmlgraph = &_prepare_html_graph(\@graph, $consensus, $maxcodelen, $param);
      &_append_graph(\@file, $htmlgraph, $maxcodelen, $param);
      for (my $i = 0; $i < scalar(@{$axis}); $i++) {
         push(@file, '<font color="'.$param->{'fg'}.'">'.$axis->[$i].'</font>');                                     # add the scale-text to the file-array
      }
   }
   else {
      $param->{'defline_len'} = length(def_alphabet);
      $axis = &_prepare_scale_axis($len, $maxcodelen, 0, $param);
      for (my $i = (scalar(@{$axis}) - 1); $i > -1; $i--) { push(@file, $axis->[$i]); }                              # add the scale-text to the file-array
      if ($prevconss) {
         push(@file, $prevconss.'</font>');
         push(@file, $axis->[0]);
      }
      push(@file, $seq.'</font>');
      my $htmlgraph = &_prepare_html_graph(\@graph, $consensus, $maxcodelen, $param);
      &_append_graph(\@file, $htmlgraph, $maxcodelen, $param);
      for (my $i = 0; $i < scalar(@{$axis}); $i++) {
         push(@file, '<font color="'.$param->{'fg'}.'">'.$axis->[$i].'</font>');                                     # add the scale-text to the file-array
      }
   }
   $param->{'breaknum'} = $breaknum;                                                                                 # remember the number of linebreaks that were introduced
   return (\@file, $seq);
}

sub _append_graph {
   my $file                            = $_[0];       # the array holding the complete HTML-file
   my $graph                           = $_[1];       # the graph in HTML code
   my $prepend                         = $_[2];       # length of graph/consensus sequence
   my $param                           = $_[3];

   for (my $y = 0; $y < $param->{'graph_height'}; $y++) {
      if ($prepend) { push(@{$file}, (' ' x ($prepend - 4)).$graph->[$y]); }
      else { push(@{$file}, $graph->[$y]); }
   }
   $file->[-1] .= '</font>';
}

sub format_consensus_sequence {
   my $aln_cons                           = $_[0];                # consensus sequence (eventually re-aligned by MAFFT)
   my $consensus                          = $_[1];                # consensus of the corresponding consensus, along with conservation-rates and stuff
   my $defline                            = $_[2];
   my $codelen                            = $_[3];
   my $param                              = $_[4];
   my ($pos, $aa, $pre_gap, $seq, $break) = (0, "", "", "", 0);   # $pos is used to count the non-gap-positions and to get the corresponding conservation rate from $consensus
   my ($precolor, $htmlcolor, $breaknum)  = ("", "", 0);

   for (my $i = 0; $i < length($aln_cons); $i++) {
      if ($param->{'linebreak'} && ($i > 0) && (($i % $param->{'linebreak'}) == 0)) {            # mark the position where a linebreak shall be introduced later
         $seq .= '</font>'.'###break###'.$defline.':   ';
         $breaknum++;
         $break = 1;
      }
      $aa = substr($aln_cons, $i, 1);
      if ($aa ne '-') {                                                                         # ignore gap positions
         $pre_gap = 0;                                                                          # remember we have *NOT* a gap now but an aminoacid
         if ($aa eq $consensus->{$pos}->{'symbol'}) {
            $htmlcolor = &_prepare_html_color($consensus->{$pos}->{'rate'}, $param);
            if (($htmlcolor eq $precolor) && (!($break))) {
               $seq .= $consensus->{$pos}->{'symbol'};
            }
            else {
               $seq .= '</font>'.$htmlcolor.$consensus->{$pos}->{'symbol'};
               $precolor = $htmlcolor;
            }
         }
         else {
            if ($param->{'fra'}) { print "[$pos] bad AA ; "; }
         }
         unless ($param->{'aligned_source'}) { $pos++; }                                        # increase the "aminoacid-counter" because we have a "no-gap" here!
      }
      else {
         $precolor = 'a'; $htmlcolor = 'b';                                                     # after we had a gap we *MUST* apply a new color at the next first aminoacid
         if ($pre_gap) { $seq .= '-'; }                                                         # if we had a gap *BEFORE* we don't need to apply a '</font>' again
         else { $seq .= '</font>-'; }                                                           # if we didn't have a gap before we must close the '<font...>'-tag to prevent gaps colored in consensus-colors
         $pre_gap = 1;                                                                          # remember we *HAVE* a gap and *NO* aminoacid
      }
      if ($param->{'aligned_source'}) { $pos++; }                                               # increase the "aminoacid-counter" at this point because we did not realign and we generated the consensi from a pre-aligned source!
   }
   $param->{'breaknum'} = $breaknum;                                                            # remember the number of linebreaks that were introduced
   return $seq;
}

sub format_consensus_comparison {
   my $allconss                  = $_[0];
   my $allcodes                  = $_[1];
   my $consensus                 = $_[2];          # keys of hash-ref $consensus are equal to member of array-ref $allcodes/keys of hash-ref $allconss
   my $param                     = $_[3];
   my ($header, $footer)         = &prep_header('Consensus comparison', $param);
   my @htmlfile                  = ();
   my ($defline, $seq)           = (0, "");
   my $codelen                   = get_max_defline_length($allcodes, $param);
   my ($axis, $filename)         = (0, "");

   if ($codelen < minscaleoffset) { $codelen = minscaleoffset; }
   $axis = &_prepare_scale_axis(length($allconss->{$allcodes->[0]}), $codelen, 1, $param);
   for (my $i = (scalar(@{$axis}) - 1); $i > -1; $i--) {
      push(@htmlfile, $axis->[$i]);                                                                         # add the scale-text to the file-array
   }
   foreach my $file (sort(@{$allcodes})) {                                                                  # sort by filename / groups
      if (length($file) < $codelen) { $defline = $file.(' ' x ($codelen - length($file))); }                # eventually fill the defline up with whitespaces if it's to short
      else { $defline = $file; }
      $filename = &VisCoSe::ConsensusTools::generate_HTML_consensus_filename($file);
      $defline = '<a href="'.$filename.'.html">'.$defline.'</a>';                                           # insert a link to the single dataset consensus-report
      if ($param->{'verbose'}) { print "   $file..."; }
      $seq = &format_consensus_sequence($allconss->{$file}, $consensus->{$file}, $defline, $codelen, $param);
      push(@htmlfile, $defline.':   '.$seq.'</font>');                                                      # add the formatted sequence to the HTMLfile-array
      if ($param->{'verbose'}) { print "done\n"; }
   }
   for (my $i = 0; $i < scalar(@{$axis}); $i++) {
      push(@htmlfile, '<font color="'.$param->{'fg'}.'">'.$axis->[$i].'</font>');                           # add the scale-text to the file-array
   }
   return \@htmlfile;
}

sub split_html_file {
   my $html                         = $_[0];                   # the alignment as "final" HTML code
   my $maxcodelen                   = $_[1];
   my $param                        = $_[2];
   my $numlines                     = scalar(@{$html});        # number of lines of the HTML file
   my (@splitlines, @splitlines2)   = ((), ());                # we will store the splitted lines here as a two-dim array
   my @splithtml                    = ();                      # new HTML file with linebreaks

   for (my $i = 0; $i < $numlines; $i++) {
      if ($html->[$i] =~ m/###break###/) {                                       # see if there have is mark for a linebreak ("###break###") in the actual line
         @splitlines = split(/###break###/, $html->[$i]);                        # if so, break it up!
         for (my $j = 0; $j < ($param->{'breaknum'} + 1); $j++) {                # for each of the number of linebreaks (+1, mind the last fragment!)
            if (($i == ($numlines - 1)) && ($j < $param->{'breaknum'})) {
               $splithtml[$i + ($numlines * $j)] = $splitlines[$j].'<BR><BR>'.('=' x ($param->{'linebreak'} + $maxcodelen + 4)).'<BR>';
            }
            else { $splithtml[$i + ($numlines * $j)] = $splitlines[$j]; }
         }
      }
      else {
         for (my $j = 0; $j < ($param->{'breaknum'} + 1); $j++) {
            $splithtml[$i + ($numlines * $j)] = $html->[$i];
         }
      }
   }
   return \@splithtml;
}
