#!/usr/bin/perl 
#
# File    : manifest.pm -- Base class for Bioinformatics Experiments
# Author  : Ohio University
# Created : July 2005
# Version : <none>
#
#     Copyright (c) 2005  Ohio University. All rights reserved.
#     This program is free software; you can redistribute it and/or
#     modify it under the same terms as Perl itself.
#
#     This is the base class for all experiments in Bio Ohio Informatics
#

use strict;
use lib '/Perl/site/lib/';
use lib '../lib/';

package manifest;


#MARKED FOR DELETION
#-------------------------
use Tk;
use Tk::Pane;
use Tk::Adjuster;
use Tk::BrowseEntry;
use Tk::DialogBox;
use Tk::Scrollbar;
use Tk::Text;
use Tk::DirSelect;
#-------------------------


use IO::File;
use ou;
use Data::Dumper;



########################################################################
# Construct a new object of this class
sub new
{
  my $class  = shift or die;
  my (%config) = @_;
  my $parent = $config{-parent};
  my $self  =
  {
    class => $class,      # remember class name
    config => {},         # hold user settable parameters
    parent => $parent,    # parent handle
  };
  my $blessed = bless ($self, $class); 
  return ($self,$class);
}



# this function tells the rest of the application what directory the
# exe file was run from so that relative paths can be converted to
# absolute when necessary
sub ExePath
{
  use Cwd qw( abs_path cwd );
  my ($volume,$exe_path,undef) = File::Spec->splitpath( File::Spec->rel2abs($0) );
  my $full_path = File::Spec->catpath( $volume,$exe_path,'');
  return $full_path;
}

########################################################################
# Virtual setup function, basically you're given a parent window... do
# what you want with it.  The default is to call generic parameter data
# entry controls using config_hash for parameters
sub setup
{
  my $this = shift or die;
  my $mw = $this->{parent};
  # create "adjuster" control for managing the two panes    
  my $side = "left";
  # The Adjuster provides a splitter between the frames on the left and
  # the right so we can resize the frames vertically
  my $lhs = $mw->Frame; # Left Hand Side
  my $adj = $mw->Adjuster(-widget => $lhs, -side => $side);
  my $rhs = $mw->Frame; # Right Hand Side
  # add process button
  $this->{run_btn} = $lhs->Button(

    -width=>40,

    -text =>"Run $this->{class} Experiment",

    ) unless $this->{run_btn};



  # add cancel button

  $this->{can_btn} = $lhs->Button(

    -state=>'disabled',

    -width=>40,

    -text =>"Cancel",

    ) unless $this->{can_btn};



  my $frm = $lhs->Scrolled("Frame",

                          -scrollbars=>'w',

                          -takefocus=>0,

                          );



    # create the Right Hand Side control

  my $txt = $rhs->Scrolled('Text',

                   -height => '1',

                   -width => '1',

                   -scrollbars => 'osoe',

                   #-wrap => 'none',

                   );

  

  # give global access to the text window

  # for the link tagging functionality

  $this->{txt} = $txt;



  $lhs->pack(qw/-side left -fill y/);

  $adj->pack(qw/-side left -fill y/);

  $rhs->pack(qw/-side right -fill both -expand 1/);

  

  $this->{run_btn}->pack(-anchor=>'n',-fill=>'x', -pady=>3 );

  $this->{can_btn}->pack(-anchor=>'n',-fill=>'x' );



  $this->{run_btn}->configure(-command => sub{ do_submit($this,$txt)});

  $this->{can_btn}->configure(-command => sub{ do_cancel($this,$txt)});



  $frm ->pack(qw/-side left -fill both -expand 1/);

  $txt ->pack(qw/-side top -fill both -expand 1/);



  # insert the initial text into the text window

  $txt -> insert('end', $this->text() );

  $this->tagit($txt);

  #add_tags( $txt, $generx, 'protein_link' );



  # configure defaults for a manifest text control

  # font and highlight color to be more readable

  $this->configure_links($txt);



  # display the parameters in editable form  

  editHash($this, $frm, $this->{config}, uc($this->{class}.' Settings: '));



  # give the scrollbar keyboard focus now and whenever the

  # pane re-gets focus

  $frm->Subwidget("yscrollbar")->focus;

  $frm->Subwidget("yscrollbar")->configure(-takefocus => 0);

}



########################################################################

# Return true or false indicating whether this object has been modified.

# This function is called by the application class to determine

# whether the object needs saving.

sub modified

{

  my $this = shift or die;

  return ou::checksum($this->{config}) ne $this->{checksum};

}



########################################################################

# return the file magic number (from the bang line) of a manifest file

# used to verify that the file is a valid manifest file

sub getFileMagic

{

  my $file = shift or die;

  if ( ou::Open( *TMP, $file, "<" ) )

  {

    my $line = <TMP>;

    close TMP;

    $line =~ /#!manifest\.(\w+)/;

    return $1;

  }

  return;

}



########################################################################

# This function writes the current object to file

# or reads in a file

#

sub serialize

{

  my $this = shift or die;

  my $optn = shift or die;

  my $file = shift or die;

  #my $skip_checksum = shift;

  

  # print "File name is ". $file. "\n";  #gud

  # serialize out (write file)

  if ( $optn eq 'out' )

  {

    if ( ou::Open (*OUT, $file, ">") )

    {

      # write a bang line to the config file ala unix et al

      print OUT "#!manifest.$this->{class}\n";

      print OUT $this->serial_comment();



      # write the actual configuration to the file

      my $dumper = Data::Dumper->new([$this->{config}],['this->{config}']);

      $dumper->Indent(1);

      $dumper->Purity(1);

      $dumper->Deepcopy(1);

      my $sav = $dumper->Dump();

      print OUT $sav;

      

      # add a perl imbedded data tag in case anyone ever

      # trys to concatenate onto this file... it'll still work

      print OUT "__DATA__\n";

      close(OUT);

    }

  }

  

  # serialize in (read file)

  elsif ($optn eq 'in' )

  {

    if ( ou::Open (*IN , $file, "<") )

    {

      my $bangline = <IN>;

      # print $bangline."\n";  #gud

   

      # if this is not a valid file

      if ( not ( $bangline =~ /\#\!manifest/ ) )

      {

        print "Not a valid Manifest File: $file\n";

        return ;

      }

   

      my $read;

      {

        local $/ = undef;

        $read = <IN>;

      }    

      eval $read;

      close(IN);

    }

  }

  # just read or wrote so mark as not needing saved... until changed

  $this->{checksum} = ou::checksum($this->{config});

}



########################################################################

# pure virtual function returns configurations for the

# current object

sub config_hash

{

  # use this line to enforce pure virtual

  # die "Pure virtual violation\n";

  

  return

  {

   message => 'Default hash',

  };

}



########################################################################
# This function Gets Defaults first from the current
# object's default file
sub GetDefaults
{
  my $this = shift or die;
  my $filename = shift || ExePath() . $this->{class}.".cfg";
  # if this filename isn't set then set it, this can
  # happen in some default cases
  $this->{filename} = $filename unless $this->{filename};
  # if there is a defaults file then load it
  # after reality check for not empty hash
  if ( $filename && -r $filename && -s $filename > 40 )
  {
    $this->serialize('in',$filename,'skip_checksum');
  }
  # else create defaults and write the file
  else
  {
    $this->{config} = $this->config_hash();
    $this->serialize('out',$filename,'skip_checksum');
  }
}

########################################################################
# this function retrieves text from an existing log file
# for this experiement, or it returns some default text
# for the experiment pane
sub text
{
  my $this = shift or die;
  my $t = '';
  # if this experiement has a previous results then
  # retrieve them
  my $found = 0;
  my $log = '';
  if ( $this->{filename} )
  {
    $log = "$this->{filename}.log";
    if ( -f $log )
    {

      $found = 1;

      if ( ou::Open (*TMP, $log, "<" ) )

      { 

        local $/ = undef;

        $t = <TMP>;

        close TMP;

      }

    }

  }

  

  # if the previous block failed to retrieve any results

  # for this experiement then retrieve some default text

  if ( !$t )

  {

    my $class = uc($this->{class});

    $t =  

    '    BioHio $Revision: 1.2 $.'.

    "

    

    This is the results panel for the $class experiment

  

    In this panel, you will see the experiment report based on

    the settings you enter in the left side panel. Enter your

    experiment settings and press the run button.";

       

    # fix windows/unix incompatability

    $t =~ s/\r\n/\n/g;

    

    return $t;

  }

  return $t;

}



########################################################################

# virtual function which returns a list of values

# for the given variable

sub list_func

{

  my ($this,$var) = @_;

  return undef; 

}



########################################################################

# virtual function which returns the type of the entry control to be

# used for a particular variable, this function doesn't indicate that the

# variable is used or even defined, but may be in the future

sub type_func

{

  return '';  # override me

}



########################################################################

# add a comment line to a serialized manifest file

sub serial_comment

{

  my $file = shift || 'saved configuration';

  my $time = localtime(time());

  return  "#\n"

        . "#  $file: time:$time\n"

        . "#\n";  

}



########################################################################

# Print an message that is somehow more visible

sub ErrorMessage

{

  my $this = shift;

  if ( my $msg = shift )

  {

    print "**\n";

    print "** $msg\n";

    print "**\n";

  }

}



############################################################################

# set some default tag configurations for the passed text control

sub configure_links

{

  my $this = shift or die;

  my $txt = shift or die;

  # a tag for comments is defined here

  $txt->tagConfigure('comment',

                     -foreground=>'blue',

                     -underline=>0,

                     );

  # the selection color is set here

  $txt->tagConfigure('sel', -background => 'yellow', -foreground=>'black');

  # our text control font...

  $txt->configure( -font => '*-courier-bold-r-*-*-14-*"' );    

}



############################################################################

# add a data entry widget, or label to the entry pane based on the settings

# from type_func and list_func, the static variables in this scope are used

# to make widgets from subsequent function calls to appear on the next line

{

  my $row = 0;        # must be static... row counter

  my $max_level = 1; # prevent too deep recursion



  # this function is passed a hash containing experiment parameters

  # It traverses the hash (recursively) and draws the controls that

  # allow the data entry

  sub editHash

  {

    my $th = shift;         # this for short

    my $pr = shift;         # parent

    my $hr = shift;         # hash reference what's being edited

    my $lb = shift || '';   # label

    my $lv = shift || 0;    # level of recursion

  

    if ($lv==0){ $row=0;}   # when recursion starts row is 0 

    

    # this is a reality check that the programmer doesn't try to

    # pass a hash structure that is too deep, so this function

    # will just not go deeper than max_level

    return if ($lv > $max_level);

  

    # draw the data entry label first

    AddLabel($th,$pr,$lb,$hr) if $lb;



    # continue traversal by processing the current hash member

    # if it is an array reference or a hash reference then special

    # special code is called, a plain old variable just gets a

    # data entry control

    foreach my $k ( sort {$th->sort_func($a,$b)} keys %$hr)

    {

        if( $th->type_func($k) eq 'settings' ){

          AddEntry($th, $pr, $k, $hr->{$k} );

        }

        # if current element is a has then recurse

        elsif ( ref($hr->{$k}) eq 'HASH' )

        {

          editHash($th,$pr, $hr->{$k},"$k : ",$lv+1);

        }

      

        # if current element is an array then process each element of it

        elsif ( ref($hr->{$k}) eq 'ARRAY' )

        {

          # for each entry in the array either recurse (for references)

          # or add the data entry widget

          for ( 0 .. @{$hr->{$k}} )

          {

            # recurse for references

            if ( ref $hr->{$k}[$_] eq 'HASH' ||

                 ref $hr->{$k}[$_] eq 'ARRAY' )

            {

              editHash($th,$pr, \$hr->{$k}[$_],"$k : ",$lv+1);

            }

            # or add the entry control

            else

            {

              if ( $hr->{$k}[$_] )

              {

                AddEntry($th,$pr,"$lb$k\[$_\]",\$hr->{$k}[$_] );

              }

            }

          }

        }



      # else current element is not a ref so just add an entry widget

      else {

        AddEntry($th,$pr,$k,\$hr->{$k} );

      }

    }

  }

  

  # add a data entry widget to the data entry pane (the left side)

  sub AddEntry

  {

    my $this = shift or die;

    my $per = shift or die; # parent

    my $lab = shift or die; # label

    my $var = shift or die; # variable reference

  

    my $type = shift;       # kind of variable being edited

  

    my $ew = 30;            # entry width (default)

    

    # resolve the specific (overridden) features of the variable

    my $control_type = $this->type_func($lab);

    my $control_values = $this->list_func($lab);

  

    # ignore parameters marked as hidden

    return if $control_type eq 'hidden';

  

    # transliterate any underscores to spaces for display

    $lab =~ tr/_/ /;

    

    # first add the text label for the variable before proceeding

    # to add the actual data entry control

    my $lbl = $per->Label(

      -text=>"$lab : ", 

      -anchor => 'e',

      -justify => 'left',

    )-> grid(-row=>$row,-column=>1,-sticky=>'ew');

  

    # The following multiple elsif statement will draw the actual

    # widget into the pane, a different widget will be drawn

    # depending on what the conrol type is

    

    # draw a checkbox entry control

    if ( $control_type eq 'checkbox' )

    {

      my $ent = $per->Checkbutton(

        -variable=>$var,

        -takefocus=>1,

      )-> grid(-row=>$row,-column=>2,-columnspan=>2,-sticky=>'w');

    }

    

    # draw a button control to browse to a file 

    elsif ( $control_type eq 'browsefile' )

    {

      my $ent = $per->Entry(

                            -width => $ew-length('...'),

                            -textvariable=>$var,

                            );

      my $but = $per->Button(-text => '...',

                             -padx=>0,

               -command => sub { fileDialog($per, $ent)});

      

      $ent->grid( -row=>$row, -column=>2, -sticky=>'ew' );

      $but->grid( -row=>$row, -column=>3, -sticky=>'e' );

    }

    

    # draw a button control to browse to a directory

    elsif ( $control_type eq 'browsedirectory' )

    {

      my $ent = $per->Entry(

                            -width => $ew-length('...'),

                            -textvariable=>$var,

                            );

      my $but = $per->Button(-text => '...',

                             -padx=>0,

               -command => sub { dirDialog($per, $ent)});

      

      $ent->grid( -row=>$row, -column=>2, -sticky=>'ew' );

      $but->grid( -row=>$row, -column=>3, -sticky=>'e' );

    }

    

    # draw a dropdown control and get the items from $control_values

    elsif ( $control_type eq 'dropdown' )

    {    

      $per->Optionmenu

      (

        -options => $control_values,

        -variable => $var,

      )

      ->grid( -row=>$row, -column=>2, -columnspan=>2, -sticky=>'ew' );

    }

    

    # gud: draw the button for display window settings

    elsif ( $control_type eq 'settings' )

    {

      my $but = $per->Button(-text => 'Edit...', -padx=>0,

                             -command => sub {

                              my $dialog = application::MW()->DialogBox(-buttons=>["Ok"]);

                              $dialog->configure( -title => "Settings" );

                              my $row = 1;

                              foreach my $key (sort keys %$var){

                              # while( my($key, $val) = each(%$var) ){

                                

                                $dialog->Label( -text=>"$key : ",

                                                     -anchor => 'e', -justify => 'left',

                                                     )-> grid(-row=>$row, -column=>1, -sticky=>'ew');

                                $dialog->Entry( -textvariable=>\$var->{$key},

                                              -width=>$ew, -takefocus=>1

                                              )-> grid(-row=>$row,-column=>2, -columnspan=>2, -sticky=>'ew');

                                $row++;

                              }

                              $dialog->Show();

                             }

                            );

      $but->grid( -row=>$row, -column=>2, -columnspan=>2, -sticky=>'ew' );

    }    

    # finally the default data entry type is a plain old entry widget

    else 

    {

      $per->Entry(

                 -textvariable=>$var,

                 -width=>$ew,

                 -takefocus=>1,

                 )

        -> grid(-row=>$row,-column=>2,-columnspan=>2,-sticky=>'ew');

    }

    

    $row++;  # remember the row (statically) so next is one below this

  }

  

  # add a label widget to the data entry pane (the left side)

  sub AddLabel

  {

      my $thi = shift or die;

      my $per = shift or die; # parent

      my $lab = shift; # label

      my $var = shift or die; # variable reference

      

      die unless ref $var;

      

      $per->Label(

                  -foreground=>'red3',

                  -text=>$lab, 

                  -anchor => 'e',

                  -justify => 'left',

                  )

      -> grid(-row=>$row,-column=>1,-sticky=>'ew');

      $row++;

  }

} # end of scope control for edit control routines



############################################################################

# called from the button in the data entry pane, this function just

# wraps the getOpenFile() function and updates the data entry control

# with the results

sub fileDialog

{

    my $w = shift;          # window

    my $ent = shift;        # entry widget



    my $file = $w->getOpenFile();

    if (defined $file and $file ne '')

    {

      $ent->delete(0, 'end');

      $ent->insert(0, $file);

      $ent->xview('end');

    }

}



############################################################################

# called from the button in the data entry pane, this function just

# wraps the DirSelect() function and updates the data entry control

# with the results

{

  my $dlg; # static for speed

  sub dirDialog

  {

      my $w = shift;    # window

      my $ent = shift;  # entry widget

      my $file;

      

      $dlg = $w->DirSelect(-title=>"Select a directory") unless $dlg;

      my $dir = $dlg->Show();

      

      if ( defined $dir and $dir )

      {

        $ent->delete(0, 'end');

        $ent->insert(0, $dir);

        $ent->xview('end');

      }

  }

}



############################################################################

# the default sorting routine for experiment paramters as edited

sub sort_func

{

  my ($this,$a,$b) = @_;



  # sort references to hashes and arrays

  # to the bottom just for clarity

  return ref $a cmp ref $b; 

}



############################################################################

# this is the function behind the main Run button for any experiment

# it launches the xxx_main.pl program that coresponds to the derived manifest

# object

sub do_submit

{

  my $this = shift or die;

  my $txt = shift or die;   # the text control what's updated



  # we'll need the main window and the experiment filename 

  my $mw = application::MW();

  my $filename = $this->{filename} || '';

  

  # the experiment object MUST be saved before we can proceed because the

  # experiment is a separate process and must read the experiment settings

  # out of this file

  if ( $this->modified() )

  {

    my $response = $mw->messageBox(-message=>"You must save before running the experiment\nSave Now?",-type=>'YesNo', -default=>'Yes', -icon=>'question');

  

    if ( $response eq 'Yes' )

    {

      application::FileSave();

    }

    else

    {

      return;

    }

  }

  

  # construct the log filename based on the experiment filename  

  my $log = $this->{filename} ? "$this->{filename}.log" : "$this->{class}.log";

  

  # construct the name of the main experiment program file based on this

  # object

  my $main;

  

  # If we are running from the command line then just use that

  if ( $this->{config}->{command_line} )

  {

    $main = $this->{config}->{command_line};

  }

  

  # if the application is to launch the experiment then construct the

  # name first looking for the compiled program and alternatively the

  # perl source program (for machines with perl)

  else

  {

    my $cfg = $this->{filename} ? $this->{filename} : "$this->{class}.cfg";

    $cfg = ou::FixFileName($cfg);

    

    my $exe = ExePath() . "$this->{class}_main.exe";

    my $pl  = ExePath() . "$this->{class}_main.pl";

    

    $main = (-f $exe) ? "$exe $cfg" : "perl $pl $cfg";

  }



  # we're starting so disable the run button and enable "cancel"

  $this->{run_btn}->configure(-state=>'disabled');

  $this->{can_btn}->configure(-state=>'active');



  # clear text display and start displaying results of new run

  $txt->delete('1.0','end');

  application::status("Running experiment...$main");

  

  # save stdout for restoring later in do_cancel or readfile

  open TMPOUT, '>&', \*STDOUT   or die "Can't dup stdout: $!\n";

  

  # open the log file what holds the results

  if ( open(STDOUT,">$log") )

  {

    # print some header information... also the checksum is

    # used to tie this experiment report to a config file

    {

      my $time = localtime(time());

      my $checksum = ou::checksum($this->{config});

      my $experiment = $this->{filename} || $this->{class};

      print "Experiment: $experiment\n" .

            "      Time: $time\n\n";

    }



    # start the new experiment process

    $this->{process} = Proc::Background->new($main);



    # Just checking... couldn't the process be done already?

    if ( not $this->{process} )

    {

      $txt->insert('end', "No Process was started : $!\n");

      $this->do_cancel();

      return;

    }

    

    # now tail the logfile by reading the log file in

    # blocks (for speed) and updating the text control

    select(STDOUT);  # force a flush on stdout

    my $fh = IO::File->new("<$log");

    if ( $fh )

    {

      $this->{tail_pos} = 0;

      $fh->autoflush(1);

      $txt->{posn} = 0;

      readfile($this,$txt,$fh,$this->{process});

      # read the file repeatedly using a timer in miliseconds

      $this->{tail_id} = $mw->repeat(1000,sub { readfile($this,$txt,$fh,$this->{process}) } );

    }

  }

}



############################################################################

# this is the function behind the main Cancel button for any experiment

sub do_cancel

{

  my $this = shift or die;

  my $txt = shift;  # text control

  my $msg = shift;  # cancelation message



  # restor saved standard output

  open STDOUT, '>&', \*TMPOUT   or die "Can't restore stdout: $!\n";

  

  # remember the state (of being canceled) in this object by

  # resetting the tail_id

  if ( $this->{tail_id} )

  {

    $this->{tail_id}->cancel;

    $this->{tail_id} = undef;

  }



  # reset the buttons back to state of "ready to run again"

  $this->{run_btn}->configure(-state=>'active');

  $this->{can_btn}->configure(-state=>'disabled');



  # finally, kill the process

  if ( $this->{process} )

  {

    $this->{process}->die();

    $this->{process} = undef;

  }



  # tell 'em we are done

  $msg = $msg ? $msg : 'Experiment canceled'; 

  if ( $txt )

  {

    $txt->insert('end', "\n\n*** $msg ***\n" );

    $txt->see('end');

  }

  application::status($msg);

}



############################################################################

# this function is called repeatedly by the timer call in do_submit()

# to read a block of data from the file being "tailed",

sub readfile

{

  # this, text control, file handle, process

  my ($this,$t,$fh,$proc) = @_;

  

  # use the file size and a pointer to the current file position

  # to calculate the length of data that is currently available for reading

  my $size = (stat($fh))[7];

  my $len  = $size - $this->{tail_pos};



  # if there is some to read...

  if ($len > 0)

  {

    # this code reads whole blocks

    my $buffer = "";

    my $got = sysread($fh,$buffer,$len);

    if ( $got )

    {

      $buffer =~ s/\r\n/\n/g;  # fix unix/window incompatiblities



      # display the read data and tag it using the virtual

      # function from this object 

      my $length = length($buffer);

      $t->insert('end',$buffer );

      $this->tagit($t,"end - $length chars" );



      # update position cursor and make sure the newly added

      # data is visible

      $this->{tail_pos} += $got;

      $t->see('end');

    }

  }

  else

  {

    # if there is nothing to read and the

    # process we're tailing has stopped

    if ( not $proc->alive() )

    {

      do_cancel($this,$t,"Experiment finished normally.");

    }

  }

}



##########################################################################

# scope for functions controling the list of open manifest objects, the

# functions within this scope support external document management by the

# user interface

{

  # static list of active manifest objects

  our $files = {};  



  # is the requested file a currently open file

  sub Available

  {

    my $file = shift or die;

    return ($file&&$files->{$file}) ? 0 : 1;

  }



  # open a file and add it to the list  

  sub Open

  {

    my $this = shift or die;

    my $file = shift;

    

    $this->GetDefaults($file);

    $this->{filename} = $file; 

    $this->{checksum} = ou::checksum( $this->{config} );

  

    $files->{$file} = 1 if $file;

    $this->setup();

  }

  

  # close a file and remove it from the list

  sub Close

  {

    my $this = shift or die;

    delete $files->{$this->{filename}} if ($this->{filename} && $files->{$this->{filename}});

  }



  # write a list of open files so that the next time the program

  # is run, they can be re-opened

  sub WriteManifestList

  {

    # dump a list of filenames for open manifests

    # so that they can be opened up again when

    # the program is restarted

    

    my @list;

    foreach my $file (keys %$files)

    {

      if ( $file && -f $file )

      {

        push @list, $file;

      }

    }

    if ( @list )

    {

      my $listfile = ExePath()."manifests.lst";

      if ( ou::Open(*MANIFESTS, $listfile, ">" ) )

      {

        foreach my $l (@list)

        {

          if ( $l )

          {

            print MANIFESTS "$l\n"

          }

        }

        close MANIFESTS;

      }

    }

  }

}



##########################################################################

# The following 3 (add_tags,tagit,rem_tags) functions are not neccesarily

# used virtually but may be used by several child classes

sub add_tags

{

  my $this = shift or die;

  my $text = shift or die;

  my $regx = shift or die;

  my $link = shift or die;

  my $start_match = shift || '0.0';

  my $end_match = shift || 'end';

  {

    my $length;

    $regx =~ s/\*/\\\*/g;

    while ( my $start = $text->search(-count => \$length, -regex, $regx, $start_match, $end_match ) )

    {

      # if exclude delimiters

      if (0)

      {

        $text->tagAdd( $link, "$start + 1 chars", "$start + $length chars - 1 chars" );

        $start_match = "$start + 1 chars";

      }

      else

      {

        $text->tagAdd( $link, $start, "$start + $length chars" );

        $start_match = "$start + 1 chars";

      }

      #$text->delete( "$start + $length chars - 1 chars" );

      #$text->delete( "$start");

    }

  }

}



##########################################################################

# this virtual function can be overridden to tag text in a text control

sub tagit

{

  # base class does nothing for now but the parameters

  # are shown below

  return;

  

  my $this = shift or die;

  my $text = shift or die;  # the text control to update

}



##########################################################################

# this  function removes ALL tags in a text control

sub rem_tags

{

  my $this = shift or die;

  my $text = shift or die;

  my $regx = shift or die;  

  my $link = shift or die;

  

  $text->tagRemove( $link, '0.0', 'end' );

}



return 1;

