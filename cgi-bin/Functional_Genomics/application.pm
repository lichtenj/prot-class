#!/usr/bin/perl -w

#

# File	  : application.pm -- application class for OU Bioinformatics Bench

# Author  : Ohio University CNS

# Created : July 2005

# Version : <none>

#

#     Copyright (c) 2005  Ohio University. All rights reserved.

#     This program is free software; you can redistribute it and/or

#

package application;



use strict;



use lib "/Perl/site/lib/";

use lib '../lib/';



use Tk;

use manifest;



use Proc::Background;

use File::Basename;

use Tk::DynaTabFrame;

use Tk::Scrollbar;

use Tk::Menu;



use ou;



our $APP_TITLE = "\(B\|Oh\)ioInformatics";

our $APP_VERSION = '1.0';



my $mw;                 # application main window

my $book;               # the notebook tabbed control

my $tabs = {};          # the tabs in the notebook

my $filesinuse = {};    # the current open files



# the main function just runs the main loop which handles all user

# and system requests, The rest of the application class is just

# responding to those requests

sub main

{

  my $this = shift or die;

  

  # co-opt the close window system call to do our close function

  $mw->protocol('WM_DELETE_WINDOW',sub{$this->FileQuit();});



  open MAINOUT, '>&', \*STDOUT   or die "Can't dup stdout: $!\n";

  open MAINERR, '>&', \*STDERR   or die "Can't dup stdout: $!\n";



  $this->MainLoop;



  close MAINOUT;

  close MAINERR;



  cleanup(); # will exit

}



# override this function to do any cleanup after the mainloop has run

sub cleanup

{

  # include any application specific cleanup code here

  exit;

}



# this function tells the rest of the application what directory the

# exe file was run from so that relative paths can be converted to

# absolute when necessary

sub ExePath

{

  return manifest::ExePath();

}



# make types static so it can be accessible by static functions

my @types = (["Experiment Manifests",'.man','TEXT'],["All Files","*"]);



sub new

{

  my $class = shift;

  $mw = MainWindow->new;



  my $self = {};

  

  bless($self,$class);



  $self->Modules();

  $self->setup();

  

  $mw -> bind('<Control-o>' => sub{$self->FileOpen});

  $mw -> bind('<Control-s>' => sub{$self->FileSave});

  $mw -> bind('<Control-a>' => sub{$self->FileSaveAs});

  $mw -> bind('<Control-q>' => sub{$self->FileQuit});



  status('BioOhio is ready');

  $self->{ShowInfoTip} = ou::SetShowInfoTip(0);

  return ($self, $class);

}



sub DESTROY

{

  my $this = shift;

  ou::SetShowInfoTip( $this->{ShowInfoTip});

}



# override this virtual function to include all the modules you want to

# "compile" into the exe program, and return the list  of modules

sub Modules

{

  return [];

}



# check whether the passed module is one that has been loaded by

# this instantiation of the main gui (in Modules())

sub ValidModule

{

  my $this = shift or die;

  my $mod = shift;

  foreach ( @{$this->Modules()} )

  {

    return 1 if $mod eq $_;

  }

  return 0;

}



# override this virtual function to create the new experiment menu specific to

# the modules included in sub Modules()

sub ModulesMenu

{

  my $this = shift or die;

  return

  [

    [

      qw/command/, '<No Experiments Defined>', -command => sub

      {

        $this->MW()

        ->messageBox(-message=>"No Experiment modules included in this build!");

      }

    ],

  ];

}



sub WriteManifestList

{

  # dump a list of filenames for open manifests

  # so that they can be opened up again when

  # the program is restarted

  my $file = ExePath()."manifests.lst";

  if ( ou::Open(*MANIFESTS, $file, ">" ) )

  {

    foreach my $file (keys %$filesinuse)

    {

      if ( $file )

      {

        print MANIFESTS "$file\n"

      }

    }

    close MANIFESTS;

  }

}



# These functions just give access to handles needed outside

# the class, 

sub MW { $mw; }

sub BOOK { $book; }



########################################################

# respond to the "File - New" menu by creating a new Document

# the manifest object knows how to create itself, the application

# just has to set the tab title, and store the object reference

sub FileNew

{

  my $this = shift or die;

  my $type = shift or die;

  my $title = "$type:<untitled>";

  

  my $tab = $book->add(-caption => $title, -label=>$title);

  my ($man, $class) = $type->new(-parent=>$tab);

  $man->Open();

  

  $tab->{man} = $man;

  $tabs->{$title} = $tab;  

}



########################################################

# respond to the "File - Close" menu by

# closing the document which has focus

sub FileClose

{

  my $this = shift or die;

  my ($caption,$man) = @_;

  ($caption,$man) = CurrentMan() unless ($caption && $man);

  

  if ( $man )

  {

    my $file = $man->{filename};

    #$file = "$man->{class}:<untitled>" unless $file;

   

    if ( $file )

    {

      # in case an experiment is running

      $man->do_cancel();

      

      if ( $man->modified() )

      {

        my $response = $mw -> messageBox(-message=>"Save $file?",

                              -type=>'YesNoCancel',-icon=>'question');

        if( $response eq 'Yes' )

        {

          $this->FileSave();

        }    

        elsif( $response eq 'Cancel' )

        {

          return;

        }

      }

    }

    delete $filesinuse->{$file} if ($file && $filesinuse->{$file});

    $man->Close();

  }

  

  delete $tabs->{$caption};

  $book->delete($caption);

}



########################################################

# respond to the "File - Open" menu

sub FileOpen

{

  my $this = shift or die;

  my $file = shift;

  my %parms = @_;



  # if the file doesn't exist then browse to

  # let the user choose it, unless quiet mode

  if ( not $file and not $parms{quiet} )

  {

    $file = $mw->getOpenFile(-filetypes => \@types,

        -initialfile=> 'Default',

        -defaultextension => '.man');

  }

  return unless $file;



  # prepare variables based on the filename and

  # contents of the bang line (magic number)

  my ($base,$path,$suffix) = fileparse($file);

  my $type = manifest::getFileMagic($file);



  # a slightly more verbose output if debugging here

  # print "$type is not a valid module\n" if (!$this->ValidModule($type));

  

  # if the type is valid, that is, if it is an experiment

  # type that has been loaded by this instantiaion of application

  # then go ahead and open the experiment

  if ( $type && $this->ValidModule($type) )

  {

    my $title = "$base";

    my $tab = $book->add(-caption => $title,);

    

    my ($man,$class);

    eval { ($man,$class) = $type->new(-parent=>$tab); };

    

    # open if created

    if ( $man )

    {

      $man->Open($file);

      $filesinuse->{$file} = 1;

      

      $tab->{man} = $man;

      $tabs->{$title} = $tab;

      return $file;

    }

  }

  # couldn't open, but only notify if being verbose

  else

  {

    if ( not $parms{quiet} )

    {

      application::MW()->messageBox(

                 -message => "Unknown file type: $file",

                 -title => 'Alert',

                 -type => 'Ok',

                 );

    }

  }

  return undef;

}



########################################################

# respond to the "File - Save" menu by

# calling the serialize method of the

# focused document

sub FileSave

{

  my $this = shift;

  my ($caption,$man) = @_;

  ($caption,$man) = CurrentMan() unless ($caption && $man);



  if ( not $man->{filename} )

  {

    $man->{filename} = FileSaveAs();

  }

  if ( $man->{filename} )

  {

    $man->serialize('out',$man->{filename});

  }

  return $man->{filename};

}



# respond to the "File - Save as" menu by

# calling the serialize method of the

# focused document after chosing a filename

sub FileSaveAs

{

  my $this = shift;

  my $tab = $book->raised();

  my $lab = $book->raised_name();

  my $man = $tab->{man};



  while(1)

  {

    my $file = $mw->getSaveFile(-filetypes => \@types,

        -initialfile => '',

        -defaultextension => '.man');

    

    if ( !$file )

    {

      last;

    }

    else

    {

      if ( manifest::Available($file) )

      {

        my ($base,$path,$suffix) = fileparse($file);

        

        $man->{filename} = $file;

        $manifest::files->{$file} = 1;

        $book->pageconfigure($lab, -label => $base);

        $man->serialize('out',$file);

        last;

      }

      else

      {

        my $response = $mw->messageBox(-message=>"The file:\n$file\nIs already open in a tab",-type=>'Ok');

      }

    }

  }

  return $man->{filename};

}



# respond to the "File - Quit menu

sub FileQuit

{

  my $this = shift or die;



  manifest::WriteManifestList();

  

  $mw->destroy;

}



# display the about box

sub HelpAbout

{

  my $this = shift or die;

	my $db = $mw->DialogBox(-title=>"ABOUT",-buttons=>["OK"],);



  $db->bind('<Escape>',       sub { $db->{'selected_button'} = 'OK' });

  

  my $text =

  "

    $APP_TITLE

    $APP_VERSION

    

    Center for Intelligent Distributed and Dependable Systems

    

    cidds\@ohio.edu

    ";

    # fix the ole windows/unix incompatability

    $text =~ s/\r\n/\n/g;



	$db->add("Label",

              -text=> $text,

              )->pack;

	$db->Show;

}



sub HelpAminoAcid

{

  my $this = shift or die;

	my $db = $mw->DialogBox(-title=>"Amino Acid Table",-buttons=>["OK"],);



  $db->bind('<Escape>',       sub { $db->{'selected_button'} = 'OK' });

  

	$db->add("Label",

    -font => '*-courier-bold-r-*-*-14-*"',

    -justify=>'left',

    -text=>

    "

    Abbreviation    Amino acid name

    ------------    ---------------

    Ala     A       Alanine

    Arg     R       Arginine

    Asn     N       Asparagine

    Asp     D       Aspartic acid (Aspartate)

    Cys     C       Cysteine

    Gln     Q       Glutamine

    Glu     E       Glutamic acid (Glutamate)

    Gly     G       Glycine

    His     H       Histidine

    Ile     I       Isoleucine

    Leu     L       Leucine

    Lys     K       Lysine

    Met     M       Methionine

    Phe     F       Phenylalanine

    Pro     P       Proline

    Ser     S       Serine

    Thr     T       Threonine

    Trp     W       Tryptophan

    Tyr     Y       Tyrosine

    Val     V       Valine

    Asx     B       Aspartic acid or Asparagine   

    Glx     Z       Glutamine or Glutamic acid.

    Xaa     X       Any amino acid.

    ",

       )->pack;

	$db->Show;

}



sub HelpProteomics

{

  my $this = shift or die;

  open (PROTEOMICS, "proteomics_help.txt");



  my $text;



  while (my $record = <PROTEOMICS>)

  {

    $text = $text.$record;

  }



  close(PROTEOMICS);

  

  # Use scrolled Tk dialog scrollbox to print $text out  



  my $mw = new MainWindow;

  $mw->title("Proteomics Help");

  my $frm_name = $mw -> Frame();

  my $textarea = $mw -> Frame();

  

  my $txt = $textarea -> Text(

      -width=>100,

      -height=>40);

      

  my $srl_y = $textarea -> Scrollbar(-orient=>'v',-command=>[yview => $txt]);

  my $srl_x = $textarea -> Scrollbar(-orient=>'h',-command=>[xview => $txt]);

  

  $txt -> configure(-yscrollcommand=>['set', $srl_y], -xscrollcommand=>['set',$srl_x]);

  $frm_name -> grid(-row=>1,-column=>1,-columnspan=>2);

  $txt -> grid(-row=>1,-column=>1);

  $srl_y -> grid(-row=>1,-column=>2,-sticky=>"ns");

  $srl_x -> grid(-row=>2,-column=>1,-sticky=>"ew");

  $textarea -> grid(-row=>5,-column=>1,-columnspan=>2);

  $txt -> insert('end',$text);

  #MainLoop;

}





# return the document that has focus

sub CurrentMan

{

  my $tab = undef;

 	my $caption = $book->raised_name();

  

  if ( $caption )

  {

    $tab = $tabs->{$caption}->{man};

  }

  return ($caption,$tab);

}



# scope control for status line variable and control

{

  my $status_line = 'Ready';

  my $status_ctrl;

  

  # print passed text on the status line

  sub status

  {

    my $txt = shift;

    if ( $txt)

    {

      $status_line = $txt;

      if ( $status_ctrl )

      {

        $status_ctrl->update();

      }

    }

  }



  sub setup

  {

    my $this = shift or die;

  

    # start out with the window invisible so the drawing

    # doesn't show on the screen... deiconify later

    $mw->withdraw();

    

    ## kindof center the window on the screen

    # this is the default setting for the non-maximized window

    # but the window is maximized later in this function by default

    # but only in windows

    {

      my ($w, $h) = ($mw->screenwidth, $mw->screenheight);

      my ($hm, $vm) = ( $w/10, $h/10 );

      my ($hs, $vs) = ((8*$w/10),(7*$h/10) );

    

      ($hm, $vm) = ( sprintf("%d", $hm), sprintf("%d", $vm) );

      ($hs, $vs) = ( sprintf("%d", $hs), sprintf("%d", $vs) );

    

      $mw->geometry($hs . 'x' . $vs . '+' . $hm .'+' . $vm );

    }

    

    # create a frame for the tabbed display on the main window

    $book = $mw->DynaTabFrame(

     #-font => 'Arial 8', 

     #-raisecmd => sub{print "raised\n";}, # \&raise_cb,

     -tabclose => sub {

        my ($obj, $tab) = @_;

        $obj->raise($tab);

        $this->FileClose();

       },

     #-tabcolor => $color,

     #-raisecolor => $raisecolor,

     #-tabside => $side,

     #-tabpadx => 3,

     #-tabpady => 3,

     #-tiptime => 600,

     #-tipcolor => 'white'

     )

     ->pack (-side => 'top', -expand => 1, -fill => 'both');

  

    # create main menu

    {

      $mw->configure(-menu => my $menubar = $mw->Menu(-tearoff=>0) );

      

      my $file = $menubar->cascade(

          -label => 'File', -tearoff=>0, -underline=>0, -menuitems => $this->file_menuitems());

      my $help = $menubar->cascade(

          -label => 'Help', -tearoff=>0, -underline=>0, -menuitems => $this->help_menuitems());

    }

  

    $status_ctrl = $mw->Label(-textvariable => \$status_line,

                -justify => 'right',

                -relief => 'sunken',

                -anchor => 'w',

                )

        ->pack(-side => 'left',

               -expand => 'x',

               -fill => 'x');

  

    # at startup, look for a file containing the list

    # of previously opened experiments and re-open them

    {

      my $count = 0;

      my $file = ExePath()."manifests.lst";

  

      # open file and load default manifests 

      if (-f $file )

      {

        if ( ou::Open(*MANIFESTS, $file, "<" ) )

        {

          while(<MANIFESTS>)

          {

            chomp;

            if ( -f $_ )

            {

              $count ++ if $this->FileOpen( $_, 'quiet' => 1 );

            }

          }

          close MANIFESTS;

        }

      }

      # if no manifests loaded then load defaults

      if ( ! $count )

      {

        foreach my $module ( @{$this->Modules()} )

        {

          $this->FileNew( $module );

        }

      }

    } 

    # maximize the window by default... in windows

    if ( $^O eq 'MSWin32' )

    {

      use Win32::API;

      use constant SW_MAXIMIZE => 3;

      my $api = new Win32::API('user32','ShowWindow',['N','N'],'N');

      my $WinID=hex($mw->frame);

      $api->Call($WinID,SW_MAXIMIZE);

    }

    # finally set application title and

    # make the window visible

    $mw->title($APP_TITLE);

    

    $mw->deiconify();

  }

}



sub help_menuitems

{

  my $this = shift;

  return

  [

    #['command', 'Contents', -underline=>0, -command => sub{$this->HelpContents}],

    ['command', 'Amino Acid Table', -underline=>11, -command => sub{$this->HelpAminoAcid}],

    ['command', 'Proteomics User Guide', -underline=>11, -command => sub{$this->HelpProteomics}],

    '',

    #['command', 'Version', -underline=>0, -command => sub{$this->help}],

    ['command', 'About',   -underline=>0, -command => sub{$this->HelpAbout} ],

  ];

}



# create the menuitems, but let the polymorphic ModulesMenu function decide

# which experiemnts to add to the "New" menu, This ensures that only loaded

# modules will have menus... and vice versa

sub file_menuitems

{

 my $this = shift;

 return

 [

     [qw/cascade New -tearoff 0 -underline 0  -menuitems/ =>

      $this->ModulesMenu()

     ],                                                      '',

     [qw/command ~Open  -accelerator Ctrl-o -command/ => sub{$this->FileOpen}],               '',

     [qw/command ~Close -command/ => sub{$this->FileClose}],

     [qw/command ~Save  -accelerator Ctrl-s -command/ => sub{$this->FileSave}],

     [qw/command/, 'S~ave As ...', qw/-accelerator Ctrl-a -command/ => sub{$this->FileSaveAs}], '',

     [qw/command ~Quit  -accelerator Ctrl-q -command/ => sub{$this->FileQuit}],

 ];



} # end file_menuitems



return 1;