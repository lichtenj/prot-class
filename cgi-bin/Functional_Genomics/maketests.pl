#!/usr/bin/perl 
#
# File	  : maketests.pl -- Make test files for Ohio University Bioinformatics Framework
# Author  : Ohio University
# Created : July 2006
#
#     Copyright (c) 2005  Ohio University. All rights reserved.
#     This program is free software; you can redistribute it and/or
#     modify it under the same terms as Perl itself.
#

# note: run this program from the <install dir>\src or <install dir>\bin directory
#       because it uses some relative paths

use strict;
use pmf;
use ou;
use Data::Dumper;
use File::Copy;

# this routine makes the manifest test file from
# the passed object and filename
{
	my $dir = '..\\tests';
	mkdir $dir unless -d $dir; 
	sub maketestfile
	{
		my $this = shift or return;
		my $file = shift or return;
		$file = "$dir\\$file.man";
		$this->{filename} = $file;
		print "Creating TestFile: $file\n";
		$this->serialize('out',$file);
	}
    sub makebogusfile
    {
        my $filename = "$dir\\CompletelyBogus.man";
        open BOGUS, ">$filename" or die $!;
		print "Creating TestFile: $filename\n";
        print BOGUS "This is not a legitimate config file man\n";
        close BOGUS;
    }
}

# record the location that the main program was run from
# and tell anyone who needs to know
sub ExePath
{
	my ($volume,$exe_path,undef)
		= File::Spec->splitpath( File::Spec->rel2abs($0) );
	my $full_path = File::Spec->catpath( $volume,$exe_path,'');
	
	if ( $^O eq 'MSWin32' )
	{
		$full_path = Win32::GetShortPathName($full_path);
	}
	return $full_path;
}


# for some experiment/tests we will create a new peaklist input file
# so that we will get a unique output file for that test...
# the source file is known to be there, this will be used to create new
# peaklist files (with different names) for subsequent tests, so that the
# output filenames are not redundant (and the files overwritten)
my $srcfile = '../dat/ms/test_well_123.pkl';
die "no test peaklist source: $!\n" unless -f $srcfile;
my $well = 0;
sub makeinputfile
{
    my $this = shift or return;
    my $file = shift or return;
    my $src  = shift || $srcfile;
    $well ++;
    $file = "../dat/ms/$file well $well.pkl";
    if ( not -f $file )
    {
        copy($src, $file) or die "File $srcfile cannot be copied.";
    }
    
    $this->{config}->{General}->{MS_Input} = $file;
}
    

# encapsulate the logging logic here
sub logit
{
    print LOG @_;
    print     @_;
}

# Basically we just create a pmf object, then make whatever parameter changes
# we want to test and save the settings file.  Be sure to call GetDefaults()
# to reset paramaters when needed...
 
my ($this,$class) = new pmf();

open LOG, ">TestLog.txt" or die $!;
logit("PMF test log\n\n");


###########################################################################
# test failure to connect
# needs white box test -- run without internet connection

###########################################################################
logit("*** Try some different 'number of hits' options\n");
$this->GetDefaults('');
for (my $i=0; $i<100; $i += 10 )
{
    my $file = "Hits_$i";
    $this->{config}->{General}->{Number_hits} = $i;
    makeinputfile($this,$file);
    maketestfile($this, $file);
}

###########################################################################
logit( "*** Try all the kinds of supported taxonomy\n" );
my @taxes = (
		'SALMONELLA',
		'H_SAPIENS',
		'M_MUSCULUS',
		'C_ELEGANS',
		'ALL',
		'S_CEREVISIAE',
		'B_SUBTILIS',
		'D_MELANOGASTER',
		'R_NORVEGICUS',
		'UNCLASSIFIED',
		'A_THALIANA',
		'O_SATIVA',
		'MAMMALS',
		'RODENT',
		'E_COLI',
        'INVALID_TAXONOMY',
        );

$this->GetDefaults('');
foreach my $t (@taxes )
{
    $t =~ tr/ /_/;
    my $file = "Tax_$t";
    
    $this->{config}->{General}->{Taxonomy} = $t;
    makeinputfile($this,$file);
    maketestfile($this, $file);
}

###########################################################################
logit( "*** Try both overwrite and not-overwrite options\n");
$this->GetDefaults('');
{
    $this->{config}->{General}->{Overwrite} = 1;
    maketestfile($this, 'Overwrite_test');
}
{
    $this->{config}->{General}->{Overwrite} = 0;
    maketestfile($this, 'NotOverwrite_test');
}

###########################################################################
logit( "*** try some different peaklist input files (normal and bogus)\n");

{
    my $file = 'DefaultPeaksFile';
    if ( open TMP, ">$file" )
    {
        # cut-n-pasted from the msfit website
       print TMP
"676.3718
825.4571
1017.5370
1026.5964
1040.6121
1105.4550
1196.7132
1218.5796
1243.5887
1255.6538
1260.5928
1266.6281
1276.5836
1288.6094
1305.6310
1388.6975
1398.6805
1413.7126
1420.6741
1482.8101
1576.8138
1932.7067
2188.0850
";
       close TMP;
    }
    if ( -f $file )
    {
        makeinputfile($this,$file, $file);
        maketestfile($this, $file);
        unlink($file);
    }
}

# empty file
$this->GetDefaults('');
{
    my $file = 'EmptyFile';
    if ( open TMP, ">$file" )
    {
       print TMP "";
       close TMP;
    }
    if ( -f $file )
    {
        makeinputfile($this,$file, $file);
        maketestfile($this, $file);
        unlink($file);
    }
}

# bogus file
$this->GetDefaults('');
{
    my $file = 'BogusFile';
    if ( open TMP, ">$file" )
    {
       print TMP "
            \"Just the place for a Snark!\" the Bellman cried,
                As he landed his crew with care;
            Supporting each man on the top of the tide
                By a finger entwined in his hair.
            
            \"Just the place for a Snark! I have said it twice:
                That alone should encourage the crew.
            Just the place for a Snark! I have said it thrice:
                What I tell you three times is true.\"
        ";
       close TMP;
    }
    if ( -f $file )
    {
        makeinputfile($this,$file, $file);
        maketestfile($this, $file);
        unlink($file);
    }
}

###########################################################################
logit( "*** try both mass calculation options\n");

$this->GetDefaults('');
{
    $this->{config}->{General}->{Mass} = 'Average';
    maketestfile($this, 'MassAverage_test');
}
{
    $this->{config}->{General}->{Mass} = 'Monoisotopic';
    maketestfile($this, 'MassMono_test');
}
{
    $this->{config}->{General}->{Mass} = 'Gibberish';
    maketestfile($this, 'MassGibberish_test');
}

###########################################################################
# try some various input spec tests
$this->GetDefaults('');
{
    $this->{config}->{General}->{MS_Input} = '../dat/ms/test*.pkl';
    maketestfile($this, 'InputTest_test');
}
{
    $this->{config}->{General}->{MS_Input} = '../dat/ms/*.pkl';
    maketestfile($this, 'InputAll_test');
}
{
    $this->{config}->{General}->{MS_Input} = '..\dat\ms\test*.pkl';
    maketestfile($this, 'InputDOSType_test');
}
{
    $this->{config}->{General}->{MS_Input} = '..\dat/ms\test*.pkl';
    maketestfile($this, 'InputDosUnixMixType_test');
}
{
    $this->{config}->{General}->{MS_Input} = 'C:\Documents and Settings\All Users\Desktop\*.*';
    maketestfile($this, 'InputBogus1_test');
}
{
    $this->{config}->{General}->{MS_Input} = '../dat/ms/bogusdata_Well123.pkl';
    maketestfile($this, 'InputBogus2_test');
}

###########################################################################
# try some invalid directory/file entries
{
    my $done = 0;
    
    my $dir = '../pmfbogusin';
    my $bogusfile = "$dir/bogusdata_Well123.pkl";
    
    # first make sure the directory doesn't exist
    rmdir($dir) if (-d $dir);
    
    if ( ! -d $dir  ) 
    {
        $this->{config}->{General}->{MS_Input} = $bogusfile;
        maketestfile($this, 'InputBogusDir_test');
        $done = 1;
    }
    logit("could not build invalid directory test: $bogusfile") unless $done; 
}
# more invalid directory/file entries
{
    my $done = 0;
    
    my $dir = '../pmfbogusin';
    my $bogusfile = "$dir/dont_exist.pkl";
    
    # first make sure the directory doesn't exist
    rmdir($dir) if (-d $dir);
    
    if ( ! -d $dir  ) 
    {
        $this->{config}->{General}->{MS_Input} = $bogusfile;
        maketestfile($this, 'InputDontExist_test');
        $done = 1;
    }
    logit("could not build invalid directory test: $bogusfile") unless $done; 
}
# test bogus output directory
{
    my $done = 0;
    
    # add invalid characters to output directory spec
    my $dir = '../pmfnon<>existout';
    
    # first make sure the directory doesn't exist
    rmdir($dir) if (-d $dir);
    
    if ( ! -d $dir  ) 
    {
        $this->{config}->{General}->{Output_Directory} = $dir;
        maketestfile($this, 'OutputBogusDir_test');
        $done = 1;
    }
    logit("could not build invalid output directory test: $dir") unless $done; 
}

###########################################################################
logit( "*** Make a completely bogus configuration file\n");
{
    makebogusfile();
}

close LOG;

# end of main... below are local functions


__DATA__
