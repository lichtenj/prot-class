#!/usr/bin/perl 


#


# File	  : ou.pm -- some common subs for Bioinformatics Experiments


# Author  : Ohio University


# Created : July 2005


# Version : <none>


#


#     Copyright (c) 2005  Ohio University. All rights reserved.


#     This program is free software; you can redistribute it and/or


#     modify it under the same terms as Perl itself.


#


package ou;


use strict;


use Tk::DialogBox;


use Digest::MD5 qw(md5_hex);


use File::Basename;


use File::Spec;





# make sure a lib is loaded if running windows


BEGIN


{


	if ( $^O eq 'MSWin32' )


	{


		require Win32::TieRegistry; 


	}


};





require Exporter;


our @ISA = qw(Exporter);


our @EXPORT_OK = qw(MakeOutputFile);





sub checksum


{


	my $val = shift or die;


	return md5_hex(Data::Dumper->Dump([$val]));


}





# take a long filename that may not work in all situations


# and turn it into a windows compatible short filename which


# should always work


sub FixFileName


{


	my $full_path = shift;


	if ( $full_path )


	{


		if ( $^O eq 'MSWin32' )


		{


			my $ret = Win32::GetShortPathName($full_path);


			return $ret if $ret;


		}


	}


	return $full_path;


}





# undo the filename fixing of FixFileName and return the full


# file name and path which should be more readable by a human


sub UnFixFileName


{


	my $full_path = shift;


	if ( $full_path )


	{


		if ( $^O eq 'MSWin32' )


		{


			my $ret = Win32::GetLongPathName($full_path);


			return $ret if $ret;


		}


	}


	return $full_path;


}








# one routine to open files after fixing the filename


sub Open


{


	my $handle = shift;


	my $file = shift;


	my $pipe = shift;


	my $ret = undef;


	if ( $handle && $file )


	{


		$file = FixFileName($file);


		$ret = open($handle,"$pipe $file");


		if ( !$ret )


		{


			print "Could not open $file: $!\n"


		}


	}


	return $ret;


}


	


# The info tips that pop up in windows explorer dialogs


# cause the application to shutdown. Therefore we will


# disable this feature while our application is running


sub SetShowInfoTip


{


	if ( $^O eq 'MSWin32' )


	{


		my $value = shift;


		use vars qw( $RegObj %_Roots %RegHash $Registry );


		my $userKey= $Registry->{"CUser/Software/Microsoft/Windows/CurrentVersion/"};


		


		my $key = $userKey->{"Explorer/Advanced/ShowInfoTip"};


		if ( $key )


		{


			$userKey->{"Explorer/Advanced/ShowInfoTip"}	=


				[ pack("L",$value), "REG_DWORD" ];	


			return hex($key);


		}


	}


}





return 1;