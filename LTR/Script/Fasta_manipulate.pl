#happye321@gmail.com
#This script is used to 
#1.delete/extract a list of fasta entries
#2.length filter
#3.header rename
#"usage:perl Fasta_manipulator.pl --method  delete/extract --length length --fasta input.fa --lis header_list.txt"				
# when method = extract/delete/rename and you specify "length" parameter, the subset data will be length filtered
# when method= mlength, the whole set data will be length filtered, removing entries  length< $length

#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use Tie::IxHash;
my($fasta_file,$fasta_line,$header,$key,$value,$list,$list_line,$method,$length,$species,$count);
my(%fasta_db,%list);

tie %fasta_db, 'Tie::IxHash' ;
tie %list, 'Tie::IxHash' ;

GetOptions('method=s'=>\$method,
		   'length:i'=>\$length,
		   'fasta=s'=>\$fasta_file,
		   'list=s'=>\$list,
		   'species=s'=>\$species);
unless($method) 
{
	die "Please specify what you want to do: extraction or deletion ?\n",
		"usage:perl Fasta_manipulator.pl --method  delete/extract --length length --fasta input.fa --list header_list.txt"				
};


open(FASTA,"$fasta_file");
open(LIST,"$list");

foreach $list_line(<LIST>)
{
	chomp $list_line;
	$list{$list_line}=1;
}


foreach $fasta_line(<FASTA>)
{
	if($fasta_line=~/>(.+)\s+/)
	{
		chomp $fasta_line;
		$fasta_db{$1}="";
		$header=$1;	
	}
	else
	{
		chomp $fasta_line;
		$fasta_db{$header}=$fasta_db{$header}.$fasta_line;
	}
}
close FASTA;

#length filter
if($method eq 'mlength')
{
	die "please give parameter to [length] options" unless ($length);
	while (($key,$value)= each %fasta_db)
	{
		delete $fasta_db{$key} if length $value < $length;
	}
	while (($key,$value)= each %fasta_db)
	{
		print ">",$key,"\n";
		$value=$value."\n";
		print $value;
	}
}

#extraction or deletion
if($method eq 'delete' or 'extract')
{
while (($key,$value)= each %fasta_db)
{
	if($method eq 'delete')
	{
			unless($list{$key})
			{
				print ">",$key,"\n";
				$value=$value."\n";
				
				print $value;
			}
	}
	elsif($method eq 'extract')
	{
		if($list{$key})
		{
			print ">",$key,"\n";
			$value=$value."\n";
			print $value;
		}
	}
}
}

#rename
$count=1;
if($method eq 'rename')
{
	die "please give parameter to [species] options" unless ($species);
	while (($key,$value)= each %fasta_db)
	{
		print ">",$species,"_TR",$count,"\n";
		$value=$value."\n";
		$value=~s/(.{60})/$1\n/g;
		print $value;
		$count++;
	}
	
}
