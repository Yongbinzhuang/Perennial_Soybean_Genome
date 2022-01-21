#!/usr/bin/perl -w
use strict;
my($blast_out,$line,@perfect,$line2);
my(@blast);
my(%ele);

$blast_out=shift @ARGV;
open(F,"$blast_out") or die $!;
@blast=<F>;
close F;

foreach $line(@blast)
{
	chomp $line;
	my @temp=split(" ",$line);
	$temp[1]=~/(.+?)_/;
	my $ref=$1;
	$temp[0]=~/(.+?)_/;
	my $length;
	if($1 eq $ref && $temp[2] ==100 && $temp[10]==0)
	{
		$temp[0]=~/(\d+)_(\d+)/;
		$length=$2-$1+1;
		$line=$line."\t".$length;
		push @perfect,$line;
	}
}

foreach $line(@perfect)
{
	$line=~/(.+?)\s+/;
	if($ele{$1})
	{
		$ele{$1}=$ele{$1}.":".$line;
	}
	else
	{
		$ele{$1}=$line;
	}
}

foreach my $key(keys %ele)
{
	open(OUTPUT,">>Nested_LTR_coordinate.txt") or die $!;
	my @array=();
	my @array2=();
	my @posA=();
	my @posB=();
	my @subA=();
	my @subB=();
	my %out=();
	my @pos;
	my $ltr='';
	my $genome='';
	my $sw=0;
	my $target=1;
	my @ltr=split(":",$ele{$key});

	foreach $line(@ltr)
	{
		@array=split(" ",$line);
		push @posA,$array[6]."_".$array[7];		
		push @posB,$array[8]."_".$array[9];
		$out{$line}=$array[6];
	}
	
	open(OUT5,">>blast.nestedLTR2genome.perfect.match.out") or die;
	foreach my $line3(sort{$out{$a}<=>$out{$b}} keys %out)
	{
		print OUT5 $line3,"\n";
	}

	foreach $line(@posA)
	{
		@array2=split("_",$line);
		push @subA,$array2[0];
		push @subB,$array2[1];
	}
	LOOP1:for(my $i=0;$i<=$#subA;$i++)
	{
		if($subA[$i]==$target)
		{
			$target=$subB[$i]+1;
			push @pos,$subA[$i]."_".$subB[$i];
			$i=-1;
			$sw++;
			next LOOP1;
		}
	}
	for(my $i=0;$i<=$#posA;$i++)
	{
		foreach $line(@pos)
		{
			if($posA[$i] eq $line)
			{
				if($ltr)
				{
					$ltr=$ltr."_".$line;
				}
				else
				{
					$ltr=$line;
				}
				if($genome)
				{
					$genome=$genome."_".$posB[$i];
				}
				else
				{
					$genome=$posB[$i];
				}
			}
		}
	}
	if($sw==0)
	{
		print OUTPUT $key, "\tStarting Sites Not Found,\n" ;
	}
	else
	{
		my @out1=split("_",$ltr);
		@out1=sort{$a<=>$b} @out1;
		my @out2=split("_",$genome);
		@out2=sort{$a<=>$b} @out2;
		my $ltr_mapped_length=$out1[$#out1]-$out1[0]+1;
		if($ltr_mapped_length == $array[12])
		{
			print OUTPUT $key,"\t",join("\t",@array[1,12]),"\t","LTR_matched:",join(",",@out1),"\t","Genome_coordinate:",join(",",@out2),"\n";
		}
		else
		{
			print OUTPUT $key,"\t",join("\t",@array[1,12]),"\t","LTR_matched:",join(",",@out1),"\t","Genome_coordinate:",join(",",@out2),"\t","Partial","\n";
	
		}
	}
	close OUT5;
	close OUTPUT;
}















