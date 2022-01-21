#!/usr/bin/perl -w
use strict;
my($cluster,$on);
my(@cluster,@pos,@element,@low);
my(%hash);
# This script is the 2nd step scripts used after extract high confient LTRs and low confident LTRs(ltrs missing pbs/ppt/5nt tsd/ but has either hit with high conifdent LTR or RTase)
# After clusting high and low confident LTRs, this script will split the file to three subfiles
$cluster=shift @ARGV;
open(F,"$cluster") or die "Can't open the input file!\n";
@cluster=<F>;

open(OUT1,">1st.low2high.has.hit") or die "Can't open the output file1,\n";
open(OUT2,">1st.low2high.no.hit.multiple") or die "Can't open the output file1,\n";
open(OUT3,">1st.low2high.no.hit.singleton") or die "Can't open the output file1,\n";

for(my $i=0;$i<=$#cluster;$i++)
{
	if($cluster[$i]=~/Cluster/)
	{
		push @pos, $i;
	}
}

for(my $i=0;$i<$#pos;$i++)
{	
	%hash=();
	@element=();
	@low=();
	$on=0;
	for(my $count=$pos[$i]+1;$count<=$pos[$i+1]-1;$count++)
	{
			#may need modify
			$cluster[$count]=~/>(.+)_(.+)_(\d+)_(\d+)/;
			my $id=$1;
			my $id2=$2;
			my $ps=$3;
			my $pe=$4;
			$hash{$id}=$1."_".$2."_".$3."_".$4;
			#May need modify
			if($cluster[$count]=~/low/)
			{
				push @low,$hash{$id};
			}
			elsif($cluster[$count]=~/>0/)
			{
				$on=1;
			}
	}
	@element=%hash;
	if($#element==1 && $on==0)
	{
		foreach my $line(@low)
		{
			print OUT3 $line,"\n";
		}
	}
	elsif($#element>2 && $on==0)
	{
		foreach my $line(@low)
		{
			print OUT2 $line,"\n";
		}
	}
	elsif($#element>2 && $on==1)
	{
		foreach my $line(@low)
		{
			print OUT1 $line,"\n";
		}
	}
}
close OUT1;
close OUT2;
close OUT3;













