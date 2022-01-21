#!/usr/bin/perl -w
use strict;
my($cluster,$key);
my(@cluster,@pos,@element,@low);
my(%hash);
# This script is the 2nd step scripts used after extract high confient LTRs and low confident LTRs(ltrs missing pbs/ppt/5nt tsd/ but has either hit with high conifdent LTR or RTase)
# After clusting high and low confident LTRs, this script will split the file to three subfiles
$cluster=shift @ARGV;
open(F,"$cluster") or die "Can't open the input file!\n";
@cluster=<F>;


for(my $i=0;$i<=$#cluster;$i++)
{
	if($cluster[$i]=~/Cluster/)
	{
		push @pos, $i;
	}
	
}
push @pos,$#cluster+1;

open(SINGLE,">1st.round.LTRepeat.high.single.list") or die $!;
open(MULTIPLE,">1st.round.LTRepeat.high.multiple.list") or die $!;

for(my $i=0;$i<$#pos;$i++)
{	
	%hash=();
	@element=();
	for(my $count=$pos[$i]+1;$count<=$pos[$i+1]-1;$count++)
	{
		print $cluster[$count];
			#may need modify
			$cluster[$count]=~/>(.+)\.\.\./;
			my $id=$1;
			$hash{$id}=1;
			#May need modify
	}

	@element=%hash; 
	if($#element<=3)
	{
		foreach $key( keys %hash)
		{
			print SINGLE $key,"\n";
 		}
	}
	elsif($#element>3)
	{
		foreach $key( keys %hash)
		{
			print MULTIPLE $key,"\n";
 		}	
	}
	
}
close SINGLE;
close MULTIPLE;













