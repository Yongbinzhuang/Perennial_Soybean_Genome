#!/usr/bin/perl -w
my($gFasta,$ltr_file,$ltr_line,$recorder,$gstart,$gend,$star5,$end5,$star3,$end3);
my(@genome,@ltr);
my(%gSeq);

$gFasta=shift @ARGV;
$ltr_file=shift @ARGV;

open(GENOME,"$gFasta") or die $!;
open(LTR,"$ltr_file") or die $!;

@genome=<GENOME>;
@ltr=<LTR>;

#read in all the fasta files in hash;
for(my $count=0;$count<$#genome;$count++)
{
	chomp $genome[$count];
	if($genome[$count]=~/>(.+)/)
	{
		my $key=$1;
		chomp $genome[$count+1];
		$gSeq{$key}=$genome[$count+1];
	}
}



#read in ltr_finder output and rextract the seq
for(my $count=0;$count<$#ltr;$count++)
{
	$ltr_line=$ltr[$count];
	chomp $ltr_line;
	if($ltr_line=~/Location/)
	{
		$ltr[$count-1]=~/\s+(.+pilon)\s+/;
		$header=$1;
		my @out=split(" ",$ltr_line);
		my $start=$out[2]-1;
		my $length=$out[4]-$start;
		my $seq=substr($gSeq{$header},$start,$length);
		print ">",$header,"_",$out[2]."_",$out[4],"\n";
	}
	elsif($ltr_line=~/^5.+LTR/)
	{
		my @out=split(" ",$ltr_line);
		$gstart=$out[4];
	}
	elsif($ltr_line=~/^3.+LTR/)
	{
		my @out=split(" ",$ltr_line);
		my $length=$out[2]-$gstart;
		my $seq=substr($gSeq{$header},$gstart,$length-1);
		print $seq,"\n";
	}
	
}










