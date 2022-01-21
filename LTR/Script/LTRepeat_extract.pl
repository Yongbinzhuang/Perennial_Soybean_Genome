#!/usr/bin/perl -w
my($gFasta,$ltr_file,$ltr_line,$recorder);
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


print("Seq_ID","\t","LTR_Length","\t","LTR_end5_Length","\t","end5_start","\t","end5_end","\t","LTR_end5_seq","\t","LTR_end3_Length","\t","end3_start","\t","end3_end","\t","LTR_end3_seq","\t","TG?","\t","CA?","\t","TSR","\t","Strand","\t","TSR1_start","\t","TSR1_end","\t","TSR2_start","\t","TSR2_end","\t","PBS_start","\t","PBS_end","\t","PPT_start","\t","PPT_end","\t","PBS","\t","PPT","\n");

#read in ltr_finder output and rextract the seq
for(my $count=0;$count<$#ltr;$count++)
{
	$ltr_line=$ltr[$count];
	chomp $ltr_line;
	if($ltr_line=~/arrow_pilon_arrow_pilon/)
	{
		my @out=split(" ",$ltr_line);
		my $ID=$out[1];
		$header=$out[1];
		print "\n";
	}
	elsif($ltr_line=~/Location/)
	{
		my @out=split(" ",$ltr_line);
		print $header,"_",$out[2],"_",$out[4],"\t",$out[6],"\t";
	}
	elsif($ltr_line=~/^5.+LTR/)
	{
		my @out=split(" ",$ltr_line);
		my $start=$out[2]-1;
		my $length=$out[4]-$start;
		my $seq=substr($gSeq{$header},$start,$length);
		print $out[6],"\t",$out[2]-1,"\t",$out[4],"\t",$seq,"\t";
	}
	elsif($ltr_line=~/^3.+LTR/)
	{
		my @out=split(" ",$ltr_line);
		my $start=$out[2]-1;
		my $length=$out[4]-$start;
		my $seq=substr($gSeq{$header},$start,$length);
		print $out[6],"\t",$out[2]-1,"\t",$out[4],"\t",$seq,"\t";
	}
	elsif($ltr_line=~/Sharpness/)
	{
		my @t1=split(" ",$ltr[$count+1]);
		print $t1[1],"\t";
		if($ltr[$count-1]=~/^TSR.+\d+/)
		{
			my @t3=split(" ",$ltr[$count-1]);
			print join("\t",@t3[2,4,,6,8]),"\t";
		}
		else
		{
			print "NA","\t","NA","\t","NA","\t","NA","\t";;
		}
		if($ltr[$count+2]=~/PBS/)
		{
			my @t2=split(" ",$ltr[$count+2]);
			print $t2[3],"\t",$t2[5],"\t";
		}
		else
		{
			print "NA","\t","NA","\t";
		}
		if($ltr[$count+3]=~/PPT/)
		{
			my @t4=split(" ",$ltr[$count+3]);
			print $t4[3],"\t",$t4[5],"\t";
		}
		elsif($ltr[$count+2]=~/PPT/)
		{
			my @t4=split(" ",$ltr[$count+2]);
			print $t4[3],"\t",$t4[5],"\t";
		}
		else
		{
			print "NA","\t","NA","\t";

		}
	}

	elsif($ltr_line=~/^5.+TG/)
	{
		my @out=split(":",$ltr_line);
		print $out[1],"\t";
	}
	elsif($ltr_line=~/^3.+CA/)
	{
		my @out=split(":",$ltr_line);
		print $out[1],"\t";
	}
	elsif($ltr_line=~/^TSR/)
	{
		if($ltr_line=~/FOUND/)
		{
			print "NA","\t";
		}
		else
		{ 
			$ltr_line=~/\[(.+)\]/;
			print $1,"\t";
		}
	}
	elsif($ltr_line=~/^tRNA/)
	{
		chomp $ltr[$count+3];
		print $ltr[$count+3],"\t";
		$recorder=1;
	}
	elsif($ltr_line=~/^Details\s+of\s+PPT/)
	{
		if($recorder==1)
		{
			chomp $ltr[$count+1];
			print $ltr[$count+1];
			$recorder=0;
		}
		elsif($recorder==0)
		{
			chomp $ltr[$count+1];
			print "NA","\t",$ltr[$count+1];
		}
	}

}










