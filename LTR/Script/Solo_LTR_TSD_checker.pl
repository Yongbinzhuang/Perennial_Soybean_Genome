#!/usr/bin/perl -w
#This will check 10bp inside plus 50bp outside the original blast output boarder,can be modified accordly 
use strict;
use Getopt::Long;
use Tie::IxHash;

my($line,$genome,$solo,$length);
my(@list,@genome);
my(%ghash);

$genome=shift @ARGV;
$solo=shift @ARGV;
GetOptions('length:i'=>\$length);
unless($length)
{
	$length=50;
}

open(F,"$genome") or die $!;
@genome=<F>;
for(my $i=0;$i<$#genome;$i++)
{
	if($genome[$i]=~/>(.+)/)
	{
		$ghash{$1}=$genome[$i+1]
	}
}

open(F2,"$solo");
@list=<F2>;

foreach $line(@list)
{
	$line=~/(.+)\s+(.+)\s+(\d+)\s+(\d+)/;
	my $left=substr($ghash{$1},$3-$length+1,$length+10);
	my $right=substr($ghash{$1},$4-10,$length+10);
	my $header=$1;
	my $family=$2;
	my $oleft=$3;
	my $oright=$4;
	my @posTG;
	my @posCA;
	my @TGs=();
	my @CAs=();
	my $posi;
	my $temp;
	push @posTG,pos($left)-2 while $left=~/TG/g;
	push @posCA,pos($right)-2 while $right=~/CA/g;
	foreach $posi(@posTG)
	{
		$temp=substr($left,$posi-5,5);
		push @TGs,$temp ;
	}
	foreach $posi(@posCA)
	{
		$temp=substr($right,$posi+2,5);
		push @CAs,$temp ;
	}
	chomp $line;
	my $recorder=0;
	#print $line,"\t",$left,"\t",$right,"\n";
	for(my $count1=0;$count1<=$#TGs;$count1++)
	{
		for(my $count2=0;$count2<=$#CAs;$count2++)
		{
			my @mafft=();
			my @match=();
			open(OUT,">temp.fasta") or die $!;
			print OUT ">forward\n",$TGs[$count1],"\n";
			print OUT ">reverse\n",$CAs[$count2],"\n";
			@mafft= `mafft --quiet --clustalout --maxiterate 1000 --localpair temp.fasta`;
			@match=$mafft[$#mafft]=~/\*/g;
			#add blast comparison with matched intact LTRepeats!!!
			#xxxx
			#xxxx
			#xxx
			if($#match==4)
			{
				open(PERFECT,">>perfect.TSD.list") or die $!;
				print PERFECT $header,"\t",$family,"\t",$oleft+$posTG[$count1]-$length+2,"\t",$oright+$posCA[$count2]-8,"\n";
				close PERFECT;
				$recorder=1;
			}
			elsif($#match==3)
			{
				open(MIS,">>1mismatch.TSD.list") or die $!;
				print MIS $header,"\t",$family,"\t",$oleft+$posTG[$count1]-$length+2,"\t",$oright+$posCA[$count2]-8,"\n";
				close MIS;
				$recorder=1;
			}
			close OUT;
		}
	}
	if($recorder==0)
	{
		open(NO,">>No_TSD_found.list") or die $!;
		print NO $line,"\n";
		close NO;
	}
}
unlink "temp.fasta";


