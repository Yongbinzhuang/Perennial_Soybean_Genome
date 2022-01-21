#!/usr/bin/perl -w
use strict;
use List::MoreUtils qw(uniq);

my($fasta,$line1,$line2,$pcount);
my(@afasta);
my(%hfasta,%cluster,@cluster);

$fasta=shift @ARGV;

`makeblastdb -in $fasta -input_type fasta -dbtype nucl -out ../BlastDB/LTR_DB`;

open(F,"$fasta") or die $!;
@afasta=<F>;

for(my $i=0;$i<$#afasta;$i++)
{
	if($afasta[$i]=~/>(.+)/)
	{
		$hfasta{$1}=$afasta[$i+1];
	}
}
my @key=keys %hfasta;
my $group=0;
foreach $line1(@key)
{	
	@cluster=();
	open(TEMP,">temp.fasta") or die $!;
	print TEMP ">",$line1,"\n",$hfasta{$line1},"\n";
	`blastn -db ../BlastDB/LTR_DB -out temp.out -query temp.fasta -num_threads 60 -evalue 1e-10 -outfmt 6`;
	open(TEMP2,"temp.out") or die $!;
	my @blast=<TEMP2>;
	foreach my $line2(@blast)
	{
		my @fields=split(" ",$line2);
		my $length1=length($hfasta{$fields[0]});
		my $length2=length($hfasta{$fields[1]});
		push @fields,$length1;
		push @fields,$length2;
		if ($fields[2]>=80 && $fields[3]/$fields[12]>=0.8)
		{ #$fields[3]/$fields[13]>=0.8 add this line to if test if need 80 coverage for both seqs
			push @cluster,$fields[0];
			push @cluster,$fields[1];	
		}
	}
	@cluster=uniq(@cluster);
	foreach my $ltr(@cluster)
	{
		if(defined $cluster{$ltr})
		{
			foreach my $ltr2(@cluster)
			{
				$cluster{$ltr2}=$cluster{$ltr};
			}
			last;
		}
		else
		{
			$group++;
			foreach my $ltr2(@cluster)
			{
				$cluster{$ltr2}=$group;
			}
		}
	}
}



my %counter;
foreach my $key(keys %cluster)
{
	if($counter{$cluster{$key}})
	{
		$counter{$cluster{$key}}++;
	}
	else
	{
		$counter{$cluster{$key}}=1;
	}
}


open(OUT,">temp.file") or die $!;

foreach my $key(sort {$counter{$b}<=>$counter{$a}} keys %counter)
{
	foreach my $key2(keys %cluster)
	{
		print OUT $key2,"\t",$cluster{$key2},"\n" if $cluster{$key2}==$key;
	}
}
close OUT;

my $count=1;
my %unique;
open(F,"temp.file") or die $!;
open(OUT,">Final.LTRepeat.family.txt") or die $!;
my @outfile=<F>;
for(my $i=0;$i<=$#outfile;$i++)
{
	my @tempin=split(" ",$outfile[$i]);
	$tempin[0]=~/(.+)_[FR]_(\d+)_(\d+)/;
	my $tt=$1."_".$2."_".$3;
	if($i==0)
	{
		unless(defined $unique{$tt})
		{
			print OUT $tempin[0],"\t","Family",$count,"\n";
			$pcount=$tempin[1];
			$unique{$tt}=1;
			next;
		}
	}
	else
	{
		if($tempin[1] == $pcount)
		{
			unless(defined $unique{$tt})
			{ 
				print OUT $tempin[0],"\t","Family",$count,"\n";
				$unique{$tt}=1;
			}
		}
		else
		{	
			unless(defined $unique{$tt})
			{
				$count++;
				print OUT $tempin[0],"\t","Family",$count,"\n";
				$pcount=$tempin[1];
				$unique{$tt}=1;
			}

		}
	}
} 
unlink "temp.file";















