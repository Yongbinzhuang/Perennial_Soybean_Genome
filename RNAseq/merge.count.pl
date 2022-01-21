#!/usr/bin/perl -w
my @file=glob("*count");
open(F,"$file[0]") or die $!;
my @file2=<F>;
close F;
my %gene;

for(my $i=0;$i<=$#file2;$i++)
{
	my @temp=split(" ",$file2[$i]);
	$gene{$temp[0]}=$i;
}


open(OUT,">merged.count.tables2") or die $!;

print OUT "gene_id\t";
foreach my $count2(@file)
{
	print OUT $count2,"\t";
}
print OUT "\n";

foreach my $key(keys %gene)
{
	print OUT "$key\t";
	foreach my $temp_file(@file)
	{
		open(F,"$temp_file") or die $!;
		my @tempfile=<F>;
		close F;
		my @temp=split(" ",$tempfile[$gene{$key}]);
		print OUT $temp[1],"\t";
	}
	print OUT "\n";
}

