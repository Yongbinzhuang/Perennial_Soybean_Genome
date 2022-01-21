#!/usr/bin/perl -w
use strict;
my $file=shift @ARGV;

open(F,"$file") or die $!;

my @file=<F>;

foreach my $line(@file)
{
        my @temp=split(" ",$line);
        my %hash;
        for(my $i=0;$i<=$#temp-1;$i++)
        {
                for(my $j=1;$j<=$#temp;$j++)
                {
                        my $id=$temp[$i].$temp[$j];
                        my $di=$temp[$j].$temp[$i];
                        unless($temp[$i] eq "\."||$temp[$j] eq "\."||$temp[$i] eq $temp[$j]||$hash{$id}||$hash{$di})
                        {
                                open(OUT,">$temp[$i].$temp[$j].list");
                                print OUT $temp[$i],"\n";
                                print OUT $temp[$j],"\n";
                                $hash{$id}=1;
                                $hash{$di}=1;
                        }
                }
        }
        my @list=glob("G*list");
        foreach my $list(@list)
        {
         `seqtk subseq Added.aa.fasta $list >$list.aa.fa`;
        `seqtk subseq Added.cds.fasta $list >$list.cds.fa`;
        `mafft --auto $list.aa.fa >$list.aa.align`;
        `~/Soft/pal2nal.v14/pal2nal.pl -nogap $list.aa.align $list.cds.fa >$list.pal2nal`;
         `~/Soft/KaKs_Calculator2.0/src/AXTConvertor $list.pal2nal $list.axt`;
        `~/Soft/KaKs_Calculator2.0/bin/Linux/KaKs_Calculator -i $list.axt -o $list.kaks -m YN`;
        `cat *kaks >>$temp[0].w`;
        unlink "$list.kaks";
        unlink "$list";
        }
}
