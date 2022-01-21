#!/usr/bin/perl -w
use strict;
my($reflist,$candilist,$line,$ltr,$ltr2,$line2);
my(@ref,@target,@uniq);
my(%hash);

$reflist=shift @ARGV;
$candilist=shift @ARGV;

open(REF,"$reflist") or die $!;
@ref=<REF>;

open(CAND,"$candilist") or die $!;
@target=<CAND>;
open(OUT2,">1st.round.LTRepeat.low2high.nested.removed.has.RTase.list") or die $!;


foreach $line(@target)
{
  my $code=0;
  $line=~/(.+?)_[FR]_(\d+)_(\d+)/;
  my $id=$1;
  my $pos1=$2;
  my $pos2=$3;
  foreach $ltr(@ref)
  {
    $ltr=~/(.+?)_[FR]_(\d+)_(\d+)/;
    if($1 eq $id)
    {
      if($3<$pos1 || $2>$pos2)
      {
        next
      }
      else
      {
        $code=1;
        $ltr2=$ltr;
        last;
      }
    }
  }
  print OUT2 $line if $code==0;
}
close OUT2;

