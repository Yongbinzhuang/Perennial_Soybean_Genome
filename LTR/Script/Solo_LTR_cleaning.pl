#!/usr/bin/perl -w 
use strict;
use List::MoreUtils qw(uniq);

my($perfect_file,$mismatch_file,$noTSD_file,$line);
my(@perfect,@mismatch,@noTSD);

$perfect_file=shift @ARGV;
$mismatch_file=shift @ARGV;
$noTSD_file=shift @ARGV;

open(F1,"$perfect_file") or die $!;
open(F2,"$mismatch_file") or die $!;
open(F3,"$noTSD_file") or die $!;

@perfect=<F1>;
@mismatch=<F2>;
@noTSD=<F3>;

open(PUNIQ,">perfect_TSD.uniq.list") or die $!;
LOOP1:for(my $i=0;$i<=$#perfect;$i++)
{
	my @a=();
	if($perfect[$i])
	{
		if($perfect[$i]=~/(\w+)\s+\w+\d+\d+/)
		{
			@a=split(" ",$perfect[$i]);
			for(my $j=0;$j<=$#perfect;$j++)
			{
				my @b=();
				if($perfect[$j])
				{
					if($perfect[$j]=~/(\w+)\s+\w+\d+\d+/)
					{
						@b=split(" ",$perfect[$j]);
						if($a[0] eq $b[0] && $a[2] == $b[2] && $a[3] ==$b[3])
						{
							$perfect[$j]=$perfect[$i] ;
							next;
						}
						elsif($a[0] eq $b[0])
						{
							if($b[2]<=$a[2] && $b[3]>=$a[3])
							{
								delete $perfect[$j];
								next;
							}
							elsif($b[2]>$a[2] && $b[2]<$a[3] && $b[3] >$a[3])
							{
								if(length($perfect[$j])>=length($perfect[$i]))
								{
									delete $perfect[$j] ;
									next;
								}
								elsif(length($perfect[$j])<length($perfect[$i]))
								{
									delete $perfect[$i];
									next LOOP1;
								}
							}
						}
					}
					else
					{
						next;
					}
				}
			}
		}
	}
	else
	{
		next;
	}
}

@perfect=uniq(@perfect);
foreach $line(@perfect)
{
	print PUNIQ $line if ($line);
}
close PUNIQ;


for(my $i=0;$i<=$#perfect;$i++)
{
	my @a=();
	if($perfect[$i])
	{
		if($perfect[$i]=~/(\w+)\s+\w+\d+\d+/)
		{
			@a=split(" ",$perfect[$i]);
			for(my $j=0;$j<=$#mismatch;$j++)
			{
				my @b=();
				if($mismatch[$j])
				{
					if($mismatch[$j]=~/(\w+)\s+\w+\d+\d+/)
					{
						@b=split(" ",$mismatch[$j]);
						if($a[0] eq $b[0] && $a[2] == $b[2] && $a[3] ==$b[3])
						{
							delete $mismatch[$j] ;
							next;
						}
						elsif($a[0] eq $b[0])
						{
							if($b[2]<=$a[2] && $b[3]>=$a[3])
							{
								delete $mismatch[$j];
								next;
							}
							elsif($b[2]>$a[2] && $b[2]<$a[3] && $b[3] >$a[3])
							{
								delete $mismatch[$j];
								next;
							}
							elsif($b[2]<$a[2] && $b[3]<$a[3] && $b[3] >$a[2])
							{
								delete $mismatch[$j];
								next;
							}
						}
					}
					else
					{
						next;
					}
				}
			}
		}
	}
	else
	{
		next;
	}
}

@mismatch=uniq(@mismatch);

LOOP2:for(my $i=0;$i<=$#mismatch;$i++)
{
	my @a=();
	if($mismatch[$i])
	{
		if($mismatch[$i]=~/(\w+)\s+\w+\d+\d+/)
		{
			@a=split(" ",$mismatch[$i]);
			for(my $j=0;$j<=$#mismatch;$j++)
			{
				my @b=();
				if($mismatch[$j])
				{
					if($mismatch[$j]=~/(\w+)\s+\w+\d+\d+/)
					{
						@b=split(" ",$mismatch[$j]);
						if($a[0] eq $b[0] && $a[2] == $b[2] && $a[3] ==$b[3])
						{
							$mismatch[$j]=$mismatch[$i] ;
							next;
						}
						elsif($a[0] eq $b[0])
						{
							if($b[2]<=$a[2] && $b[3]>=$a[3])
							{
								delete $mismatch[$j];
								next;
							}
							elsif($b[2]>$a[2] && $b[2]<$a[3] && $b[3] >$a[3])
							{
								if(length($mismatch[$j])>=length($mismatch[$i]))
								{
									delete $mismatch[$j] ;
									next;
								}
								elsif(length($mismatch[$j])<length($mismatch[$i]))
								{
									delete $mismatch[$i];
									next LOOP2;
								}
							}
						}
					}
					else
					{
						next;
					}
				}
			}
		}
	}
	else
	{
		next;
	}
}

@mismatch=uniq(@mismatch);
open(MUNIQ,">1mismatch_TSD.uniq.list") or die $!;

foreach $line(@mismatch)
{
	print MUNIQ $line if ($line);
}
close MUNIQ;


LOOP3:for(my $i=0;$i<=$#noTSD;$i++)
{
	my @a=();
	if($noTSD[$i])
	{
		if($noTSD[$i]=~/(\w+)\s+\w+\d+\d+/)
		{
			@a=split(" ",$noTSD[$i]);
			for(my $j=0;$j<=$#noTSD;$j++)
			{
				my @b=();
				if($noTSD[$j])
				{
					if($noTSD[$j]=~/(\w+)\s+\w+\d+\d+/)
					{
						@b=split(" ",$noTSD[$j]);
						if($a[0] eq $b[0] && $a[2] == $b[2] && $a[3] ==$b[3])
						{
							$noTSD[$j]=$noTSD[$i] ;
							next;
						}
						elsif($a[0] eq $b[0])
						{
							if($b[2]<=$a[2] && $b[3]>=$a[3])
							{
								delete $noTSD[$j];
								next;
							}
							elsif($b[2]>$a[2] && $b[2]<$a[3] && $b[3] >$a[3])
							{
								if(length($noTSD[$j])>=length($noTSD[$i]))
								{
									delete $noTSD[$j] ;
									next;
								}
								elsif(length($noTSD[$j])<length($noTSD[$i]))
								{
									delete $noTSD[$i];
									next LOOP3;
								}
							}
						}
					}
					else
					{
						next;
					}
				}
			}
		}
	}
	else
	{
		next;
	}
}

@noTSD=uniq(@noTSD);
open(NUNIQ,">noTSD_TSD.uniq.list") or die $!;

foreach $line(@noTSD)
{
	print NUNIQ $line if ($line);
}
close NUNIQ;


























