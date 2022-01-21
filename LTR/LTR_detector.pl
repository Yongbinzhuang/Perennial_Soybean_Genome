#!/usr/bin/perl -w
use strict;
use Cwd;
use Array::Utils qw(:all);
use List::MoreUtils qw(uniq);
my($method,$LTR_length,$config,$num_threads,$fasta_file,$list,$line,$count,$LTR_finder,$LTR_detector,$cd_hit,$blast);
my($key,$RTase,$rounds,$job, $temp);
my(@fasta,@config,@LTR_repeat,@blastx,@genome,@blastn,@temp,@temp2,@nested);
my(%temp_hash,%gSeq);

$config=shift @ARGV;
unless(defined $config)
{
	print "You must provide config file containing configuration info for LTR identification!\n";
}

# readin config file
open(CONFIG,"$config") or die $!;
@config=<CONFIG>;
foreach $line(@config)
{
	$line=~/(.+)\=(.+)/;
	if($1 eq "fasta_file")
	{
		$fasta_file=$2
	}
	elsif($1 eq "LTR_finder")
	{
		$LTR_finder=$2
	}
	elsif($1 eq "LTR_detector")
	{
		$LTR_detector=$2
	}
	elsif($1 eq "method")
	{
		$method=$2;
	}
	elsif($1 eq "LTR_length")
	{
		$LTR_length=$2
	}
	elsif($1 eq "cd-hit")
	{
		$cd_hit=$2;
	}
	elsif($1 eq "blast")
	{
		$blast=$2;
	}
	elsif($1 eq "RTase")
	{
		$RTase=$2;
	}
	elsif($1 eq "num_threads")
	{
		$num_threads=$2;
	}
	elsif($1 eq "rounds")
	{
		$rounds=$2
	}
}


# Build temp folder and rename fasta file header
my $dir=getcwd;
if( -d "Temp" )
{
	system("rm -rf Temp");
	mkdir "Temp";
}
else
{
	mkdir "Temp";
}

open(FASTA,"$fasta_file") or die $!;
open(SINGLE,">Temp/formated.fasta") or die $!;
open(DICTION,">Temp/dictionary.list") or die $!;
@fasta=<FASTA>;

$count=0;
for(my $i=0;$i<=$#fasta;$i++)
{
	if($fasta[$i]=~/>(.+)/)
	{
		if($i==0)
		{
			$count++;
			print SINGLE ">","000",$count,"F_arrow_pilon_arrow_pilon","\n";
			print DICTION  $1,"\t","000",$count,"F_arrow_pilon_arrow_pilon\n";
		}
		else
		{
			$count++;
			print SINGLE "\n",">","000",$count,"F_arrow_pilon_arrow_pilon","\n";
			print DICTION  $1,"\t","000",$count,"F_arrow_pilon_arrow_pilon\n";			
		}
	}
	else
	{
		chomp $fasta[$i];
		print SINGLE $fasta[$i];
	}
}
close FASTA;
close SINGLE;
close DICTION;


#readin formated fasta file and start LTR identification
chdir "$dir/Temp";
if($method eq "LTR_finder")
{
	`$LTR_finder/ltr_finder  -s $LTR_finder/tRNAdb/Athal-tRNAs.fa formated.fasta >1st.round.LTR.out`;
	`$LTR_detector/Script/LTRepeat_extract.pl formated.fasta 1st.round.LTR.out >1st.round.LTRepeat.txt`;
	`$LTR_detector/Script/LTR_full_extract.pl formated.fasta 1st.round.LTR.out >1st.round.LTR.full.fasta`;
	`$LTR_detector/Script/LTR_internal_extract.pl formated.fasta 1st.round.LTR.out >1st.round.LTR.internal.fasta`;
}

open(F,"1st.round.LTRepeat.txt") or die $!;
open(OUT1,">1st.round.LTRepeat.fasta") or die $!;
@LTR_repeat=<F>;
foreach $line(@LTR_repeat)
{	unless($line=~/^$/)
	{
		my @F=split(" ",$line);
		if($F[0]=~/(.+)\_arrow\_pilon\_arrow\_pilon\_(\d+)_(\d+)/)
		{
			print OUT1 ">",$1,"_F_",$2,"_",$3,"\n",$F[5],"\n" if length($F[5]) >=$LTR_length && length($F[9]) >=$LTR_length;
			print OUT1 ">",$1,"_R_",$2,"_",$3,"\n",$F[9],"\n" if length($F[5]) >=$LTR_length && length($F[9]) >=$LTR_length;
		}
	}
}
close F;
close OUT1;


#Extract high confident LTRepeat list
open(F,"1st.round.LTRepeat.txt") or die $!;
open(OUT1,">1st.round.LTRepeat.high.list") or die $!;
@LTR_repeat=<F>;
foreach $line(@LTR_repeat)
{
	unless($line=~/^$/)
	{
		unless($line=~/NA/)
		{
			my @F=split(" ",$line);
			if(length($F[16])==5)
			{
				$F[0]=~/(.+?)_arrow_pilon_arrow_pilon_(\d+)_(\d+)/;
				print OUT1 $1,"_F_",$2,"_",$3,"\n";
	 			print OUT1 $1,"_R_",$2,"_",$3,"\n";
			}
		}
	}
}
close F;
close OUT1;

#3 split fasta file from #4 into two sets of subset data
`perl $LTR_detector/Script/Fasta_manipulate.pl --fasta 1st.round.LTRepeat.fasta --list 1st.round.LTRepeat.high.list --method extract >1st.round.LTRepeat.high.fasta`;
`perl $LTR_detector/Script/Fasta_manipulate.pl --fasta 1st.round.LTRepeat.fasta --list 1st.round.LTRepeat.high.list --method delete >1st.round.LTRepeat.low.fasta`;

open(F,"1st.round.LTRepeat.low.fasta") or die $!;
open(LOW,">1st.round.LTRepeat.low.temp") or die $!;
@LTR_repeat=<F>;
foreach $line(@LTR_repeat)
{
	if($line=~/>(.+)/)
	{
		print LOW $1,"\n";
	}
}
close F;
close LOW;

open(F,"1st.round.LTRepeat.high.fasta") or die $!;
open(HIGHT,">1st.round.LTRepeat.high.temp") or die $!;
@LTR_repeat=<F>;
foreach $line(@LTR_repeat)
{
	if($line=~/>(.+)/)
	{
		print HIGHT $1,"\n";
	}
}
close F;
close HIGHT;

my($reflist,$candilist,$ltr,$ltr2,$line2);
my(@ref,@target,@uniq);

open(HIGHT,"1st.round.LTRepeat.high.temp") or die $!;
open(LOW,"1st.round.LTRepeat.low.temp") or die $!;
open(LOW2,">1st.round.LTRepeat.low.temp2") or die $!;
@ref=<HIGHT>;
@target=<LOW>;

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
      if($3<=$pos1 || $2>=$pos2)
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
  print LOW2 $line if $code==0;
}
close LOW2;


`perl  $LTR_detector/Script/Fasta_manipulate.pl --method extract --list 1st.round.LTRepeat.low.temp2 --fasta 1st.round.LTRepeat.fasta >1st.round.LTRepeat.low.temp.fasta`;
open(F,"1st.round.LTRepeat.low.temp.fasta");
@genome=<F>;
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

foreach $key(keys %gSeq)
{
	if($key=~/(\w+)_F_(.+)/)
	{
		my $seqid=$1."_R_".$2;
		open(SCORE,">temp.score.fasta") or die $!;
		print SCORE ">F\n",$gSeq{$key},"\n";
		print SCORE ">R\n",$gSeq{$seqid},"\n";
		my $seq_length=length($gSeq{$key});
		my @align1=`mafft --quiet --maxiterate 1000 --clustalout --localpair --thread 1 temp.score.fasta`;
		my @match=();
		my $sum=0;
		my $similarity=0;
		foreach my $alignment(@align1)
		{
			if($alignment=~/\*/)
			{
				@match=$alignment=~/\*/g;
				$sum=$sum+$#match;
			}

		}
		$similarity=($sum/$seq_length)*log_N($seq_length,10);
		open(MARK,">>1st.round.LTRepeat.low.temp3") or die $!;
		print MARK $key,"_",$similarity,"\n";
		close MARK;
		close SCORE;
	}
}
%gSeq=();

open(PUNIQ,">1st.round.LTRepeat.low.list") or die $!;
open(CLEAN,"1st.round.LTRepeat.low.temp3") or die $!;
my @perfect=<CLEAN>;
@perfect=sort{$a cmp $b} @perfect;
LOOP1:for(my $i=0;$i<=$#perfect;$i++)
{
	my @a=();
	if($perfect[$i])
	{
		if($perfect[$i]=~/(\w+)_F_(\d+)_(\d+)_(\d+)/)
		{
			@a=split("_",$perfect[$i]);
			for(my $j=0;$j<=$#perfect;$j++)
			{
				my @b=();
				if($perfect[$j])
				{
					if($perfect[$j]=~/(\w+)_F_(\d+)_(\d+)_(\d+)/)
					{
						@b=split("_",$perfect[$j]);
						if($a[0] eq $b[0] && $a[2] == $b[2] && $a[3] ==$b[3])
						{
							$perfect[$j]=$perfect[$i] ;
							next;
						}
						elsif($a[0] eq $b[0])
						{
							if($b[2]<=$a[2] && $b[3]>=$a[3])
							{
								next;
							}
							elsif($b[2]>=$a[2] && $b[2]<=$a[3] && $b[3] >=$a[3])
							{
								if($a[4]>=$b[4])
								{
									delete $perfect[$j] ;
									next;
								}
								elsif($a[4]<$b[4])
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
	if(defined $line && $line=~/(\w+)_F_(.+)_(\d+)/)
	{
		print PUNIQ $1,"_F_",$2,"\n";
		print PUNIQ $1,"_R_",$2,"\n";
	}
}
close PUNIQ;

`$cd_hit/cd-hit-est -i 1st.round.LTRepeat.high.fasta -o 1st.round.LTRepeat  -c 0.8 -n 5 -d 0 -M 80000 -T 50 -aS 0.8 -aL 0.8 -sc 1 -b 40`;
`perl  $LTR_detector/Script/LTRepeat_high_confi_cluster_split.pl 1st.round.LTRepeat.clstr`;
`perl  $LTR_detector/Script/Fasta_manipulate.pl --method extract --list 1st.round.LTRepeat.high.multiple.list --fasta 1st.round.LTRepeat.fasta >1st.round.LTRepeat.high.M.fasta`;
`cat 1st.round.LTRepeat.high.single.list 1st.round.LTRepeat.low.list >1st.round.LTRepeat.high_single.and.L.list`;
`perl  $LTR_detector/Script/Fasta_manipulate.pl --method extract --list 1st.round.LTRepeat.high_single.and.L.list --fasta 1st.round.LTRepeat.fasta >1st.round.LTRepeat.high_single.and.L.fasta`;

open(F,"1st.round.LTRepeat.high_single.and.L.fasta") or die $!;
open(OUT,">1st.round.LTRepeat.high_single.and.L.renamed.fasta") or die $!;
@LTR_repeat=<F>;
foreach $line(@LTR_repeat)
{
	if($line=~/>/)
	{
		$line=~s/>/>low/;
		print OUT $line;
	}
	else
	{
		print OUT $line;
	}
	
}
close F;
close OUT;

`cat 1st.round.LTRepeat.high.M.fasta 1st.round.LTRepeat.high_single.and.L.renamed.fasta >1st.round.LTRepeat.all.low.renamed.fasta`;
`$cd_hit/cd-hit-est -i 1st.round.LTRepeat.all.low.renamed.fasta -o 1st.round.LTRepeat.low2high  -c 0.8 -n 5 -d 0 -M 80000 -T 50 -aS 0.8 -aL 0.8 -sc 1 -b 40`;
`perl  $LTR_detector/Script/LTRepeat_low2high_cluster.process.pl 1st.round.LTRepeat.low2high.clstr`;
`cat 1st.low2high.no.hit.multiple 1st.low2high.no.hit.singleton >1st.low2high.nohit.list`;

open(F,"1st.round.LTRepeat.high.multiple.list") or die $!;
open(OUT,">1st.round.LTRepeat.high.multiple.long.list") or die $!;
@LTR_repeat=<F>;
%temp_hash=();
foreach $line(@LTR_repeat)
{
	$line=~/(.+)\_[FR]_(\d+)_(\d+)/;
	my $temp=$1."_arrow_pilon_arrow_pilon_".$2."_".$3;
	$temp_hash{$temp}=1;
}
foreach $key(keys %temp_hash)
{
	print OUT $key,"\n";
}

close F;
close OUT;
`perl  $LTR_detector/Script/Fasta_manipulate.pl --method extract --list 1st.round.LTRepeat.high.multiple.long.list --fasta 1st.round.LTR.internal.fasta >1st.round.LTR_internal.high_M.fasta`;


%temp_hash=();
open(F,"1st.low2high.nohit.list") or die $!;
open(OUT,">1st.low2high.nohit.long.list") or die $!;
@LTR_repeat=<F>;
foreach $line(@LTR_repeat)
{
	$line=~/low(.+)\_[FR]_(\d+)_(\d+)/;
	my $temp=$1."_arrow_pilon_arrow_pilon_".$2."_".$3;
	$temp_hash{$temp}=1;
}
foreach $key(keys %temp_hash)
{
	print OUT $key,"\n";
}

close F;
close OUT;
`perl  $LTR_detector/Script/Fasta_manipulate.pl --method extract --list 1st.low2high.nohit.long.list --fasta 1st.round.LTR.internal.fasta >1st.low2high.nohit.internal.fasta`;



mkdir "$dir/BlastDB";
`$blast/makeblastdb  -in $RTase -input_type fasta -dbtype prot -out $dir/BlastDB/Arabidopsis_RTase`;
`$blast/blastx -query 1st.round.LTR.internal.fasta -db $dir/BlastDB/Arabidopsis_RTase  -evalue 1e-10 -num_threads $num_threads -out 1st.round.high.M.LTR_Arabidopsis_domain_blast.out`;


open(F,"1st.round.high.M.LTR_Arabidopsis_domain_blast.out") or die $!;
open(OUT1,">temp1.txt") or die $!;
@blastx=<F>;
$count=1;
for (my $i=0;$i<=$#blastx;$i++)
{
	if($blastx[$i]=~/(\>\s+((Copia_RT)|(Gypsy_RT)))/)
	{
		print OUT1 ">",$2,$count++,"\n";
	}
	elsif($blastx[$i]=~/Query/ && $blastx[$i]!~/arrow_pilon/)
	{
		@temp=split(" ",$blastx[$i]);
		$temp[2]=~s/-//g;
		print OUT1 $temp[2],"\n";
	}
}
@temp=();
close F;
close OUT1;

$count=0;
open(F,"temp1.txt") or die $!;
open(OUT2,">temp2.txt") or die $!;
@blastx=<F>;
for (my $i=0;$i<=$#blastx;$i++)
{
	unless ($blastx[$i]=~/^$/)
	{
		chomp $blastx[$i];
		if($blastx[$i]=~/>/)
		{
			$count++;
			if($count==1)
			{
				print OUT2 $blastx[$i],"\n"
			}
			else
			{
				print OUT2 "\n",$blastx[$i],"\n"
			}
		}
		else
		{
			print OUT2 $blastx[$i];
		}
	}
}

close F;
close OUT2;


open(F,"temp2.txt") or die $!;
open(OUT3,">1st.round.high.M.RTase.domain.fasta") or die $!;
@blastx=<F>;
for(my $i=0;$i<$#blastx;$i+=2)
{
	if($blastx[$i]=~/Copia/ && length($blastx[$i+1])>=100 && length($blastx[$i+1])<=108)
	{
		print OUT3 $blastx[$i],$blastx[$i+1]
	}
	elsif($blastx[$i]=~/Gypsy/ && length($blastx[$i+1])>=145 && length($blastx[$i+1])<=155)
	{
		print OUT3 $blastx[$i],$blastx[$i+1]
	}
} 
close F;
close OUT3;
unlink "temp1.txt";
unlink "temp2.txt";

# Finish first round DB1 RTase construction
`$blast/makeblastdb  -in 1st.round.high.M.RTase.domain.fasta -input_type fasta -dbtype prot -out $dir/BlastDB/RTase`;

#For high confident singleton and low confident LTRs, compare them with high confident multiple LTR/RTase database, pick out those with hits to either database and make 2nd LTR DB

open(F,"1st.low2high.has.hit") or die $!;
open(OUT,">1st.low2high.has.hit.renamed.list") or die $!;
@LTR_repeat=<F>;
foreach $line(@LTR_repeat)
{

	$line=~s/low//;
	print OUT $line;
}
close F;
close OUT;

`$blast/blastx -query 1st.low2high.nohit.internal.fasta -db $dir/BlastDB/RTase  -evalue 1e-10 -num_threads $num_threads -out 1st.low2high.nohit.RTase_domain_blast.out -outfmt 6`;
open(F,"1st.low2high.nohit.RTase_domain_blast.out") or die $!;
open(OUT,">1st.low2high.nohit.has.RTase_domain.list") or die $!;
@LTR_repeat=<F>;
%temp_hash=();
foreach $line(@LTR_repeat)
{
	$line=~/(.+?)\s+/;
	$temp_hash{$1}=$line;
}
foreach $key(keys %temp_hash)
{
	$key=~/(.+?)_arrow_pilon_arrow_pilon_(\d+)_(\d+)/;
	print OUT $1,"_F_",$2,"_",$3,"\n";
	print OUT $1,"_R_",$2,"_",$3,"\n";
}
close F;
close OUT;


`cat 1st.low2high.has.hit.renamed.list 1st.round.LTRepeat.high.multiple.list >1st.round.DBI.LTRepeat.list`;
`perl  $LTR_detector/Script/Remove_low2high_nested.LTR.pl 1st.round.DBI.LTRepeat.list 1st.low2high.nohit.has.RTase_domain.list`;
`cat 1st.round.LTRepeat.low2high.nested.removed.has.RTase.list 1st.round.DBI.LTRepeat.list >1st.round.LTR.DB.list`;
`perl  $LTR_detector/Script/Fasta_manipulate.pl --method extract --list 1st.round.LTR.DB.list --fasta 1st.round.LTRepeat.fasta >1st.round.LTRepeat.DB.fasta`;

open(F,"1st.round.LTR.DB.list") or die $!;
open(OUT,">1st.round.LTR.DB.long.list") or die $!;
@LTR_repeat=<F>;
%temp_hash=();
foreach $line(@LTR_repeat)
{
	$line=~/(.+)\_[FR]_(\d+)_(\d+)/;
	my $temp=$1."_arrow_pilon_arrow_pilon_".$2."_".$3;
	$temp_hash{$temp}=1;
}
foreach $key(keys %temp_hash)
{
	print OUT $key,"\n";
}

close F;
close OUT;
`perl  $LTR_detector/Script/Fasta_manipulate.pl --method extract --list 1st.round.LTR.DB.long.list --fasta 1st.round.LTR.full.fasta >1st.round.LTR.DB.full.fasta`;


if( -d "$dir/Results" )
{
    `rm -rf $dir/Results`;
    `mkdir $dir/Results`;
}
else
{
    `mkdir $dir/Results`;

}
`mv 1st.round.high.M.RTase.domain.fasta $dir/Results/`;
`cp 1st.round.LTRepeat.DB.fasta $dir/Results` ;

open(F,"1st.round.LTR.full.fasta") or die $!;
open(OUT1,">1st.round.LTR.list") or die $!;
@LTR_repeat=<F>;
foreach $line(@LTR_repeat)
{
	if($line=~/>(.+)/)
	{
		print OUT1 $1,"\n";
	}
}
close F;
close OUT1;

#This step will masked all LTRs identified in 1st round, including thoese excluded false ones

if($rounds>1)
{
	open(F,"$dir/Temp/formated.fasta") or die $!;
	open(F2,"1st.round.LTR.list") or die $!;
	open(OUT,">1st.round.ltr.X.masked.genome.fasta") or die $!;
	@genome=<F>;
	%temp_hash=();
	for(my $i=0;$i<$#genome;$i+=2)
	{
		chomp $genome[$i];
		chomp $genome[$i+1];
		$genome[$i]=~/\>(.+?)\_arrow/;
		$temp_hash{$1}=$genome[$i+1]
	}

  	@LTR_repeat=<F2>;
	foreach $line(@LTR_repeat)
	{
		unless($line=~/^$/)
		{
			$line=~/(.+)_arrow_pilon_arrow_pilon_(\d+)_(\d+)/;
			my $sub=substr($temp_hash{$1},$2-1,$3-$2+1);
			my $n= "n" x length($sub);
			substr($temp_hash{$1},$2-1,$3-$2+1)=~s/$sub/$n/;
		}
	}


	foreach $key(keys %temp_hash)
	{
		print OUT ">",$key,"_arrow_pilon_arrow_pilon","\n",$temp_hash{$key},"\n";
	} 
	%temp_hash=();
	@genome=();
	@LTR_repeat=();
	close F;
	close F2;
	close OUT;


	open(F,"1st.round.ltr.X.masked.genome.fasta") or die $!;
	open(OUT,">1st.round.ltr.X.removed.genome.fasta") or die $!;
	@genome=<F>;
	foreach $line(@genome)
	{
		if($line=~/>/)
		{
			print OUT $line;
		}
		else
		{
			$line=~s/n//g;
			print OUT $line;
		}
	}
	@genome=();
	close F;
	close OUT;
	if($method eq "LTR_finder")
	{
		`$LTR_finder/ltr_finder  -s $LTR_finder/tRNAdb/Athal-tRNAs.fa 1st.round.ltr.X.removed.genome.fasta >2nd.round.LTR.out`;
		`$LTR_detector/Script/LTRepeat_extract.pl 1st.round.ltr.X.removed.genome.fasta 2nd.round.LTR.out >2nd.round.LTRepeat.txt`;
		`$LTR_detector/Script/LTR_full_extract.pl 1st.round.ltr.X.removed.genome.fasta 2nd.round.LTR.out >2nd.round.LTR.full.fasta`;
	}

	open(F,"2nd.round.LTRepeat.txt") or die $!;
	open(OUT1,">2nd.round.LTRepeat.fasta") or die $!;
	@LTR_repeat=<F>;
	foreach $line(@LTR_repeat)
	{	
		unless($line=~/^$/)
		{
			my @F=split(" ",$line);
			if($F[0]=~/(.+)\_arrow\_pilon\_arrow\_pilon\_(\d+)_(\d+)/)
			{
				print OUT1 ">",$1,"_F_",$2,"_",$3,"\n",$F[5],"\n" if length($F[5]) >=$LTR_length && length($F[9]) >=$LTR_length;
				print OUT1 ">",$1,"_R_",$2,"_",$3,"\n",$F[9],"\n" if length($F[5]) >=$LTR_length && length($F[9]) >=$LTR_length;
			}
		}
	}
	@LTR_repeat=();
	close F;
	close OUT1;

	open(F,"2nd.round.LTRepeat.fasta") or die $!;
	open(OUT,">2nd.round.LTRepeat.renamed.fasta") or die $!;
	@LTR_repeat=<F>;
	foreach $line(@LTR_repeat)
	{
		if($line=~/>/)
		{
			$line=~s/>/>low/;
			print OUT $line;
		}
		else
		{
			print OUT $line;
		}
	}
	@LTR_repeat=();
	close F;
	close OUT1;	

	`cat $dir/Results/1st.round.LTRepeat.DB.fasta 2nd.round.LTRepeat.renamed.fasta >1st.round.DB.and.2nd.round.candidate.fasta`;
	`$cd_hit/cd-hit-est -i 1st.round.DB.and.2nd.round.candidate.fasta -o 1st.round.LTRepeat.2nd.round.candidates  -c 0.8 -n 5 -d 0 -M 80000 -T 50 -aS 0.8 -aL 0.8 -sc 1 -b 40`;
	`perl  $LTR_detector/Script/LTRepeat_low2high_cluster.process2.pl 1st.round.LTRepeat.2nd.round.candidates.clstr`;
	
	open(F,"2nd.low2high.has.hit") or die $!;
	open(OUT,">2nd.nested.list") or die $!;
	@LTR_repeat=<F>;
	foreach $line(@LTR_repeat)
	{

		$line=~s/low//;
		print OUT $line;
	}
	@LTR_repeat=();
	close F;
	close OUT;

	`perl  $LTR_detector/Script/Fasta_manipulate.pl --method extract --list 2nd.nested.list --fasta 2nd.round.LTRepeat.fasta >2nd.round.nested.fasta`;

	open(F,"2nd.round.nested.fasta") or die $!;
	open(OUT,">2nd.round.nested.renamed.fasta") or die $!;
	@LTR_repeat=<F>;
	foreach $line(@LTR_repeat)
	{
		if($line=~/>(.+)/)
		{
			print OUT ">nested1_",$1,"\n";
		}	
		else
		{
			print OUT $line;
		}
	}
	@LTR_repeat=();
	close F;
	close OUT;
}



#Round3
if($rounds>2)
{
	open(F,"2nd.round.LTR.full.fasta") or die $!;
	open(OUT1,">2nd.round.LTR.list") or die $!;
	@LTR_repeat=<F>;
	foreach $line(@LTR_repeat)
	{
		if($line=~/>(.+)/)
		{
			print OUT1 $1,"\n";
		}
	}
	close F;
	close OUT1;

	open(F,"1st.round.ltr.X.removed.genome.fasta") or die $!;
	open(F2,"2nd.round.LTR.list") or die $!;
	open(OUT,">2nd.round.ltr.X.masked.genome.fasta") or die $!;
	@genome=<F>;
	%temp_hash=();
	for(my $i=0;$i<$#genome;$i+=2)
	{
		chomp $genome[$i];
		chomp $genome[$i+1];
		$genome[$i]=~/\>(.+?)\_arrow/;
		$temp_hash{$1}=$genome[$i+1]
	}

 	@LTR_repeat=<F2>;
	foreach $line(@LTR_repeat)
	{
		unless($line=~/^$/)
		{
			$line=~/(.+)_arrow_pilon_arrow_pilon_(\d+)_(\d+)/;
			my $sub=substr($temp_hash{$1},$2-1,$3-$2+1);
			my $n= "n" x length($sub);
			substr($temp_hash{$1},$2-1,$3-$2+1)=~s/$sub/$n/;
		}
	}

	foreach $key(keys %temp_hash)
	{
		print OUT ">",$key,"_arrow_pilon_arrow_pilon","\n",$temp_hash{$key},"\n";
	} 
	%temp_hash=();
	@genome=();
	@LTR_repeat=();
	close F;
	close F2;
	close OUT;


	open(F,"2nd.round.ltr.X.masked.genome.fasta") or die $!;
	open(OUT,">2nd.round.ltr.X.removed.genome.fasta") or die $!;
	@genome=<F>;
	foreach $line(@genome)
	{
		if($line=~/>/)
		{
			print OUT $line;
		}
		else
		{
			$line=~s/n//g;
			print OUT $line;
		}
	}
	@genome=();
	close F;
	close OUT;

	if($method eq "LTR_finder")
	{
		`$LTR_finder/ltr_finder  -s $LTR_finder/tRNAdb/Athal-tRNAs.fa 2nd.round.ltr.X.removed.genome.fasta >3rd.round.LTR.out`;
		`$LTR_detector/Script/LTRepeat_extract.pl 2nd.round.ltr.X.removed.genome.fasta 3rd.round.LTR.out >3rd.round.LTRepeat.txt`;
		`$LTR_detector/Script/LTR_full_extract.pl 2nd.round.ltr.X.removed.genome.fasta 3rd.round.LTR.out >3rd.round.LTR.full.fasta`;
	}

	open(F,"3rd.round.LTRepeat.txt") or die $!;
	open(OUT1,">3rd.round.LTRepeat.fasta") or die $!;
	@LTR_repeat=<F>;
	foreach $line(@LTR_repeat)
	{	
		unless($line=~/^$/)
		{
			my @F=split(" ",$line);
			if($F[0]=~/(.+)\_arrow\_pilon\_arrow\_pilon\_(\d+)_(\d+)/)
			{
				print OUT1 ">",$1,"_F_",$2,"_",$3,"\n",$F[5],"\n" if length($F[5]) >=$LTR_length && length($F[9]) >=$LTR_length;
				print OUT1 ">",$1,"_R_",$2,"_",$3,"\n",$F[9],"\n" if length($F[5]) >=$LTR_length && length($F[9]) >=$LTR_length;
			}
		}
	}
	@LTR_repeat=();
	close F;
	close OUT1;

	open(F,"3rd.round.LTRepeat.fasta") or die $!;
	open(OUT,">3rd.round.LTRepeat.renamed.fasta") or die $!;
	@LTR_repeat=<F>;
	foreach $line(@LTR_repeat)
	{
		if($line=~/>/)
		{
			$line=~s/>/>low/;
			print OUT $line;
		}
		else
		{
			print OUT $line;
		}
	}
	@LTR_repeat=();
	close F;
	close OUT1;	

	`cat $dir/Results/1st.round.LTRepeat.DB.fasta 3rd.round.LTRepeat.renamed.fasta >1st.round.DB.and.3rd.round.candidate.fasta`;
	`$cd_hit/cd-hit-est -i 1st.round.DB.and.3rd.round.candidate.fasta -o 1st.round.LTRepeat.3rd.round.candidates  -c 0.8 -n 5 -d 0 -M 80000 -T 50 -aS 0.8 -aL 0.8 -sc 1 -b 40`;
	`perl  $LTR_detector/Script/LTRepeat_low2high_cluster.process3.pl 1st.round.LTRepeat.3rd.round.candidates.clstr`;
	
	open(F,"3rd.low2high.has.hit") or die $!;
	open(OUT,">3rd.nested.list") or die $!;
	@LTR_repeat=<F>;
	foreach $line(@LTR_repeat)
	{

		$line=~s/low//;
		print OUT $line;
	}
	@LTR_repeat=();
	close F;
	close OUT;

	`perl  $LTR_detector/Script/Fasta_manipulate.pl --method extract --list 3rd.nested.list --fasta 3rd.round.LTRepeat.fasta >3rd.round.nested.fasta`;

	open(F,"3rd.round.nested.fasta") or die $!;
	open(OUT,">3rd.round.nested.renamed.fasta") or die $!;
	@LTR_repeat=<F>;
	foreach $line(@LTR_repeat)
	{
		if($line=~/>(.+)/)
		{
			print OUT ">nested2_",$1,"\n";
		}	
		else
		{
			print OUT $line;
		}
	}
	@LTR_repeat=();
	close F;
	close OUT;
}

# add more rounds above 
# Solo LTR identification
open(F,"2nd.nested.list") or die $!;
open(OUT,">2nd.nested.long.list") or die $!;
@LTR_repeat=<F>;
%temp_hash=();
foreach $line(@LTR_repeat)
{
	$line=~/(.+)\_[FR]_(\d+)_(\d+)/;
	my $temp=$1."_arrow_pilon_arrow_pilon_".$2."_".$3;
	$temp_hash{$temp}=1;
}
foreach $key(keys %temp_hash)
{
	print OUT $key,"\n";
}
@LTR_repeat=();
%temp_hash=();
close F;
close OUT;

open(F,"3rd.nested.list") or die $!;
open(OUT,">3rd.nested.long.list") or die $!;
@LTR_repeat=<F>;
%temp_hash=();
foreach $line(@LTR_repeat)
{
	$line=~/(.+)\_[FR]_(\d+)_(\d+)/;
	my $temp=$1."_arrow_pilon_arrow_pilon_".$2."_".$3;
	$temp_hash{$temp}=1;
}
foreach $key(keys %temp_hash)
{
	print OUT $key,"\n";
}
@LTR_repeat=();
%temp_hash=();
close F;
close OUT;





`cat $dir/Results/1st.round.LTRepeat.DB.fasta *nested.renamed.fasta >All.LTRepeat.identified.fasta`;
`cp All.LTRepeat.identified.fasta $dir/Results/`;
`$blast/makeblastdb  -in All.LTRepeat.identified.fasta -input_type fasta -dbtype nucl -out $dir/BlastDB/LTRepeat`;


`perl  $LTR_detector/Script/Fasta_manipulate.pl --method extract --list 1st.round.LTR.DB.long.list --fasta 1st.round.LTR.full.fasta >1st.round.LTR.DB.full.fasta`;
`perl  $LTR_detector/Script/Fasta_manipulate.pl --method extract --list 2nd.nested.long.list --fasta 2nd.round.LTR.full.fasta >2nd.round.nested.LTR.full.fasta`;
`perl  $LTR_detector/Script/Fasta_manipulate.pl --method extract --list 3rd.nested.long.list --fasta 3rd.round.LTR.full.fasta >3rd.round.nested.LTR.full.fasta`;
`perl  $LTR_detector/Script/Fasta_manipulate.pl --method extract --list 1st.round.LTR.DB.long.list --fasta 1st.round.LTR.internal.fasta >1st.round.LTR.DB.internal.fasta`;

`cat 1st.round.LTR.DB.full.fasta 2nd.round.nested.LTR.full.fasta 3rd.round.nested.LTR.full.fasta >All.LTR.identified.full.fasta`;
`$blast/makeblastdb  -in All.LTR.identified.full.fasta -input_type fasta -dbtype nucl -out $dir/BlastDB/LTRfull`;
`$blast/makeblastdb  -in 1st.round.LTR.DB.internal.fasta -input_type fasta -dbtype nucl -out $dir/BlastDB/LTRinternal`;
`$blast/blastn -db $dir/BlastDB/LTRfull -out soloLTR.genome.blast.out -query  1st.round.ltr.X.masked.genome.fasta -num_threads $num_threads -evalue 1e-10 -outfmt 6`;


open(F,"All.LTRepeat.identified.fasta") or die $!;
open(OUT,">All.LTRepeat.identified.length") or die $!;
@LTR_repeat=<F>;
foreach $line(@LTR_repeat)
{
	if ($line=~/>(.+)/)
	{
		print OUT $1,"\t" 
	}
	else
	{
		print OUT length($line),"\n";
	}
}
@LTR_repeat=();
close F;
close OUT;


open(F,"All.LTRepeat.identified.length") or die $!;
open(F2,"soloLTR.genome.blast.out") or die $!;
open(OUT,">solo.LTR.candidate.txt") or die $!;
@LTR_repeat=<F>;
@blastn=<F2>;
foreach $line(@LTR_repeat)
{
	$line=~/(.+)\s+(.+)/;
	$temp_hash{$1}=$2;
}

foreach $line(@blastn)
{
	chomp $line;
	@temp=split(" ",$line);
	$temp[1]=~/(.+)_arrow_pilon_arrow_pilon_(\d+)_(\d+)/;
	$temp=$1."_F_".$2."_".$3;
	if(defined $temp_hash{$temp})
	{
		print OUT $line,"\t",$temp_hash{$temp},"\n" if $temp[3]>$temp_hash{$temp}*0.8 && $temp[3]<$temp_hash{$temp}*1.2 && $temp[9]<$temp_hash{$temp}*1.2;
		next;
	}
	else
	{
		$temp=$1."_R_".$2."_".$3;
		if(defined $temp_hash{$temp})
		{
			print OUT $line,"\t",$temp_hash{$temp},"\n" if $temp[3]>$temp_hash{$temp}*0.8 && $temp[3]<$temp_hash{$temp}*1.2 && $temp[9]<$temp_hash{$temp}*1.2;
		}
	}
}
@LTR_repeat=();
@blastn=();
%temp_hash=();
@temp=();
close F;
close F2;
close OUT;


open(F,"solo.LTR.candidate.txt") or die $!;
open(OUT,">solo.LTR.candidates.txt") or die $!;
@LTR_repeat=<F>;
foreach $line(@LTR_repeat)
{
	unless($line=~/^$/)
	{
		chomp $line;
		$temp_hash{$line}=1;
	}
}

foreach $key(keys %temp_hash)
{
	@temp=split(" ",$key);
	print OUT $temp[0],"\t",$temp[1],"\t",$temp[6],"\t",$temp[7],"\n";
}
%temp_hash=();
@temp=();
@LTR_repeat=();
close F;
close OUT;
unlink "solo.LTR.candidate.txt";


open(F,"solo.LTR.candidates.txt") or die $!;
open(OUT,">solo.LTRs.candidates.txt") or die $!;
@LTR_repeat=<F>;
foreach $line(@LTR_repeat)
{
	@temp=split(" ",$line);
	if(defined $temp_hash{$temp[2]})
	{
		if($temp_hash{$temp[2]} eq $temp[0])
		{
			next;
		}
		else
		{
			$temp_hash{$temp[2]}=$temp[0];
			$temp_hash{$temp[3]}=$temp[0];
			print OUT $line;
			next;
		}
	}
	elsif(defined $temp_hash{$temp[3]} )
	{
		if($temp_hash{$temp[3]} eq $temp[0])
		{
			next;
		}
		else
		{
			$temp_hash{$temp[2]}=$temp[0];
			$temp_hash{$temp[3]}=$temp[0];
			print OUT $line;
			next;

		}
	}
	else
	{
			$temp_hash{$temp[2]}=$temp[0];
			$temp_hash{$temp[3]}=$temp[0];
			print OUT $line;
			next;
	}
	
}
%temp_hash=();
@temp=();
@LTR_repeat=();
close F;
close OUT;
unlink "solo.LTR.candidates.txt";


#Check TSD
`split solo.LTRs.candidates.txt -l380`;
`cp $LTR_detector/Script/*sh .`;

open(F,"Run.TSD.checker.sh") or die $!;
open(OUT,">Run.TSDs.checker.sh") or die $!;
@LTR_repeat=<F>;
foreach $line(@LTR_repeat)
{
	if($line=~/LTR_detector/)
	{
		$line=~s/LTR_detector/$LTR_detector/;
		$line=~s/Genome_seq/$dir/;		
		print OUT $line;
	}
	else
	{
		print OUT $line;
	}
}
@LTR_repeat=();
close F;
close OUT;
unlink "Run.TSD.checker.sh";
`bash  Run.TSDs.checker.sh; bash Run.TSD.checker.batch.sh`;

$job=1;
while($job)
{
	sleep (60);
	my @job=`ps -U root -u root --deselect -F`;
	foreach $line(@job)
	{
		if($line=~/Solo_LTR_TSD_checker/)
		{
			$job=2;
		}
	}
	if($job==1)
	{
		last;
	}
}

`bash Combine.TSD.check.sh`;
`perl  $LTR_detector/Script/Solo_LTR_cleaning.pl perfect.TSD.list 1mismatch.TSD.list No_TSD_found.list`;

#find nested LTR coordinate and remove solo ltr overlap with those nested LTRs



`cat 2nd.round.nested.LTR.full.fasta 3rd.round.nested.LTR.full.fasta >nested.LTR.full.fasta`;
`rm -rf *.dir;rm Run.x*;`;
`$blast/makeblastdb  -in 1st.round.ltr.X.masked.genome.fasta -input_type fasta -dbtype nucl -out $dir/BlastDB/Genome`;
`$blast/blastn -db $dir/BlastDB/Genome -out Nested.LTR.genome.blast.out -query  nested.LTR.full.fasta -num_threads $num_threads -evalue 1e-10 -outfmt 6`;
`perl $LTR_detector/Script/Nested_LTRs_coordinate_finder.pl Nested.LTR.genome.blast.out`;


open(F,"Nested_LTR_coordinate.txt") or die $!;
open(F2,"perfect_TSD.uniq.list") or die $!;
open(OUT2,">perfect_TSD.uniq.nested.removed.list") or die $!;

@nested=<F>;
@LTR_repeat=<F2>;
foreach $line(@LTR_repeat)
{
	my $counter=0;
	$line=~/(.+)\s+(.+)\s+(\d+)\s+(\d+)/;
	my $chr=$1;
	my $start=$3;
	my $end=$4;
	foreach my $nest(@nested)
	{
		@temp=split(" ",$nest);
		$temp[4]=~/(.+)\:(.+)/;
		@temp2=split(",",$2);
		if($chr eq $temp[1])
		{
			unless($start >$temp2[1]||$end < $temp2[0])
			{
				$counter=1;
			}
			if(defined $temp2[2])
			{
				unless($start >$temp2[3]||$end < $temp2[2])
				{
					$counter=1;
				}
			}
			if(defined $temp2[4])
			{
				unless($start >$temp2[5]||$end < $temp2[4])
				{
					$counter=1;
				}
			}
			if(defined $temp2[6])
			{
				unless($start >$temp2[7]||$end < $temp2[6])
				{
					$counter=1;
				}
			}
		}
	}
	print OUT2 $line if $counter==0;
}

@nested=();
@LTR_repeat=();
@temp=();
@temp2=();
close F;
close F2;
close OUT2;


open(F,"Nested_LTR_coordinate.txt") or die $!;
open(F3,"1mismatch_TSD.uniq.list") or die $!;
open(OUT3,">1mismatch_TSD.uniq.nested.removed.list") or die $!;

@nested=<F>;
@LTR_repeat=<F3>;
foreach $line(@LTR_repeat)
{
	my $counter=0;
	$line=~/(.+)\s+(.+)\s+(\d+)\s+(\d+)/;
	my $chr=$1;
	my $start=$3;
	my $end=$4;
	foreach my $nest(@nested)
	{
		@temp=split(" ",$nest);
		$temp[4]=~/(.+)\:(.+)/;
		@temp2=split(",",$2);
		if($chr eq $temp[1])
		{
			unless($start >$temp2[1]||$end < $temp2[0])
			{
				$counter=1;
			}
			if(defined $temp2[2])
			{
				unless($start >$temp2[3]||$end < $temp2[2])
				{
					$counter=1;
				}
			}
			if(defined $temp2[4])
			{
				unless($start >$temp2[5]||$end < $temp2[4])
				{
					$counter=1;
				}
			}
			if(defined $temp2[6])
			{
				unless($start >$temp2[7]||$end < $temp2[6])
				{
					$counter=1;
				}
			}
		}
	}
	print OUT3 $line if $counter==0;
}

@nested=();
@LTR_repeat=();
@temp=();
@temp2=();
close F;
close F3;
close OUT3;

open(F,"Nested_LTR_coordinate.txt") or die $!;
open(F4,"noTSD_TSD.uniq.list") or die $!;
open(OUT4,">noTSD_TSD.uniq.nested.removed.list") or die $!;
@nested=<F>;
@LTR_repeat=<F4>;
foreach $line(@LTR_repeat)
{
	my $counter=0;
	$line=~/(.+)\s+(.+)\s+(\d+)\s+(\d+)/;
	my $chr=$1;
	my $start=$3;
	my $end=$4;
	foreach my $nest(@nested)
	{
		@temp=split(" ",$nest);
		$temp[4]=~/(.+)\:(.+)/;
		@temp2=split(",",$2);
		if($chr eq $temp[1])
		{
			unless($start >$temp2[1]||$end < $temp2[0])
			{
				$counter=1;
			}
			if(defined $temp2[2])
			{
				unless($start >$temp2[3]||$end < $temp2[2])
				{
					$counter=1;
				}
			}
			if(defined $temp2[4])
			{
				unless($start >$temp2[5]||$end < $temp2[4])
				{
					$counter=1;
				}
			}
			if(defined $temp2[6])
			{
				unless($start >$temp2[7]||$end < $temp2[6])
				{
					$counter=1;
				}
			}
		}
	}
	print OUT4 $line if $counter==0;
}

@nested=();
@LTR_repeat=();
@temp=();
@temp2=();
close F;
close F4;
close OUT4;


open(F,"1st.round.ltr.X.masked.genome.fasta") or die $!;
open(F2,"perfect_TSD.uniq.nested.removed.list") or die $!;
open(OUT,">SoloLTR_perfect_TSD.uniq.nested.removed.fasta") or die $!;
@genome=<F>;
@LTR_repeat=<F2>;

for(my $i=0;$i<$#genome;$i+=2)
{
	chomp $genome[$i];
	chomp $genome[$i+1];
	$genome[$i]=~/\>(.+)/;
	$temp_hash{$1}=$genome[$i+1];
}

foreach $line(@LTR_repeat)
{
	@temp=split(" ",$line);
	$temp=substr($temp_hash{$temp[0]},$temp[2]-1,$temp[3]-$temp[2]+2);
	print OUT ">",$temp[0],"_",$temp[2],"_",$temp[3],"\n";
	print OUT $temp,"\n";
}

@genome=();
@temp=();
@LTR_repeat=();
%temp_hash=();
close F;
close F2;
close OUT;

open(F,"1st.round.ltr.X.masked.genome.fasta") or die $!;
open(F2,"1mismatch_TSD.uniq.nested.removed.list") or die $!;
open(OUT,">SoloLTR_1mismatch_TSD.uniq.nested.removed.fasta") or die $!;
@genome=<F>;
@LTR_repeat=<F2>;

for(my $i=0;$i<$#genome;$i+=2)
{
	chomp $genome[$i];
	chomp $genome[$i+1];
	$genome[$i]=~/\>(.+)/;
	$temp_hash{$1}=$genome[$i+1];
}

foreach $line(@LTR_repeat)
{
	@temp=split(" ",$line);
	$temp=substr($temp_hash{$temp[0]},$temp[2]-1,$temp[3]-$temp[2]+2);
	print OUT ">",$temp[0],"_",$temp[2],"_",$temp[3],"\n";
	print OUT $temp,"\n";
}

@genome=();
@temp=();
@LTR_repeat=();
%temp_hash=();
close F;
close F2;
close OUT;

open(F,"1st.round.ltr.X.masked.genome.fasta") or die $!;
open(F2,"noTSD_TSD.uniq.nested.removed.list") or die $!;
open(OUT,">SoloLTR_noTSD_TSD.uniq.nested.removed.fasta") or die $!;
@genome=<F>;
@LTR_repeat=<F2>;

for(my $i=0;$i<$#genome;$i+=2)
{
	chomp $genome[$i];
	chomp $genome[$i+1];
	$genome[$i]=~/\>(.+)/;
	$temp_hash{$1}=$genome[$i+1];
}

foreach $line(@LTR_repeat)
{
	@temp=split(" ",$line);
	$temp=substr($temp_hash{$temp[0]},$temp[2]-1,$temp[3]-$temp[2]+2);
	print OUT ">",$temp[0],"_",$temp[2],"_",$temp[3],"\n";
	print OUT $temp,"\n";
}

@genome=();
@temp=();
@LTR_repeat=();
%temp_hash=();
close F;
close F2;
close OUT;

`cp SoloLTR*uniq.nested.removed.fasta $dir/Results/`;
`perl $LTR_detector/Script/LTR_cluster.80.80.80.pl  All.LTRepeat.identified.fasta`;


open(F,"Final.LTRepeat.family.txt") or die $!;
open(F1,"1st.round.LTRepeat.txt") or die $!;
open(F2,"2nd.round.LTRepeat.txt") or die $!;
open(F3,"3rd.round.LTRepeat.txt") or die $!;
open(OUT,">Final.LTRepeat.family.detailed.txt") or die $!;

my @LTR_list=<F>;
my @round1=<F1>;
my @round2=<F2>;
my @round3=<F3>;
my $temp_id;
foreach $line(@LTR_list)
{
	my $printer=0;
	@temp=split(" ",$line);
	$temp[0]=~/(.+)_[FR]_(\d+)_(\d+)/;
	$temp=$1."_arrow_pilon_arrow_pilon_".$2."_".$3;
	if(defined $temp_hash{$temp})
	{
		$temp_hash{$temp}=2;
	}
	else
	{
		$temp_hash{$temp}=1;
	}
	if($temp=~/nested1/)
	{
		foreach my $line11(@round2)
		{
			$line11=~/(.+?)\s+/;
			$temp_id=$1;
			$temp=~/nested1_(.+)/;
			if($1 eq $temp_id)
			{
				if($temp_hash{$temp}==1)
				{
					unless($printer==1)
					{
						chomp $line11;
						chomp $line;
						print OUT $temp[1],"\t",$line11,"\t","nested1","\n";
						$printer=1;
					}
				}
			}
		}
	}
	elsif($temp=~/nested2/)
	{
		foreach my $line11(@round3)
		{
			$line11=~/(.+?)\s+/;
			$temp_id=$1;
			$temp=~/nested2_(.+)/;
			if($1 eq $temp_id)
			{
				if($temp_hash{$temp}==1)
				{
					unless($printer==1)
					{
						chomp $line11;
						chomp $line;
						print OUT $temp[1],"\t",$line11,"\t","nested2","\n";
						$printer=1;
					}

				}
			}
		}
	}
	else
	{
		foreach my $line11(@round1)
		{
			$line11=~/(.+?)\s+/;
			$temp_id=$1;
			if($1 eq $temp)
			{
				if($temp_hash{$temp}==1)
				{
					unless($printer==1)
					{
						chomp $line;
						chomp $line11;
						print OUT $temp[1],"\t",$line11,"\n";
						$printer=1;
					}
				}
			}
		}
	}

}
@LTR_list=();
@round1=();
@round2=();
@round3=();
@temp=();
%temp_hash=();
close F;
close F1;
close F2;
close F3;
close OUT;
`cp Final.LTRepeat.family.detailed.txt ../Results/`;

open(F,"dictionary.list") or die $!;
@LTR_list=<F>;
my %hash_name1;
my %hash_name2;
my $file;
foreach $line(@LTR_list)
{
	@temp=split(" ",$line);
	$temp[1]=~/(.+?_)arrow_pilon_arrow_pilon/;
	$hash_name1{$temp[1]}=$temp[0];
	$hash_name2{$1}=$temp[0]."_";

}

chdir "$dir/Results";
my @result_files=glob("*");
foreach $file(@result_files)
{
	open(F2,"$file") or die $!;
	open(OUT,">final.$file") or die $!;
	@temp=<F2>;
	foreach $line(@temp)
	{
		if($line=~/(000.+_arrow_pilon_arrow_pilon)/)
		{	
			$temp=$1;
			$line=~s/$temp/$hash_name1{$temp}/;
			print OUT $line;

		}
		elsif($line=~/(000.+?_)/)
		{
			$temp=$1;
			$line=~s/$temp/$hash_name2{$temp}/;
			print OUT $line;

		}
		else
		{
			print OUT $line;
		}
	}
}


sub log_N
{
  my $num = shift;
  my $base = shift;
  return log($num)/log($base);
}
