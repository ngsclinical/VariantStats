## created by Arun Rawat Version 2.1
## Modified on Aug 8th 2012
## This version does not generate other files like mean, std dev.
## #!c:\perl\bin\perl.exe
## #!/usr/bin/perl -w

use strict;
#use warnings;
use POSIX qw(ceil floor);

my $ARGCNT = $#ARGV + 1;
if ($ARGCNT != 5)
{
        print "usage: perl <perl_script> <MINOR VARIANT FILE> <UPPER LIMIT> <INTERVAL> <OUTPUT FILE> <TOTAL ERROR CALLS (y/n)> \n";
            exit(1);
}

(open(minorFH, "<$ARGV[0]")) or die "Can't open file : $ARGV[0] \n";
my $upperRange=$ARGV[1];
my $interval=$ARGV[2];
(open(LogOF,"> $ARGV[3]_DISTRIBUTION.tsv")) or die "Cannot open file $ARGV[3]_STATISTICS";
my $agg_call=$ARGV[4];

my $header=<minorFH>;
my @inFileArr="";
my (@nucA,@nucC,@nucG,@nucT)=();
my (@cf_A,@cf_C,@cf_G,@cf_T)=();
my (@Nuc_2d,@cmlFreqArr)=();
my ($sumA,$sumC,$sumG,$sumT)=0;
my ($covA,$covC,$covG,$covT)=0;
my ($arrCntA,$arrCntC,$arrCntG,$arrCntT)=0;
my $flag="N";

while (<minorFH>)
{
	@inFileArr = split('\t',$_);
	
	if ($agg_call eq "y")
	{
		if ($inFileArr[1] ne undef)
		{
			if ($inFileArr[1]<$upperRange)
			{
			push(@nucA,$inFileArr[1]);
			$sumA=$sumA+$inFileArr[1];
			$covA=$covA+$inFileArr[6];
			$arrCntA++;
			}
		} 
	}
	elsif ($agg_call eq "n")
	{
		if ($inFileArr[1] ne undef)
		{
			if ($inFileArr[1]<$upperRange)
			{
			push(@nucA,$inFileArr[1]);
			$sumA=$sumA+$inFileArr[1];
			$covA=$covA+$inFileArr[6];
			$arrCntA++;
			}
		} 
		if ($inFileArr[2] ne undef)
		{
			if ($inFileArr[2]<$upperRange)
			{
			push(@nucC,$inFileArr[2]);
			$sumC=$sumC+$inFileArr[2];
			$covC=$covC+$inFileArr[6];
			$arrCntC++;
			}
		}
		if ($inFileArr[3] ne undef)
		{
			if ($inFileArr[3]<$upperRange)
			{
			push(@nucG,$inFileArr[3]);
			$sumG=$sumG+$inFileArr[3];
			$covG=$covG+$inFileArr[6];
			$arrCntG++;
			}
		}
		if ($inFileArr[4] ne undef)
		{
			if ($inFileArr[4]<$upperRange)
			{
			push(@nucT,$inFileArr[4]);
			$sumT=$sumT+$inFileArr[4];
			$covT=$covT+$inFileArr[6];
			$arrCntT++;
			}
		}
	}
}

&CumulativeCount;


sub CumulativeCount
	{
	
		if ($agg_call eq "y")
		{
			my $tmpMeanA=0;
			if ($arrCntA>0){$tmpMeanA= sprintf("%.3f",$sumA/$arrCntA);}
			my @sortArrA=sort{$a<=>$b}@nucA;
			&distribution(\@nucA,\@cf_A);
		}
		elsif ($agg_call eq "n")
		{
			my ($tmpMeanA,$tmpMeanC,$tmpMeanG,$tmpMeanT)=0;
			if ($arrCntA>0){$tmpMeanA= sprintf("%.3f",$sumA/$arrCntA);}
			if ($arrCntC) {$tmpMeanC= sprintf("%.3f",$sumC/$arrCntC);}
			if ($arrCntG) {$tmpMeanG= sprintf("%.3f",$sumG/$arrCntG); }
			if ($arrCntT) {$tmpMeanT= sprintf("%.3f",$sumT/$arrCntT);}
		
			my @sortArrA=sort{$a<=>$b}@nucA;
			my @sortArrC=sort{$a<=>$b}@nucC;
			my @sortArrG=sort{$a<=>$b}@nucG;
			my @sortArrT=sort{$a<=>$b}@nucT;
		
			&distribution(\@nucA,\@cf_A);
			&distribution(\@nucC,\@cf_C);
			&distribution(\@nucG,\@cf_G);
			&distribution(\@nucT,\@cf_T);
		}
	}

sub distribution
 {
   my ($tmp_nucArr,$cf_Nuc)=@_;
   my @nucArr=@{$tmp_nucArr};

   my %cFreq;
   $cFreq{ceil($_ /$interval)}++ for @nucArr;
   my $max=0;
   
   $max=$upperRange/$interval;
   my $j=0;
   #print "lowerRange=" . $lowerRange. "\t". $max."\n";
   for (my $i = 1; $i <= $max; $i++)
   {

     my $cmlRange       =  ($i) * $interval;
     my $frequency = $cFreq{$i} || 0;

    
     if ($flag eq "N")
     {
     	@cmlFreqArr[$i-1]=$cmlRange;
     }
     $cf_Nuc->[$i-1]=$frequency;
   }
   
	$flag="Y";

 }

my @transposeArr=(); 
if ($agg_call eq "y")
{
	@transposeArr=(\@cf_A);
	
}
elsif ($agg_call eq "n")
{	
	@transposeArr=(\@cf_A,\@cf_C,\@cf_G,\@cf_T);
}
	for my $row(@transposeArr) 
	{
		for my $col(0..$#{$row})
		{
		push(@{$Nuc_2d[$col]}, $row->[$col]);
		}
	}

my $i=0;

if ($agg_call eq "y")
	{
	print LogOF "Interval" ."\t". "AGGREGATE_ERROR" ."\n";
	}
elsif ($agg_call eq "n")
	{
	print LogOF "Interval" . "\t". "A" . "\t". "C"."\t"."G"."\t". "T"."\n";
	}

foreach my $row(@Nuc_2d) {
   #print LogOF $nuc[$i] . "\t";
   print LogOF $cmlFreqArr[$i] . "\t";
   foreach my $val(@$row) {
      print LogOF $val ."\t";
   }
   print LogOF"\n";
   $i++;
}