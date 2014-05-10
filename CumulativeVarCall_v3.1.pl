## Run the executable as described in Usage from the Directory. 
## Created By Arun Rawat Version 1.3.1
## Modified By Arun Rawat Dated:08/03/12
#!/usr/bin/perl -w
### perl CumulativeVarCall.pl Ames1Percent.txt RV500wtplasmid_fixed.fas 1 50 100 Plasmid y 

use strict;
no warnings;
use POSIX qw(ceil floor);

my $ARGCNT = $#ARGV + 1;
if ($ARGCNT != 7)
{
        print "usage: perl <perl_script> <GATK Output File> <Reference Genome> <Minor Variant % Lower Cutoff> <Minor Variant % Upper Cutoff> <Minimum Coverage> <Project Prefix> <Generate Log File (y/n)>\n";
            exit(1);
}


(open(gatkFH, "<$ARGV[0]")) or die "Can't open file : $ARGV[0] \n";
(open(refFH, "<$ARGV[1]")) or die "Can't open file : $ARGV[1] \n";
my $minVarPercent=$ARGV[2];
my $maxVarPercent=$ARGV[3];
my $minCov=$ARGV[4];
(open(OF,">$ARGV[5]_ALL")) or die "Cannot open file $ARGV[5]_ALL";
(open(LogOF,"> $ARGV[5]_SUMMARY")) or die "Cannot open file $ARGV[5]_SUMMARY";

my $generateLog=$ARGV[6];
if ($generateLog eq "y")
{
(open(CumulOF,"> $ARGV[5]_CUMULATIVE_STATS")) or die "Cannot open file $ARGV[5]_CUMULATIVE_STATS";
(open(ZOF,">$ARGV[5]_ZERO_COVERAGE")) or die "Cannot open file $ARGV[5]_ZERO_COVERAGE";
(open(LOF,">$ARGV[5]_LOW_COVERAGE_CUTOFF")) or die "Cannot open file $ARGV[5]_LOW_COVERAGE_CUTOFF";
(open(MVOF,">$ARGV[5]_MINOR_VARIANT")) or die "Cannot open file $ARGV[5]_MINOR_VARIANT";
(open(MVOF_AVGERR,">$ARGV[5]_MINOR_VARIANT_TOTALERR")) or die "Cannot open file $ARGV[5]_MINOR_VARIANT_TOTALERR";
(open(MVOF_GATK,">$ARGV[5]_MINOR_VARIANT_GATK")) or die "Cannot open file $ARGV[5]_MINOR_VARIANT_GATK";
(open(NMVOF,">$ARGV[5]_NON_MINOR_VARIANT")) or die "Cannot open file $ARGV[5]_NON_MINOR_VARIANT";
(open(AMBOF,">$ARGV[5]_AMBIGUITY_SITES")) or die "Cannot open file $ARGV[5]_AMBIGUITY_SITES";
(open(EXOF,">$ARGV[5]_MINOR_MATCHING_REF")) or die "Cannot open file $ARGV[5]_MINOR_MATCHING_REF";
(open(MNMROF,">$ARGV[5]_MAJOR_NONMATCHING_REF")) or die "Cannot open file $ARGV[5]_MAJOR_NONMATCHING_REF";
#(open(test,">>$ARGV[5]_TEST")) or die "Cannot open file $ARGV[5]_TEST";
}

my $fastaHeader=<refFH>;
my $refSeq="";
my @rseqArr="";

#local $/="\n";
while (my $tmprefSeq=<refFH>)
{
$tmprefSeq=~ s/\r?\n$//;
chomp($tmprefSeq);
$refSeq=$refSeq . $tmprefSeq;
#$refSeq=$refSeq . $_;
}
#print length($refSeq) . "\n";
@rseqArr=split(//,$refSeq);

	my @arr=();
	my $arrSize=0;
	my $tmparrSize=0;
    my $posCnt=0;
    my $header=<gatkFH>;
	my ($cntZ,$cntL,$cntX,$cntM,$cntN,$cntW,$cntMNMR)=0;
	my ($sumA,$sumT,$sumG,$sumC)=0;
	my ($nucA,$nucC,$nucG,$nucT)=0;
	my ($nucA_Val,$nucC_Val,$nucG_Val,$nucT_Val)=0;
	my (@arrA,@arrC,@arrG,@arrT)=();
	my ($arrCntA,$arrCntC,$arrCntG,$arrCntT)=0;
	my ($avgCovA,$avgCovC,$avgCovG,$avgCovT)=0;
	my $gblCnt=0;
	my $gblFlag;
	my @nucCount=();
	
	print OF "Position" . "\t" . "A" . "\t" . "C" . "\t" . "G" . "\t" . "T" ."\t" . "Reference" . "\t" . "Coverage" . "\t" ."Flag" . "\n" ; 
	if ($generateLog eq "y")
	{
	print LOF "Position" . "\t" . "A" . "\t" . "C" . "\t" . "G" . "\t" . "T" ."\t" . "Reference" . "\t" . "Coverage" . "\t" ."Flag" . "\n" ;
	print MVOF "Position" . "\t" . "A" . "\t" . "C" . "\t" . "G" . "\t" . "T" ."\t" . "Reference" . "\t" . "Coverage" . "\t" ."Flag" . "\n" ;
	print MVOF_GATK "Position" . "\t" . "A" . "\t" . "C" . "\t" . "G" . "\t" . "T" ."\t" . "Reference" . "\t" . "Coverage" . "\t" ."Flag" . "\n" ;
	print MVOF_AVGERR "Position"."\t"."Aggregate_Error"."\t". "Reference" . "\t" . "Coverage" . "\t" ."Flag"."\n";
	print ZOF "Position" . "\t" . "A" . "\t" . "C" . "\t" . "G" . "\t" . "T" ."\t" . "Reference" . "\t" . "Coverage" . "\t" ."Flag" . "\n" ;
	print NMVOF "Position" . "\t" . "A" . "\t" . "C" . "\t" . "G" . "\t" . "T" ."\t" . "Reference" . "\t" . "Coverage" . "\t" ."Flag" . "\n" ;
	print AMBOF "Position" . "\t" . "A" . "\t" . "C" . "\t" . "G" . "\t" . "T" ."\t" . "Reference" . "\t" . "Coverage" . "\t" ."Flag" . "\n" ;
	print EXOF "Position" . "\t" . "A" . "\t" . "C" . "\t" . "G" . "\t" . "T" ."\t" . "Reference" . "\t" . "Coverage" . "\t" ."Flag" . "\n" ;
	print MNMROF "Position" . "\t" . "A" . "\t" . "C" . "\t" . "G" . "\t" . "T" ."\t" . "Reference" . "\t" . "Coverage" . "\t" ."Flag" . "\n" ;
	}
while (<gatkFH>)
{
	@arr = split('\t',$_);
	@nucCount=split(/\:|\s/,$arr[4]);
	
	
	if ($arr[1]>0)
	{
		$nucA=sprintf("%.3f",$nucCount[1]*100/$arr[1]);
	    $nucC=sprintf("%.3f",$nucCount[3]*100/$arr[1]);
	    $nucG=sprintf("%.3f",$nucCount[5]*100/$arr[1]);
	    $nucT=sprintf("%.3f",$nucCount[7]*100/$arr[1]);

		if ($arr[1]>$minCov)
		{

			if ($rseqArr[$posCnt] =~ m/(A|T|G|C|a|t|g|c)/)
			{
				&AssignNucVal;
			}
			else
			{
				print OF $posCnt+1 . "\t"."-3"."\t". "-3"."\t"."-3". "\t"."-3". "\t".$rseqArr[$posCnt] ."\t". $arr[1]."\t". "X"."\n";
				if ($generateLog eq "y")
				{
				print AMBOF $posCnt+1 ."\t".$nucA ."\t". $nucC ."\t" . $nucG ."\t" .$nucT."\t".$rseqArr[$posCnt] ."\t". $arr[1]."\t". "X"."\n";
				}
				$cntX++;
			}
		}
		else
		{
		print OF $posCnt+1 . "\t"."-2"."\t"."-2"."\t"."-2"."\t"."-2". "\t".$rseqArr[$posCnt] ."\t". $arr[1]. "\t"."L"."\n";
		if ($generateLog eq "y")
				{
				print LOF $posCnt+1 . "\t".$nucA ."\t". $nucC ."\t" . $nucG ."\t" .$nucT."\t".$rseqArr[$posCnt] ."\t". $arr[1]. "\t"."L"."\n";
				}
		$cntL++;
		}
	}
	else
	{
		print OF $posCnt+1 . "\t"."-1". "\t"."-1"."\t"."-1"."\t"."-1". "\t".$rseqArr[$posCnt] ."\t". $arr[1]."\t"."Z"."\n";
		if ($generateLog eq "y")
				{
				print ZOF $posCnt+1 . "\t"."0" ."\t". "0" ."\t" . "0" ."\t" ."0"."\t".$rseqArr[$posCnt] ."\t". $arr[1]."\t".$rseqArr[$posCnt] ."\t". $arr[1]."\t"."Z"."\n";
				}
		$cntZ++;
	}
$posCnt++;
}
#my ($meanA,$meanC,$meanG,$meanT)=0;
if ($generateLog eq "y")
{
	&CumulativeCount;
}

#if ($generateLog eq "y")
#{
#&distribution(.1,@arrA);
#}

print "Please see the log file generated for statistics" . "\n";
print LogOF "(W) Total number of positions falling within range of minor variant cutoff:" . $cntW . "\n";
print LogOF "(N) Total number of positions falling out of range of minor variant cutoff :" . $cntN . "\n";
print LogOF "(Z) Total number of positions with 0 coverage:".$cntZ . "\n";
print LogOF "(L) Total number of positions with coverage less than ". $minCov . " :". $cntL . "\n";
print LogOF "(X) Total number of positions with ambiguity call (Ambiguity call less than user defined coverage is ignored):".$cntX . "\n";
print LogOF "(M) Total number of positions with minor variant matching Reference (Exception):".$cntM . "\n";
print LogOF "(R) Total number of positions with major variant non-matching Reference:" . $cntMNMR++ . "\n";	
###Now close the files 
close LogOF; 
close OF;
if ($generateLog eq "y")
{
close CumulOF;
close ZOF;
close LOF;
close MVOF;
close NMVOF;
close AMBOF;
close EXOF;
close MNMROF;
}
#@sortArrA=();@sortArrC=();@sortArrG=();@sortArrT=();

	sub CumulativeCount
	{
	#print "calculate cumulative count=".$arrA[0];	
	#print "arrCntA=$arrCntA, \t, arrCntC=$arrCntC,\t,arrCntG=$arrCntG,\t,arrCntT=$arrCntT,\n";
	if ($arrCntA==0 || $arrCntC==0 || $arrCntG==0 || $arrCntT==0) {print "No Alignment Statistics found, please correct the reads or reference files  \n"; exit(1);}
	my $tmpMeanA= sprintf("%.3f",$sumA/$arrCntA);#&Mean(@arrA); 
	my $tmpMeanC= sprintf("%.3f",$sumC/$arrCntC);#&Mean(@arrC);
	my $tmpMeanG= sprintf("%.3f",$sumG/$arrCntG); #&Mean(@arrG);
	my $tmpMeanT= sprintf("%.3f",$sumT/$arrCntT);#&Mean(@arrT);
	print $gblCnt;
	#print $sumA ."\t". $sumC . "\t" . $sumG . "\t" . $sumT . "\t". $gblCnt . "\n";
	my $gblMean= sprintf("%.2f",($sumA+$sumC+$sumG+$sumT)/$gblCnt);
	#print $sumA ."\t". $sumC . "\t" . $sumG . "\t" . $sumT . "\t". $gblCnt. "\t" . $gblMean . "\n";
		
	my @sortArrA=sort{$a<=>$b}@arrA;
	
	my @sortArrC=sort{$a<=>$b}@arrC;
	
	my @sortArrG=sort{$a<=>$b}@arrG;
	
	my @sortArrT=sort{$a<=>$b}@arrT;
	
	print CumulOF "Nucleotide"."\t"."First Quartile" . "\t"."Median" . "\t". "Third Quartile". "\t". "Mean". "\t". "Std. Dev" . "\t". "Variables" . "\t". "Avg. Coverage"."\n";
	
	print CumulOF "A"."\t".&FirstQuartile(\@sortArrA)."\t".&Median(\@sortArrA) ."\t".&ThirdQuartile(\@sortArrA). "\t". $tmpMeanA. "\t". &StdDev($tmpMeanA,\@sortArrA) ."\t".scalar(@sortArrA) ."\t".sprintf("%.2f",$avgCovA/$arrCntA) ."\n";
	
	print CumulOF "C"."\t".&FirstQuartile(\@sortArrC)."\t".&Median(\@sortArrC). "\t".&ThirdQuartile(\@sortArrC)."\t". $tmpMeanC. "\t". &StdDev($tmpMeanC,\@sortArrC) ."\t". scalar(@sortArrC) . "\t".sprintf("%.2f",$avgCovC/$arrCntC)."\n";
	
	print CumulOF "G"."\t".&FirstQuartile(\@sortArrG)."\t".&Median(\@sortArrG). "\t".&ThirdQuartile(\@sortArrG)."\t". $tmpMeanG. "\t". &StdDev($tmpMeanG,\@sortArrG) ."\t".scalar(@sortArrG). "\t".sprintf("%.2f",$avgCovG/$arrCntG)."\n";
	
	print CumulOF "T"."\t".&FirstQuartile(\@sortArrT)."\t".&Median(\@sortArrT). "\t".&ThirdQuartile(\@sortArrT)."\t". $tmpMeanT. "\t". &StdDev($tmpMeanT,\@sortArrT) ."\t".scalar(@sortArrT). "\t".sprintf("%.2f",$avgCovT/$arrCntT)."\n";
	
	#print "$gblCnt";
	my $avgMean=($tmpMeanA+$tmpMeanC+$tmpMeanG+$tmpMeanT)/4;
	#my $avgMean=($meanA+$meanC+$meanG+$meanT)/4;
	print CumulOF "Average Mean" ."\t". $avgMean ."\n";	
	print CumulOF "Global Mean" . "\t" . $gblMean ."\n";
	#@sortArrA=();@sortArrC=();@sortArrG=();@sortArrT=();
	
	}
	
	sub AssignNucVal
	{

	my $minorFlag="F";
	my $windowFlag="F";
	my $majorFlag="F";
	$gblFlag="F";
	($nucA_Val,$nucC_Val,$nucG_Val,$nucT_Val)=0;
		
		if (($nucA >= $minVarPercent) && ($nucA < $maxVarPercent))
		{
		#this is a minor variant
		#print 
			if ($rseqArr[$posCnt] eq $nucCount[0])
			{
			#this is major variant falling in minor variant window
			$minorFlag="T";
			}
			else
			{
			$windowFlag="T";
			@arrA[$arrCntA]=$nucA;
			$avgCovA=$avgCovA+$arr[1];
			$sumA=$nucA+$sumA;
			$arrCntA++;
			
			
			if ($gblFlag eq "F") { &validateGblCount;}
			}
			
		}
		else 
		{
			if ($rseqArr[$posCnt] eq $nucCount[0])
			{
			$nucA_Val=$nucA;
			$nucA=undef;
			}
			else
			{
				if ($nucA >= $maxVarPercent)
				{
				$majorFlag="T";
				}
				else 
				{
				#@arrA[$arrCntA]=$nucA;
				#$avgCovA=$avgCovA+$arr[1];
				#$sumA=$nucA+$sumA;
				#$arrCntA++;
				}
			}
			
		}
		
		if ($nucC >= $minVarPercent && $nucC < $maxVarPercent)
		{
			if ($rseqArr[$posCnt] eq $nucCount[2])
			{
			$minorFlag="T";
			}
			else
			{
			$windowFlag="T";
			@arrC[$arrCntC]=$nucC;
			$avgCovC=$avgCovC+$arr[1];
			$sumC=$nucC+$sumC;
			$arrCntC++;
						
			if ($gblFlag eq "F") { &validateGblCount;}
			
			}

		}
		else 
		{
			if ($rseqArr[$posCnt] eq $nucCount[2])
			{
			$nucC_Val=$nucC;
			$nucC=undef;
			}
			else
			{
				if ($nucC >= $maxVarPercent)
				{
				$majorFlag="T";
				}
				else
				{
				#@arrC[$arrCntC]=$nucC;
				#$avgCovC=$avgCovC+$arr[1];
				#$sumC=$nucC+$sumC;
				#$arrCntC++;
				}
			}
		}
		
		if ($nucG >= $minVarPercent && $nucG < $maxVarPercent)
		{
			if ($rseqArr[$posCnt] eq $nucCount[4])
			{
			$minorFlag="T";
			}
			else
			{
			$windowFlag="T";
			@arrG[$arrCntG]=$nucG;
			$avgCovG=$avgCovG+$arr[1];
			$sumG=$nucG+$sumG;
			$arrCntG++;
			
			if ($gblFlag eq "F") { &validateGblCount;}
			
			}
			
		}	
		else 
		{		
			if ($rseqArr[$posCnt] eq $nucCount[4])
			{
			$nucG_Val=$nucG;
			$nucG=undef;
			}
			else
			{
				if ($nucG >= $maxVarPercent)
				{
				$majorFlag="T";
				}
				else
				{
				#@arrG[$arrCntG]=$nucG;
				#$avgCovG=$avgCovG+$arr[1];
				#$sumG=$nucG+$sumG;
				#$arrCntG++;
				}
			}
		}
		
		if ($nucT >= $minVarPercent && $nucT < $maxVarPercent)
		{
			if ($rseqArr[$posCnt] eq $nucCount[6])
			{
			$minorFlag="T";
			}
			else
			{
			$windowFlag="T";
			@arrT[$arrCntT]=$nucT;
			$avgCovT=$avgCovT+$arr[1];
			$sumT=$nucT+$sumT;
			$arrCntT++;
			
			if ($gblFlag eq "F") { &validateGblCount;}
			
			}

		}
		else 
		{
			if ($rseqArr[$posCnt] eq $nucCount[6])
			{
			$nucT_Val=$nucT;
			$nucT=undef;
			}
			else
			{
			
				if ($nucT >= $maxVarPercent)
				{
				$majorFlag="T";
				}
				else
				{
				#@arrT[$arrCntT]=$nucT;
				#$avgCovT=$avgCovT+$arr[1];
				#$sumT=$nucT+$sumT;
				#$arrCntT++;
				}
			}
		}
	
	#print "\n";
		
	if ($minorFlag eq "T")
	{
	print OF $posCnt+1 . "\t". $nucA ."\t". $nucC ."\t" . $nucG ."\t". $nucT ."\t".$rseqArr[$posCnt] ."\t". $arr[1]."\t"."M"."\n";
	if ($generateLog eq "y")
				{
				print EXOF $posCnt+1 . "\t". $nucA ."\t". $nucC ."\t" . $nucG ."\t". $nucT ."\t".$rseqArr[$posCnt] ."\t". $arr[1]."\t"."M"."\n";
				}
	$cntM++;
	}
	elsif ($majorFlag eq "T")
	{
	print OF $posCnt+1 . "\t". $nucA ."\t". $nucC ."\t" . $nucG ."\t". $nucT ."\t".$rseqArr[$posCnt] ."\t". $arr[1]."\t"."M"."\n";
	if ($generateLog eq "y")
				{
				print MNMROF $posCnt+1 . "\t". $nucA ."\t". $nucC ."\t" . $nucG ."\t". $nucT ."\t".$rseqArr[$posCnt] ."\t". $arr[1]."\t"."M"."\n";
				}
	$cntMNMR++;
	}
	else
	{
		if ($windowFlag eq "T")
		{
		print OF $posCnt+1 . "\t".$nucA ."\t". $nucC ."\t" . $nucG ."\t" .$nucT."\t".$rseqArr[$posCnt] ."\t". $arr[1]."\t"."W"."\n";
		if ($generateLog eq "y")
				{
				print MVOF $posCnt+1 . "\t".$nucA ."\t". $nucC ."\t" . $nucG ."\t" .$nucT."\t".$rseqArr[$posCnt] ."\t". $arr[1]."\t"."W"."\n";
				print MVOF_GATK $posCnt+1 . "\t".($nucA_Val>$minVarPercent ?$nucA_Val:$nucA ) ."\t".($nucC_Val>$minVarPercent ?$nucC_Val:$nucC ) ."\t" . ($nucG_Val>$minVarPercent ?$nucG_Val:$nucG) ."\t" .($nucT_Val>$minVarPercent ?$nucT_Val:$nucT )."\t".$rseqArr[$posCnt] ."\t". $arr[1]."\t"."W"."\n";
				print MVOF_AVGERR $posCnt+1 . "\t". ($nucA + $nucC+ $nucG+$nucT)."\t".$rseqArr[$posCnt] ."\t". $arr[1]."\t"."W"."\n";
				}
		
		$cntW++;
		}
		else
		{
		print OF $posCnt+1 . "\t".$nucA ."\t". $nucC ."\t" . $nucG ."\t" .$nucT."\t".$rseqArr[$posCnt] ."\t". $arr[1]."\t"."N"."\n";
		if ($generateLog eq "y")
				{
				print NMVOF $posCnt+1 . "\t".$nucA ."\t". $nucC ."\t" . $nucG ."\t" .$nucT."\t".$rseqArr[$posCnt] ."\t". $arr[1]."\t"."N"."\n";
				}
		$cntN++;
		}
	}

	}
	sub validateGblCount
	{
	$gblCnt++;
	$gblFlag="T";
	#print "inside gblFunc"."\n";
	}
	
	sub Median 
	{
    my ($myTmpArr)=@_;
    my $medianVal;
	
	my @myArr=@{$myTmpArr};

    if( (@myArr % 2) == 1 ) 
    {
        $medianVal = $myArr[((@myArr+1) / 2)-1];
    } 
    else 
    {
        $medianVal = ($myArr[(@myArr / 2)-1] + $myArr[@myArr / 2]) / 2;
    }

    return $medianVal;
	}
	
	sub FirstQuartile 
	{
    my ($myTmpArr)=@_;
    my $q1Val;
	my @myArr=@{$myTmpArr};


    if( (@myArr % 2) == 1 ) 
    {
        $q1Val = $myArr[((@myArr+1) / 4)-1];
    } 
    else 
    {
        $q1Val = ($myArr[(@myArr / 4)-1] + $myArr[@myArr / 4]) / 2;
    }

    return $q1Val;
	}
	
	sub ThirdQuartile 
	{
    my ($myTmpArr)=@_;
    my $q3Val;
	my @myArr=@{$myTmpArr};

    if( (@myArr % 2) == 1 ) 
    {
        $q3Val = $myArr[(3*(@myArr+1) / 4)-1];
    } 
    else 
    {
        $q3Val = ($myArr[(3*(@myArr / 4))-1] + $myArr[(3*(@myArr / 4))]) / 2;
    }

    return $q3Val;
	}	
	
	
	sub Mean
	{
	my (@myArr)=@_;

    my $meanVal=0;
	my ($total,$arrCnt)=0;
		foreach my $val(@myArr)
		{
		$total += $val;
	
		}
	$arrCnt=scalar @myArr;
	if ($total>0 && $arrCnt)
	{
		$meanVal=sprintf("%.3f",$total/$arrCnt);
	}
	#$$tmpMean=$meanVal;
	return $meanVal;
	}
	
	sub StdDev
	{
	my ($tmpMean,$myTmpArr)=@_;
	my @myArr=@{$myTmpArr};
	my $total=0;
	#my $meanVal;
	my ($meanVal,$stdMean,$stdDev,$arrCnt)=0;
	
	#$meanVal=&Mean(@myArr);
	$meanVal=$tmpMean;
	foreach my $val(@myArr)
	{
		$total += (($meanVal-$val)*($meanVal-$val));
	}
	$arrCnt=scalar @myArr;
	if ($arrCnt>0 && $total>0)
	{
		$stdMean=$total/$arrCnt;
	}
	$stdDev=sprintf("%.3f",sqrt($stdMean));
	return $stdDev;
	}


	

