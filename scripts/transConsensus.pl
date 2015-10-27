#!/usr/bin/perl

my $snpFile = $ARGV[0];
my $consensusFile = $ARGV[1];

open (SNP, "$snpFile") || die "cant load file";
open (DEST, ">$consensusFile") || die "cant load file";

#my @tracking = ("-","-","-");
$| = 1;
while (<SNP>)
{
	chomp;
	@input = split(" ", $_);	
	$scaffold = $input[0];
	$position = $input[1];
	$ref = $input[2];
	$N = 0;
	print DEST $scaffold;
	print DEST " ";
	print DEST $position;
	print DEST " ";
	print DEST $ref;
	$Acons = 0;
	$Ccons = 0;
	$Gcons = 0;
	$Tcons = 0;
	$allTotal = 0;
	$totalLength = scalar(@input);
	
	$x = 4;
	$MAtotal = 0;
	$A=$input[$x];
	$a=$input[$x+1];
    $C=$input[$x+2];
    $c=$input[$x+3];
    $G=$input[$x+4];
    $g=$input[$x+5];
    $T=$input[$x+6];
    $t=$input[$x+7];
		
	$Atotal = $A+$a;
    $Ctotal = $C+$c;
    $Gtotal = $G+$g;
    $Ttotal = $T+$t;
	$curNuc = "";

	$allTotal = $Atotal+$Ctotal+$Gtotal+$Ttotal;
	
	if ($allTotal eq 0)
	{
#		if (($tracking[$position % 3] eq "-") || ($tracking[$position % 3] eq "N"))
#		{
#			$curNuc = "-";
#		}else{
#			$curNuc = "N";
#		}
		$curNuc = "-";
	}else{
		$curNuc = $ref;
	}
	print DEST " ";
	print DEST $curNuc;
#	$tracking[$position % 3] = $curNuc;
	print DEST "\n";
}
