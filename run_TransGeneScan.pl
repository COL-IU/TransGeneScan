#!/usr/bin/perl -w
use strict;
use Getopt::Long;

my $genome_file = "";
my $FGS_result = "";
my $FGS_whole = "1";
my $FGS_train_file = "trans";
my $command;
my $program = $0;
my $dir = substr($0, 0, length($0)-20);
my $train_file;
my $threadnum = 1;
my $starttime = time();
my $endtime;

GetOptions(
           'in=s' => \$genome_file,
           'out=s' => \$FGS_result,
           'thread=s' => \$threadnum,
           );

if (length($genome_file)==0){
    print "ERROR: An input genome file was not specified.\n";
    print_usage();
    exit;
}elsif (! -e $genome_file){
    print "ERROR: The input genome file [$genome_file] does not exist.\n";
    print_usage();
    exit;
}

if (length($FGS_result) == 0 ){
    print "ERROR: An output file name was not specified.\n";
    print_usage();
    exit;
}

if ($threadnum < 1)
{
   print "ERROR: Thread number [$threadnum] error,\n";
   print_usage();
   exit;
}

$command = $dir."TransGeneScan";
$command .= " -s ".$genome_file;
$command .= " -o ".$FGS_result;
$command .= " -w ".$FGS_whole ;
$command .= " -t ".$FGS_train_file;
$command .= " -p ".$threadnum;
print "$command\n";
system($command); 

if ($FGS_whole eq "1"){
    system($dir."post_process.pl -genome=".$genome_file." -pre=".$FGS_result." -post=".$FGS_result.".out");
    system("rm ".$FGS_result);
}else{
    system("mv ".$FGS_result." ".$FGS_result.".out");
}

system($dir."FGS_gff.py ".$FGS_result.".out ".$FGS_result.".gff");
system($dir."processFragOut.py ".$FGS_result.".out ".$genome_file." ".$FGS_result);
$endtime = time();
getElapsedTime($endtime - $starttime);

sub print_usage{

    print "USAGE: ./run_TransGeneScan.pl -in=[seq_file_name] -out=[output_file_name] (-thread=[number of thread; default 1])\n";
    print "       [seq_file_name]:    sequence file name including the full path\n";
    print "       [output_file_name]: output file name including the full path\n";
    print "       [num_thread]:       number of thread used in TransGeneScan. Default 1.\n";
}

sub getElapsedTime
{
	my $input = $_[0];
	my $hour;
	my $min;
	my $sec;
	my $str;

	$sec = $input % 60;
	$input = int($input / 60);
	$min = $input % 60;
	$input = int($input / 60);
	$hour = $input;

	$str = $hour . " hours " . $min . " minutes and " . $sec . " seconds.\n";
	print $str;
}

