#!/usr/bin/perl

#This script aligns paired end NGS data to the illumina TruSeq or Nextera Exome or sureselect_EZ_Exome_v4 + UTR
#written by richard bagnall (r.bagnall@centenary.org.au)
#Usage: Control_exome_alignment.v1.pl -fastq [path/2/user] -exome [truseq|seqcapv4|vcrome2.1]
#save raw reads as [name1].1.fastq.gz, [name1].2.fastq.gz , [name2].1.fastq.gz , [name2].2.fastq.gz, etc...
#save the raw reads in a directory called Rawdata e.g. /home/shared/NGS/human/[USER]/Rawdata
#This script removes a space from the fastq header, which is added by NCBI

use strict; use warnings;
use File::Copy;
use Getopt::Long;
use Parallel::ForkManager;

####################################################
# Setting up command line options and usage errors #
####################################################

# set a command line option variable (-fastq); is required and is the path to the rawdata
my $input ='';
# set a command line option variable (-exome); is required and is the exome enrichment file: truseq vcrome2.1 sureselectv2 seqcapv4
my $exome ='';

GetOptions 	("fastq=s" => \$input,
            "exome=s" => \$exome);

my $offset = length($input); # use this extensively in subroutines for getting names of samples, files and folders

# -fastq is required, else print usage and die
if ($input eq '') {
	print "\n\n\t*** ERROR: You need to define the path to the raw fastq files with the -fastq option\n";
	usage();
    exit;
}
# pre-empt common error
elsif ($input =~ m/\/$/) {
	print "\n\n\t*** ERROR: Please remove the last / from the -fastq option\n";
	usage();
    exit;
}
# -exome is required, else print usage and die
if ($exome !~ m/^(truseq|seqcapv4|vcrome2.1|sureselectv2)$/) {
	print "\n\n\t*** ERROR: You need to define the exome enrichment kit used with the -exome option\n";
	usage();
    exit;
}

print "\n\n\n\t\t---------CONTROL EXOME SEQUENCING ALIGNMENT-----------\n";
print "\t\tAligning Paired-end sequence reads (Illumina 1.8+)\n";
print "\t\tUsing the $exome exome enrichment kit\n";
print "\t\tUsing BWA MEM for read alignment\n";
print "\t\tUsing Novosort to remove duplicate reads\n";
print "\t\t-------------------------------------------------------------\n\n";
print "\n";

# create temporary folder
mkdir ("$input/tempBAMs") or die "Unable to create tempBAMs directory: <$!>\n";

# paths
my $path2rawdata = "$input/Rawdata";
my $path2bam = "$input/tempBAMs";
my $path2gatk = '/home/shared/NGS/human/richardb/Applications/GenomeAnalysisTK-3.3-0/GenomeAnalysisTK.jar';
my $path2targetregions = "/home/shared/NGS/human/richardb/TargetRegions/$exome.targetregions.bed";
my $path2ref = '/home/groups/cardio/References/Bwa_b37/hs37d5.fa';

#####################
#   Start BWA mem   #
#####################

my @fq = glob("$path2rawdata/*.fastq.gz"); # make array of fastq.gz

# check that there are the same number of forward and reverse files, and they have correct extensions
my @f_fq = grep(/1.fastq.gz$/i, @fq);
my @r_fq = grep(/2.fastq.gz$/i, @fq);
if (scalar(@f_fq) != scalar(@r_fq)) {
    rmdir "$path2bam";
    die "*** error: There must be a 1.fastq.gz and 2.fastq.gz file for each sample\n\n";
}

my $sixfork_manager = Parallel::ForkManager->new(3);

for (my $i = 0; $i < @fq; $i = $i ++) {
    my @fq_pair = splice(@fq, $i, 2);
    $sixfork_manager->start and next;
    bwa_mem(join(" ", @fq_pair));
    $sixfork_manager->finish;
}

$sixfork_manager-> wait_all_children;

print "\n\n*** BWA mapping complete ***\n";
clock();

#################################
# Depth of coverage calculation #
#################################

my $sndsixfork_manager = Parallel::ForkManager->new(6);

my @recalibrated = glob("$path2bam/*calibrated.bam"); #make array of .(re/un)calibrated.bam files
for (my $i = 0; $i < @recalibrated; $i++) { # loop through array of bamfiles
    $sndsixfork_manager->start and next;
	coverage($recalibrated[$i]); # pass to coverage sub
    $sndsixfork_manager->finish;
}

$sndsixfork_manager-> wait_all_children;

print "\n\n*** Depth of coverage calculation complete ***\n";
clock();

################## SUBROUTINES ###################

sub clock {
	my @months = qw(Jan Feb Mar Apr May Jun Jul Aug Sep Oct Nov Dec);
	my @weekDays = qw(Sun Mon Tue Wed Thu Fri Sat Sun);
	my ($second, $minute, $hour, $dayOfMonth, $month, $yearOffset, $dayOfWeek, $dayOfYear, $daylightSavings) = localtime();
	my $year = 1900 + $yearOffset;
	my $theTime = "$hour:$minute:$second, $weekDays[$dayOfWeek] $months[$month] $dayOfMonth, $year";
	print "*** $theTime ***\n\n";
}

sub usage {
	print "\n\n\n\t---------Control exome alignment ------------------\n\n";
	print "\tUsage:\tControl_exome_alignment_v1.pl [-fastq -exome]\n";
	print "\t-fastq - is required. The path to the directory containing the fastq files. NB: do not add trailing / \n";
    print "\t-exome - is required. Define an exome used from truseq, seqcapv4, vcrome2.1, sureselectv2 / \n";
	print "\tUse this script to align exome sequencing data. Requires paired end fastq.gz compressed reads\n";
	print "\tReads must be stored in your /path/2/Rawdata folder and labelled:\n\n";
	print "\t\tsamplename1.f.fastq.gz\n";
	print "\t\tsamplename1.r.fastq.gz\n";
	print "\t\tsamplename2.f.fastq.gz\n";
	print "\t\tsamplename2.r.fastq.gz\n";
	print "\t\tetc..  \n\n";
	print "\tReplace samplename with a unique ID (e.g. blood code)\n";
	print "\tContact r.bagnall\@centenary.org.au\n";
	print "\t-------------------------------------------------------------\n\n";
	exit;
}

sub bwa_mem {
    my ($current_read_pair) = shift(@_); # get pair of read names
    my @read_pair = split(" ", $current_read_pair);
    my $current_samplename = substr $read_pair[0], ($offset + 9), -11; # i.e. get IO2 from $input/Rawdata/IO2.1.fastq.gz
    print "\n\n*** Aligning $current_samplename reads ***\n";
    clock();
    
    my $return_value = system(
    "bash", "-c", q!bwa mem -M -t 6 -R ! . "'\@RG\tID:$current_samplename\tSM:$current_samplename\tPL:ILLUMINA\tLB:$current_samplename'" . " $path2ref " .
    q!<! . "(zcat $read_pair[0] | tr ' ' '_'  ) " .
    q!<! . "(zcat $read_pair[1] | tr ' ' '_'  ) " .
    q! | ! . " sed -e 's/CL:bwa.*//' " .
    q! | ! .
    q! samtools view -Su - ! .
    q! | ! . "novosort - --removeDuplicates --threads 6 --ram 10G --index --output $path2bam/$current_samplename.uncalibrated.bam"
    );
    print "\n\n*** Finished aligning $current_samplename reads ***\n";
    clock();
}

sub coverage {
	my ($current_bam) = shift(@_);# $input/tempBAMs/sample1.recalibrated.bam
	my $current_name = substr $current_bam, ($offset + 10), -17;
	print "\n\n*** Calculating coverage for $current_name ***\n";
	clock();
    
    my $coverage = system("coverageBed -abam $current_bam -b $path2targetregions -d | gzip > $input/tempBAMs/$current_name.per.base.coverage.gz");
    my $coverage1 = system("Rscript /home/shared/NGS/human/richardb/Applications/Rscripts/exome_coverage.R $input/tempBAMs/$current_name.per.base.coverage.gz");
    print "\n\n*** Finished $current_name ***\n\n";
}