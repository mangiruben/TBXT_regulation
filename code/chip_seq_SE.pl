#!/usr/bin/perl

##################################################################################################################
##														##
##	Reuben_UKDRI, 2022											##
##														##
## Automated analysis of single end ChIP-seq data ##
## Runs in the terminal as perl chip_seq_SE.pl   ##

##################################################################################################################

use strict;
use warnings;
use Getopt::Long;
use File::Temp qw(tempdir tempfile);
use File::Basename;
use IO::File;
use List::Util qw(shuffle max sum);
use Math::Trig;
use File::Path qw(make_path);
use Term::ANSIColor qw(:constants);

my ($data,$GENOME_IDX, $genome, $result_path, $fastqc, $mapfiles);#declaring variables
#defaults
my $trim5=0;
my $trim3=0;
my $threads=6;
 
GetOptions(
    "g=s" => \$genome,
    "p=s" => \$threads,
    "w=s" => \$data,
    "3=s" => \$trim3,
    "5=s" => \$trim5,
    "o=s" => \$result_path,
          );
            
if (!$data ) {help();}

sub help 
{
    print GREEN"\nChiPseq pipeline on making\n", RESET;
    print GREEN "\nUsage: ChiPseq [options]\n", RESET;
    print "NOTE: This script uses Bowtie2, fastQC, Samtools and HOMER  to map reads to the target genome\n";
    print "      All the modules and genome index must be installed and available from \$PATH environment variable.\n\n";
    print "Provide the below options to start analysis:\n\n";
    print "-g         Name the genome build; eg hg19, mm10\n";
    print "-p         number of threads to use, default is 6\n";
    print "-w         state the path of the fastq files directory \n";
    print "-5         trim bases from 5'/left end reads default is 0\n";
    print "-3         trim bases from 3'/Right end reads default is 0\n";
    print "-o         results output directory, if not analysis sent to working directory \n\n";
    exit;
}

#checking Ref_genome path: replace the path to bowtie2 indices, will do a control for this
#this indices are available at epinott at RDS Imperial

my $hg19="/rds/general/user/ryaa/projects/epinott/live/Ref_genome/human/hg19_index/hg19";
my $hg38="/rds/general/user/ryaa/projects/epinott/live/Ref_genome/human/hg38_index/hg38";
my $mm10="/rds/general/user/ryaa/projects/epinott/live/Ref_genome/mouse/mm10_index/mm10";
if($genome eq 'hg19'){$GENOME_IDX =$hg19}
else{ if($genome eq 'hg38'){ $GENOME_IDX =$hg38}
else{if($genome eq 'mm10'){ $GENOME_IDX =$mm10 } 
else { 
    die print GREEN "\n\nAnalysis stopped, specify correct genome build, either hg38, hg19 or mm10, and re-submit \n\n", RESET;
     }
}
}

if($data) { 
    if(-d $data) { 
        print "\n################################################################################" .
"\n##                                                                            ##" .
"\n##                         ChIP-Seq Analysis Initiated!                       ##" .
"\n##                                                                            ##" .
"\n################################################################################\n";
        print "\n### DIRECTORY TO BEGIN ANALYSIS:\n";

        print "data directory=$data\n";
    } else { # $data is NOT a real directory
        print "\n### WARNING: DIRECTORY $data DOESN'T EXIST. Checking data in the working directory:\n";
        $data = `pwd`; 
        chomp($data);
        print "data directory=$data\n";
    }
} else {
    print "\n### Begin analysis in working directory:\n";
    $data = `pwd`; 
    chomp($data);
    print "data directory=$data\n";
}
print "Checking QCs and mapping\n";

my ($file , $sample_names, $sample_name);

my @sample_names;
my %sample_names;
opendir DIR, $data or die "cannot open folder $data check the fq please \n";
while (@sample_names = readdir DIR) {
foreach  (@sample_names) { 
my $sample_name= $_;
    if ($_ =~ /fastq.gz/) {
fastqc("$data/$_");
bowtiemap("$data/$_");
}
       }     
    }
closedir DIR;

sub fastqc 
    {
    my $file= @_;
    make_path ("$result_path/chip_results/fastqc");
    my $tmp = File::Temp->new();
    print "    $_\n";
    my $out = "$result_path/chip_results/fastqc";
    system ( "fastqc $data/$_ -o $out >> $tmp 2>&1");
    }
sub bowtiemap 
    {
    my $file= @_;
   make_path ("$result_path/chip_results/mapfiles");
   make_path ("$result_path/chip_results/bowtie_log");
    my $tmp = File::Temp->new();
    my $log ="$result_path/chip_results/bowtie_log";
    my $out ="$result_path/chip_results/mapfiles";
    #print "mapping metrics >> $_\n";
    system ("bowtie2 -p $threads -x $GENOME_IDX -U $data/$_  -5 $trim5 -3 $trim3 2> $log/$_\_log.out | samtools sort  -o $out/$_.bam ");
    }

print "making Tag directories and Bigwigs \n";

my @bamfiles = <$result_path/chip_results/mapfiles/*.bam>;
my $tmp; my $name;
# looping through all the result folders and making Tag directories
foreach  my $bamfiles(@bamfiles)
    {
        my $name = fileparse($bamfiles, '\..*');
        print "$name\n";
        make_path ("$result_path/chip_results/tags");
        my $tags  = "$result_path/chip_results/tags/$name\_tag";
        $tmp = File::Temp->new();
        system ("makeTagDirectory $tags  $bamfiles  -genome $genome >> $tmp 2>&1");
    }

my @tags = <$result_path/chip_results/tags/*_tag>;

# looping through all the result folders and making bigwig
foreach  my $tags(@tags)
    {
        my $name = basename($tags);
        make_path ("$result_path/chip_results/bedGraphs");
        my $bedgraph ="$result_path/chip_results/bedGraphs/$name\.bedGraph";
        my $c_peaks = "$result_path/chip_results/tags/$name/";
        my $tmp = File::Temp->new();
        system ("makeUCSCfile  $c_peaks  -o $bedgraph  -name $name -res 1  -fsize 5e7 -color 0,0,0 >> $tmp 2>&1");
    }

my $dat = basename($result_path);
print   "\n\nAnalysis completed, see results in ", GREEN "$data", RESET " output folder you provided \n\n";

