#!/usr/bin/perl
package geneCov;

sub geneCov{
use strict;
use warnings;
use Getopt::Std;
use vars qw($opt_h $opt_n $opt_t);
getopts("hn:t:");
my $usage="\nUsage: cpan ccov [options] normal_bam_directory tumor_bam_directory sample_list output_directory gene_annotation 

cpan ccov is used to caculate gene coverages of each gene of both normal sample and tumor sample.

The script will call samtools.

Necessary input description:
    
  normal_bam_directory       <string>      This directory should contain many sub-directories
                                      	   named by normal sample names, such as Normal_1, 
                                      	   Normal_2,etc. In each sub-directory, mapping result,
                                      	   a sorted and index bam file, should exist.

  tumor_bam_directory        <string>      This directory should contain many sub-directories
                                      	   named by tumor sample names, such as Tumor_1, 
                                      	   Tumor_2,etc. In each sub-directory, mapping result,
                                      	   a sorted and index bam file, should exist. 

  sample_list                <string>      File contains the relationship of normal sample names 
                                           and tumor sample names. Each row has two columns 
                                           sepereted by \"\t\", with the first column defines 
                                           normal sample name and the second column defines tumor
                                           sample name, in one individual.

  output_directory           <string>      Results will be output to this directory.To avoid 
                                           overwriting of existing files. We kindly request
                                           that the output_directory should not exist. It is
                                           to say, this directory will be created by the 
                                           script itself.

  gene_annotation            <string>      Gene annotations in a single gtf file.



Options: 

      -h                             Print this usage page.

      -n               <int>         Minimal depth threshold to consider position covered
                                     in normal sample [default: 1].

      -t               <int>         Minimal depth threshold to consider position covered
                                     in tumor sample [default: 1].
 
";

die $usage if @ARGV!=5;
die $usage if defined($opt_h);
my ($normal_dir,$tumor_dir,$sample_list,$out_dir,$gtf_file)=@ARGV;

#Read the minimal depth threshold to consider position covered 
my $normal_minDep=1;
$normal_minDep=$opt_n if defined($opt_n);
my $tumor_minDep=1;
$tumor_minDep=$opt_t if defined($opt_t);

#Check existence of output directory
if(-e $out_dir){
    die("Error: output directory \"$out_dir\" already exists. To avoid overwriting of existing files. We kindly request that the output directory should not exist.");
}

#Adjust directory names and create output directory
$normal_dir.="/" unless($normal_dir=~/\/$/);
$tumor_dir.="/" unless($tumor_dir=~/\/$/);
$out_dir.="/" unless($out_dir=~/\/$/);

#Create output directory and sub-directories
mkdir($out_dir);
my $out_data=$out_dir."data/";
mkdir($out_data);
my $err_data=$out_dir."err/";
mkdir($err_data);

#detect perl script to caculate the coverages
my $exe_ccov="ccov.pl";
my @path=split /:/,$ENV{PATH};
my $fpflag=0;
foreach my $p (@path){
    $p.="/".$exe_ccov;
    if(-e $p && -x $p){
        $fpflag=1;
        last;
    }
}
die("Executable ccov.pl cannot be found in your PATH!\n") unless($fpflag);

#Read normal samples
opendir(DATA,$normal_dir) || die("Error03: cannot open data directory: $normal_dir\n");
my @normal_samples=readdir(DATA);
closedir DATA;

#Read tumor samples
opendir(DATA,$tumor_dir) || die("Error03: cannot open data directory: $tumor_dir\n");
my @tumor_samples=readdir(DATA);
closedir DATA;

#Process each individual
open(IN,$sample_list)||die("Error: cannot read the file: $sample_list.\n");
while(my $line=<IN>){
	chomp $line;
	my ($normal,$tumor)=split "\t",$line;
	print STDERR "Process samples $normal and $tumor\n";
	my $normal_sample_dir=$normal_dir.$normal."/";
	my $tumor_sample_dir=$tumor_dir.$tumor."/";
	unless (-d $normal_sample_dir) {
		print STDERR "Warnings: cannot find the directory of normal sample: $normal in $normal_dir.\n";
		next;
	}
	unless (-d $tumor_sample_dir) {
		print STDERR "Warnings: cannot find the directory of tumor sample: $tumor in $tumor_dir.\n";
                next;
	}
	my $normal_bam_file=$normal_sample_dir.$normal.".bam";
        unless(-e $normal_bam_file){
        	print STDERR "Warnings: cannot find normal bam file($normal_bam_file) in $normal_sample_dir: skip this individual\n";
        	next; 
        }
    	my $tumor_bam_file=$tumor_sample_dir.$normal.".bam";
    	unless(-e $tumor_bam_file){
        	print STDERR "Warnings: cannot find tumor bam file($tumor_bam_file) in $tumor_sample_dir: skip this individual\n";
        	next; 
    	}
    	my $sample_out=$out_data.$normal.".".$tumor."/";
    	mkdir($sample_out) unless(-e $sample_out);
    	my $sample_err=$err_data.$normal.".".$tumor."/";;
    	my $sta_file=$sample_out.$normal.".".$tumor.".sta";
	#generate command
    	my $com="$exe_ccov $gtf_file $normal_bam_file $tumor_bam_file $sta_file $normal_minDep $tumor_minDep\n";
    	system($com);
}
}
1;
