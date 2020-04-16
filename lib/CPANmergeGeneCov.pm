#!/usr/bin/perl
package mergeGeneCov;

sub mergeGeneCov{
use strict;
use warnings;
use Getopt::Std;
use vars qw($opt_h);
getopts("h");

my $usage="\nUsage: cpan mergeGeneCov  <data_directory> <output_directory>

cpan mergeGeneCov is used to merge the gene and CDS coverage from multiple individual genomes.

Necessary input description:

  data_directory          <string>    This directory should contain many sub-directories
                                      named by sample names, such as Normal1.Tumor1, 
                                      Normal2.Tumor2,etc. In each sub-directory, coverage
                                      file ended by \".sta\", should exist.

  output_directory        <string>    Results will be output to this directory.To avoid 
                                      overwriting of existing files. We kindly request
                                      that the output_directory should not exist. It is
                                      to say, this directory will be created by the 
                                      script itself.

Options:
     -h                               Print this usage page.

";

die $usage if @ARGV!=2;
die $usage if defined($opt_h);
my ($data_dir,$out_dir)=@ARGV;

#Check existence of output directory
if(-e $out_dir){
    die("Error: output directory \"$out_dir\" already exists. To avoid overwriting of existing files. We kindly request that the output directory should not exist.");
}

#Adjust directory names and create output directory
$data_dir.="/" unless($data_dir=~/\/$/);
$out_dir.="/" unless($out_dir=~/\/$/);
mkdir($out_dir);

my $normal_gene_out=$out_dir."Normal.merged.gene.cov.txt";
my $normal_cds_out=$out_dir."Normal.merged.cds.cov.txt";
my $tumor_gene_out=$out_dir."Tumor.merged.gene.cov.txt";
my $tumor_cds_out=$out_dir."Tumor.merged.cds.cov.txt";

#Read samples
opendir(DATA,$data_dir) || die("Error03: cannot open data directory: $data_dir\n");
my @samples=readdir(DATA);
closedir DATA;

my %normal_gene;
my %normal_cds;
my %tumor_gene;
my %tumor_cds;
my $header="Gene";
my %gene_index;
my $k=1;
#process each sample
foreach my $sample (@samples){
    next if $sample=~/^\./;
    my $sample_dir=$data_dir.$sample."/";
    next unless(-d $sample_dir);
    
    #obtain *.sta file within the sample directory
    my $sta_file=$sample_dir.$sample.".sta";
    unless(-e $sta_file){
        print STDERR "Warnings: cannot find coverage file($sta_file) in $sample_dir: skip this sample\n";
        next; 
    }
    $header.="\t".$sample;
    open(IN,$sta_file)||die("Error: cannot read the file: $sta_file.\n");
    while(my $line=<IN>){
    	chomp $line;
    	next if $line=~/^\#/;
        my @string=split /\t/,$line;
        if(exists($normal_gene{$string[0]})){
        	$normal_gene{$string[0]}.="\t".$string[4];
        	$tumor_gene{$string[0]}.="\t".$string[5];
        	$normal_cds{$string[0]}.="\t".$string[8];
        	$tumor_cds{$string[0]}.="\t".$string[9];
        }else{
        	$normal_gene{$string[0]}=$string[4];
        	$tumor_gene{$string[0]}=$string[5];
        	$normal_cds{$string[0]}=$string[8];
        	$tumor_cds{$string[0]}=$string[9];
        	$gene_index{$string[0]}=$k;
        	$k+=1;
        }
    }
    close IN;
}

open(OUT1,">$normal_gene_out")||die("Error: cannot write the file: $normal_gene_out.\n");
open(OUT2,">$tumor_gene_out")||die("Error: cannot write the file: $tumor_gene_out.\n");
open(OUT3,">$normal_cds_out")||die("Error: cannot write the file: $normal_cds_out.\n");
open(OUT4,">$tumor_cds_out")||die("Error: cannot write the file: $tumor_cds_out.\n");
print OUT1 $header."\n";
print OUT2 $header."\n";
print OUT3 $header."\n";
print OUT4 $header."\n";

foreach my $key ( sort { $gene_index{$a} <=> $gene_index{$b} } keys %gene_index ) {
	print OUT1 $key."\t".$normal_gene{$key}."\n";
	print OUT2 $key."\t".$tumor_gene{$key}."\n";
	print OUT3 $key."\t".$normal_cds{$key}."\n";
	print OUT4 $key."\t".$tumor_cds{$key}."\n";
}
close OUT1;
close OUT2;
close OUT3;
close OUT4;
}
1;
