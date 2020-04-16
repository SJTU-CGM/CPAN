#!/usr/bin/perl
package geneExist;

sub geneExist{
use strict;
use warnings;
use Getopt::Std;
use vars qw($opt_h);
getopts("h");

my $usage="\nUsage: cpan geneExist <data_directory> <output_directory> <min_gene_cov> <min_cds_cov>

cpan geneExist is used to determine gene presence-absence based on gene body coverage and CDS coverage.

Necessary input description:

  data_directory      <string>        This directory should contain four files: 
                                      \"Normal.merged.gene.cov.txt\", \"Normal.merged.cds.cov.txt\",
                                      \"Tumor.merged.gene.cov.txt\" and \"Tumor.merged.cds.cov.txt\".

  output_directory    <string>        Results will be output to this directory.To avoid overwriting 
                                      of existing files. We kindly request that the output_directory
                                      should not exist. It is to say, this directory will be created
                                      by the script itself.

  <min_gene_cov>      <float>         Minimum gene body coverage.

  <min_cds_cov>       <float>         Minimum CDS coverage.

Options:
     -h                               Print this usage page.

";

die $usage if @ARGV!=4;
die $usage if defined($opt_h);
my ($data_dir,$out_dir,$min_gene_cov,$min_cds_cov)=@ARGV;

#Check existence of output directory
if(-e $out_dir){
    die("Error: output directory \"$out_dir\" already exists. To avoid overwriting of existing files. We kindly request that the output directory should not exist.");
}

#Adjust directory names and create output directory
$data_dir.="/" unless($data_dir=~/\/$/);
$out_dir.="/" unless($out_dir=~/\/$/);

#Create output directory
mkdir($out_dir);

my $normal_gene_file=$data_dir."Normal.merged.gene.cov.txt";
my $normal_cds_file=$data_dir."Normal.merged.cds.cov.txt";
my $tumor_gene_file=$data_dir."Tumor.merged.gene.cov.txt";
my $tumor_cds_file=$data_dir."Tumor.merged.cds.cov.txt";

unless (-e $normal_gene_file) {
  die("Error: cannot find \"Normal.merged.gene.cov.txt\" in the data_dir.\n");
}
unless (-e $tumor_gene_file) {
  die("Error: cannot find \"Tumor.merged.gene.cov.txt\" in the data_dir.\n");
}
unless (-e $normal_cds_file) {
  die("Error: cannot find \"Normal.merged.cds.cov.txt\" in the data_dir.\n");
}
unless (-e $tumor_cds_file) {
  die("Error: cannot find \"Tumor.merged.cds.cov.txt\" in the data_dir.\n");
}

my $normal_geneExist_file=$out_dir."Normal.geneExist.txt";
my $normal_distributed_file=$out_dir."Normal.distributeGene.txt";
my $tumor_geneExist_file=$out_dir."Tumor.geneExist.txt";
my $tumor_distributed_file=$out_dir."Tumor.distributeGene.txt";

open(NGF,$normal_gene_file)||die("Error: cannot read the file: $normal_gene_file");
open(NCF,$normal_cds_file)||die("Error: cannot read the file: $normal_cds_file");
open(OUT1,">$normal_geneExist_file")||die("Error: cannot write the file: $normal_geneExist_file.\n");
open(OUT2,">$normal_distributed_file")||die("Error: cannot write the file: $normal_distributed_file.\n");
my $n=0;
while(my $gene_line=<NGF>){
    chomp $gene_line;
    my $cds_line=<NCF>;
    chomp $cds_line;
    $n++;
    if($n==1){
        print OUT1 $gene_line,"\n";
        print OUT2 $gene_line,"\n";
        next;
    }
    else{
        my @t1=split /\t/,$gene_line;
        my @t2=split /\t/,$cds_line;
        my $record=$t1[0];
        my $sum=0;
        for(my $i=1;$i<@t1;$i++){
            if($t1[$i]>=$min_gene_cov && $t2[$i]>=$min_cds_cov){
                $record.="\t1";
                $sum+=1;
            }
            else{
                $record.="\t0";
            }
        }
        print OUT1 $record."\n";
        if($sum!=$#t1){
          print OUT2 $record."\n";
        }
    }
}
close NGF;
close NCF;
close OUT1;
close OUT2;

open(TGF,$tumor_gene_file)||die("Error: cannot read the file: $tumor_gene_file");
open(TCF,$tumor_cds_file)||die("Error: cannot read the file: $tumor_cds_file");
open(OUT3,">$tumor_geneExist_file")||die("Error: cannot write the file: $tumor_geneExist_file.\n");
open(OUT4,">$tumor_distributed_file")||die("Error: cannot write the file: $tumor_distributed_file.\n");
my $k=0;
while(my $gene_line=<TGF>){
    chomp $gene_line;
    my $cds_line=<TCF>;
    chomp $cds_line;
    $k++;
    if($k==1){
        print OUT3 $gene_line,"\n";
        print OUT4 $gene_line,"\n";
        next;
    }
    else{
        my @t1=split /\t/,$gene_line;
        my @t2=split /\t/,$cds_line;
        my $record=$t1[0];
        my $sum=0;
        for(my $i=1;$i<@t1;$i++){
            if($t1[$i]>=$min_gene_cov && $t2[$i]>=$min_cds_cov){
                $record.="\t1";
                $sum+=1;
            }
            else{
                $record.="\t0";
            }
        }
        print OUT3 $record."\n"; 
        if($sum!=$#t1){
          print OUT4 $record."\n";
        }
    }
}
close TGF;
close TCF;
close OUT3;
close OUT4;
}
1;
