#!/usr/bin/perl
#Created by Duan Zhongqu at 2019/4/5;
#This script is designed to merge and clust the novel sequences from multiple samples according to the anchored positions.
package mergeNovSeq;

sub mergeNovSeq{
    use strict;
    use warnings;
    use Getopt::Std;
    use vars qw($opt_h);
    getopts("h");
    my $usage="\nUsage: cpan mergeNovSeq [options] <data_directory> <output_directoty>

cpan mergeNovSeq is used to merge the novel sequences from multiple samples.\n.

Necessary input description:

  data_directory       <string>      This directory should contain many sub-directories
                                     named by sample names, such as Sample1, Sample2,etc.
                                     In each sub-directory, there should be the files of 
                                     novel sequences (include *.two_end.fa and *.one_end.fa).
    
  output_directory     <string>      Final output files will be found in this directory. 
                                     To avoid overwriting of existing files. We kindly 
                                     request that the output_directory should not exist.
                                     It is to say, this directory will be created by the 
                                     script itself.

Options:

  -h                                 Print this usage page.
";

    die $usage if @ARGV!=2;
    die $usage if defined($opt_h);
    my ($data_dir,$out_dir)=@ARGV;

    #check existance of output directory
    if(-e $out_dir){
        die("Error: output directory \"$out_dir\" already exists. To avoid overwriting of existing files, we kindly request that the output directory should not ex
ist.\n");
    }
    
    #adjust directory names and create output directory
    $data_dir.="/" unless($data_dir=~/\/$/);
    $out_dir.="/" unless($out_dir=~/\/$/);
    mkdir($out_dir);
 
    #merge the novel sequence from multiple individuals
    my $fa_dir=$out_dir."fa/";
    mkdir($fa_dir);
    my $one_end_fa=$fa_dir."Merged_one_end.fa";
    my $two_end_fa=$fa_dir."Merged_two_end.fa";
    open(FA1,">$one_end_fa")||die("Error10: cannot write the file: $one_end_fa.\n");
    open(FA2,">$two_end_fa")||die("Error10: cannot write the file: $two_end_fa.\n");
    
    my $bed_dir=$out_dir."bed/";
    mkdir($bed_dir);
    my $one_end_bed=$bed_dir."Merged_one_end.bed";
    my $two_end_bed=$bed_dir."Merged_two_end.bed";
    open(BED1,">$one_end_bed")||die("Error10: cannot write the file: $one_end_bed.\n");
    open(BED2,">$two_end_bed")||die("Error10: cannot write the file: $two_end_bed.\n");
    
    #Read samples
    opendir(DATA,$data_dir)||die("Error03: cannot open data directory: $data_dir.\n");
    my @sample=readdir(DATA);
    close DATA;

    #process each sample
    foreach my $s (@sample){
        next if $s=~/^\./;
        print STDERR "Process sample $s ...\n";
        my $s_dir=$data_dir.$s."/";
        unless(-d $s_dir){
            print STDERR "Warning: cannot find $s in $data_dir. Not processing.\n";
            next;
        }
        opendir(SAM,$s_dir)||die("Warning: cannot open data directory: $s_dir.\n");
        my @files=readdir(SAM);
        close SAM;
        my $one_end_file="";
        my $two_end_file="";
        foreach my $f (@files){
            next if $s=~/^\./;
            if($f=~/\.one_end\.fa$/){
                $one_end_file=$s_dir.$f;
            }
            if($f=~/\.two_end\.fa$/){
                $two_end_file=$s_dir.$f;
            }
        }
        if($one_end_file eq ""){
            die("Error: cannot find the file of one-end placed novel sequences in $s_dir. Please check the directory and rerun.\n");
        }
        if($two_end_file eq ""){
            die("Error: cannot find the file of two-end placed novel sequences in $s_dir. Please check the directory and rerun.\n");
        }
        open(IN,$one_end_file)||die("Error09: cannot read the file: $one_end_file.\n");
        while(my $line=<IN>){
            chomp $line;
            print FA1 $line."\n";
            if($line =~ /^>/){
            	$line =~ s/>//g;
            	my @arr = split(/[:|\-|\t|_]/,$line);
            	if($arr[2] eq "right"){
                    $arr[2] = $arr[1]+1;
                    print BED1 "$arr[0]\t$arr[1]\t$arr[2]\t$line\n";
                }
                elsif($arr[2] eq "left"){
                    $arr[2] = $arr[1]-1;
                    print BED2 "$arr[0]\t$arr[2]\t$arr[1]\t$line\n";
                }
            }
        }
        close IN;
        open(IN,$two_end_file)||die("Error09: cannot read the file: $two_end_file.\n");
        while(my $line=<IN>){
            chomp $line;
            print FA2 $line."\n";
            if($line =~ /^>/){
            	$line =~ s/>//g;
            	my @arr = split(/[:|\-|\t|_]/,$line);
            	if($arr[1] > $arr[2]){
                    print BED2 "$arr[0]\t$arr[2]\t$arr[1]\t$line\n";
                }
                else{
                    print BED2 "$arr[0]\t$arr[1]\t$arr[2]\t$line\n";
                }
            }
        }
        close IN;
    }
    close FA1;
    close FA2;
    close BED1;
    close BED2;
}
1;
