#!/usr/bin/perl
#Created by Duan Zhongqu at 2019/4/5;
#This script is designed to cluster the novel sequences according to the anchored positions.
package clusterNovSeq;

sub clusterNovSeq{
    use strict;
    use warnings;
    use Getopt::Std;
    use vars qw($opt_h);
    getopts("h");
    my $usage="\nUsage: cpan clusterNovSeq [options] <data_directory> <output_directoty>

cpan clusterNovSeq is used to cluster the novel sequences according to the anchored positions.\n.

Necessary input description:

  data_directory       <string>      This directory should contain two sub-directories named
                                     \"fa\" and \"bed\", respectively. The \"fa\" directory 
                                     contains two files named \"Merged_one_end.fa\" and 
                                     \"Merged_two_end.fa\", respectively. The \"bed\" 
                                     directory contains two files named \"Merged_one_end.bed\"
                                     and \"Merged_two_end.bed\". 
    
  output_directory     <string>      Both final output files and intermediate results will
                                     be found in this directory. To avoid overwriting of 
                                     existing files, we kindly request that the output 
                                     directory should not exist. It is to say, this 
                                     directory will be created by the script itself.

Options:

  -h                                 Print this usage page.
";

    die $usage if @ARGV!=2;
    die $usage if defined($opt_h);
    my ($data_dir,$out_dir)=@ARGV;

    #check existance of output directory
    if(-e $out_dir){
        die("Error: output directory \"$out_dir\" already exists. To avoid overwriting of existing files, we kindly request that the output directory should not exist.\n");
    }
    
    #adjust directory names and create output directory
    $data_dir.="/" unless($data_dir=~/\/$/);
    $out_dir.="/" unless($out_dir=~/\/$/);
    mkdir($out_dir);
 
    my $seq_dir=$data_dir."fa/";
    my $bed_dir=$data_dir."bed/";
    unless(-d $seq_dir){
        die("Error: cannot find the \"fa\" directory in: $data_dir.\n");
    }
    unless(-d $bed_dir){
        die("Error: cannot find the \"bed\" directory in: $data_dir.\n");
    }

    my $one_end_fa=$seq_dir."Merged_one_end.fa";
    my $two_end_fa=$seq_dir."Merged_two_end.fa";
    my $one_end_bed=$bed_dir."Merged_one_end.bed";
    my $two_end_bed=$bed_dir."Merged_two_end.bed";
    unless(-e $one_end_fa){
        die("Error: cannot find the \"Merged_one_end.fa\" directory in: $seq_dir.\n");
    }
    unless(-e $two_end_fa){
        die("Error: cannot find the \"Merged_two_end.fa\" directory in: $seq_dir.\n");
    }
    unless(-e $one_end_bed){
        die("Error: cannot find the \"Merged_one_end.bed\" directory in: $bed_dir.\n");
    }
    unless(-e $two_end_bed){
        die("Error: cannot find the \"Merged_two_end.bed\" directory in: $bed_dir.\n");
    }
    
    my %two_end_seq;
    my $name;
    open(IN,$two_end_fa)||die("Error: cannot read the file: $two_end_fa.\n");
    while(my $line=<IN>){
        chomp $line;
        if($line=~/>/){
            $name=$line;
        }else{
            $two_end_seq{$name}=$line;
        }
    }
    close IN;
    my %one_end_seq;
    open(IN,$one_end_fa)||die("Error: cannot read the file: $one_end_fa.\n");
    while(my $line=<IN>){
        chomp $line;
        if($line=~/>/){
            $name=$line;
        }else{
            $one_end_seq{$name}=$line;
        }
    }
    close IN;
    
    #two_end placed novel sequences
    my $sorted_two_bed=$out_dir."Merged_two_end.sorted.bed";
    my $merged_two_bed=$out_dir."Merged_two_end.unique.bed";
    system("sort -k1,1 -k2n,2 $two_end_bed >$sorted_two_bed");
    system("bedtools merge -d 10 -i $sorted_two_bed >$merged_two_bed");

    open(IN,$merged_two_bed)||die("Error09: cannot read the file: $merged_two_bed.\n");
    my %two_bed_position;
    while(<IN>){
        chomp;
        my @string=split /\s+/,$_;
        $two_bed_position{$string[0]."_".$string[1]."_".$string[2]}=0;
    }
    close IN;

    my $two_end_fa_dir=$out_dir."two_end/";
    mkdir($two_end_fa_dir);
    my $fa_dir=$two_end_fa_dir."fa/";
    mkdir($fa_dir);
    while((my $key, my $value)=each %two_end_seq){
        my @string = split(/[:|\-|_]/,$key);
        my $chr=$string[0];
        my $start=$string[1];
        my $end=$string[2];
        while((my $k,my $v)=each %two_bed_position){
            my @tmp=split "_",$k;
            if(($tmp[0] eq $chr)&&($tmp[1]<=$start)&&($tmp[2]>=$end)){
                my $file=$fa_dir.$tmp[0].".".$tmp[1]."-".$tmp[2].".fa";
                open(OUT,">>$file")||die("Error10: cannot write the file: $file.\n");
                my $seq=substr($value,100,(length($value)-200));
                print OUT ">".$key."\n".$seq."\n";
                close OUT;
            }
        }  
    }

    opendir(FILE,$fa_dir)||die("Error03: cannot open data directory: $fa_dir.\n");
    my @fa_files=readdir(FILE);
    close FILE;
    foreach my $sample (@fa_files){
        next if $sample=~/^\./;
        my $out=substr($sample,0,length($sample)-3);
        my $input_file=$fa_dir.$sample;
        open(IN,$input_file)||die("Error09: cannot read the file: $input_file.\n");
        my %seqs;
        while(my $line=<IN>){
            chomp $line;
            if($line=~/^>/){
                $line=~s/^>//g;
                my $seq=<IN>;
                chomp $seq;
                $seqs{$line}=$seq;
            }
        }
        close IN;
        my $max_seq_name="";
        my $max_len=0;
        while((my $key,my $value)=each %seqs){
            if(length($value)>$max_len){
                $max_seq_name=$key;
                $max_len=length($value);
            }
        }
        my $ref_dir=$two_end_fa_dir."ref/";
        mkdir($ref_dir);
        my $ref_file=$ref_dir.$out.".ref.fa";
        open(OUT,">$ref_file")||die("Error10: cannot write the file: $ref_file.\n");
        print OUT ">".$max_seq_name."\n".$seqs{$max_seq_name}."\n";
        my $nucmer_dir=$two_end_fa_dir."nucmer/";
        mkdir($nucmer_dir);
        my $nucmer_out=$nucmer_dir.$out;
        system("nucmer -coords -maxmatch -p $nucmer_out $ref_file $input_file");
        close OUT;
        my $coords_file=$nucmer_dir.$out.".coords";
        my %align;
        my %noalign;
        open(IN,$coords_file)||die("Error09: cannot read the file: $coords_file.\n");
        my $i=0;
        while(my $line=<IN>){
            $i+=1;
            if($i>5){
                chomp $line;
                $line=~s/^ +//g;
                my @string=split /\s+/,$line;
                $align{$string[12]}=1;
            }
        }
        close IN;
        my $cluster_dir=$two_end_fa_dir."cluster/";
        mkdir($cluster_dir);
        my $cluster_file=$cluster_dir.$out.".cluster.fa";
        my $nocluster_dir=$two_end_fa_dir."nocluster/";
        mkdir($nocluster_dir);
        my $no_cluster_file=$nocluster_dir.$out.".nocluster.fa";
        open(OUT1,">$cluster_file")||die("Error10: cannot write the file: $cluster_file.\n");
        open(OUT2,">$no_cluster_file")||die("Error10: cannot write the file: $no_cluster_file.\n");
        while((my $key, my $value)=each %seqs){
            if(exists($align{$key})){
                print OUT1 ">".$key."\n".$two_end_seq{$key}."\n";
            }
            else{
                print OUT2 ">".$key."\n".$two_end_seq{$key}."\n";
            }
        }
        close OUT1;
        close OUT2;
    } 

    my $sorted_one_bed=$out_dir."Merged_one_end.sorted.bed";
    my $merged_one_bed=$out_dir."Merged_one_end.unique.bed";
    system("sort -k1,1 -k2n,2 $one_end_bed >$sorted_one_bed");
    system("bedtools merge -d 10 -i $sorted_one_bed >$merged_one_bed");

    #delete the one-end placed novel sequences located within 100bp of a two-end placed novel sequence.
    open(IN,$merged_one_bed)||die("Error09: cannot read the file: $merged_two_bed.\n");
    my %one_bed_position;
    while(<IN>){
        chomp;
        my @string=split /\s+/,$_;
        $one_bed_position{$string[0]."_".$string[1]."_".$string[2]}=0;
    }
    close IN;
    my $final_one_bed=$out_dir."Merged_one_end.final.bed";
    open(OUT,">$final_one_bed")||die("Error10: cannot write the file: $final_one_bed.\n");
    while((my $key,my $value)=each %one_bed_position){
        my ($chr1,$start1,$end1)=split /_/,$key;
        my $flag=0;
        while(my ($k,$v)=each %two_bed_position){
            my ($chr2,$start2,$end2)=split /_/,$k;
            if(($chr1 eq $chr2)&&($start1>=$start2-100)&&($end1<=$end2+100)){
                $flag=1;
            }
        }
        if($flag==0){
            print OUT "$chr1\t$start1\t$end1\n";
        }
        else{
            delete $one_bed_position{$key};
        }
    }
     
    my $one_end_fa_dir=$out_dir."one_end/";
    mkdir($one_end_fa_dir);
    my $one_fa_dir=$one_end_fa_dir."fa/";
    mkdir($one_fa_dir);
    while((my $key, my $value)=each %one_end_seq){
        my @string = split(/[:|\-|_]/,$key);
        my $chr1=$string[0];
        my $position=$string[1];
        my $label=$string[2];
        while((my $k,my $v)=each %one_bed_position){
            my ($chr2,$start,$end)=split "_",$k;
            if(($chr2 eq $chr1)&&($start<=$position)&&($end>=$position)){
                my $file=$one_fa_dir.$chr1.".".$start."-".$end.".fa";
                open(OUT,">>$file")||die("Error10: cannot write the file: $file.\n");
                my $seq;
                if($label eq "right"){
                    $seq=substr($value,100);
                }else{
                    $seq=substr($value,0,length($value)-100); 
                }
                print OUT ">".$key."\n".$seq."\n";
                close OUT;
            }
        }
    }
    
    opendir(FILE,$one_fa_dir)||die("Error03: cannot open data directory: $one_fa_dir.\n");
    my @one_end_fa_files=readdir(FILE);
    close FILE;
    foreach my $sample (@one_end_fa_files){
        next if $sample=~/^\./;
        #my $out=substr($sample,0,length($sample)-3);
        my $input_file=$one_fa_dir.$sample;
        open(IN,$input_file)||die("Error09: cannot read the file: $input_file.\n");
        my %seqs;
        while(my $line=<IN>){
            chomp $line;
            if($line=~/^>/){
                $line=~s/^>//g;
                my $seq=<IN>;
                chomp $seq;
                $seqs{$line}=$seq;
            }
        }
        close IN;
        my $max_seq_name="";
        my $max_len=0;
        while((my $key,my $value)=each %seqs){
            if(length($value)>$max_len){
                $max_seq_name=$key;
                $max_len=length($value);
            }
        }
        my @array = split(/[:|\-|_]/,$max_seq_name);
        my $out=$array[0].".".$array[1];
        my $ref_dir=$one_end_fa_dir."ref/";
        mkdir($ref_dir);
        my $ref_file=$ref_dir.$out.".ref.fa";
        open(OUT,">$ref_file")||die("Error10: cannot write the file: $ref_file.\n");
        print OUT ">".$max_seq_name."\n".$seqs{$max_seq_name}."\n";
        my $nucmer_dir=$one_end_fa_dir."nucmer/";
        mkdir($nucmer_dir);
        my $nucmer_out=$nucmer_dir.$out;
        system("nucmer -coords -maxmatch -p $nucmer_out $ref_file $input_file");
        close OUT;
        my $coords_file=$nucmer_dir.$out.".coords";
        my %align;
        my %noalign;
        open(IN,$coords_file)||die("Error09: cannot read the file: $coords_file.\n");
        my $i=0;
        while(my $line=<IN>){
            $i+=1;
            if($i>5){
                chomp $line;
                $line=~s/^ +//g;
                my @string=split /\s+/,$line;
                $align{$string[12]}=1;
            }
        }
        close IN;
        my $cluster_dir=$one_end_fa_dir."cluster/";
        mkdir($cluster_dir);
        my $cluster_file=$cluster_dir.$out.".cluster.fa";
        my $nocluster_dir=$one_end_fa_dir."nocluster/";
        mkdir($nocluster_dir);
        my $no_cluster_file=$nocluster_dir.$out.".nocluster.fa";
        open(OUT1,">$cluster_file")||die("Error10: cannot write the file: $cluster_file.\n");
        open(OUT2,">$no_cluster_file")||die("Error10: cannot write the file: $no_cluster_file.\n");
        while((my $key, my $value)=each %seqs){
            if(exists($align{$key})){
                print OUT1 ">".$key."\n".$one_end_seq{$key}."\n";
            }
            else{
                print OUT2 ">".$key."\n".$one_end_seq{$key}."\n";
            }
        }
        close OUT1;
        close OUT2;
        if(-s $no_cluster_file){
            open(IN,$no_cluster_file)||die("Error09: cannot read the file: $no_cluster_file.\n");
            my %seqs;
            while(my $line=<IN>){
                chomp $line;
                if($line=~/^>/){
                    $line=~s/^>//g;
                    my $seq=<IN>;
                    chomp $seq;
                    $seqs{$line}=$seq;
                }
            }
            close IN;
            my $max_seq_name="";
            my $max_len=0;
            while((my $key,my $value)=each %seqs){
                if(length($value)>$max_len){
                    $max_seq_name=$key;
                    $max_len=length($value);
                }
            }
            my @array = split(/[:|\-|_]/,$max_seq_name);
            my $out=$array[0].".".$array[1];
            my $ref_file=$ref_dir.$out.".ref.fa";
            open(OUT,">$ref_file")||die("Error10: cannot write the file: $ref_file.\n");
            print OUT ">".$max_seq_name."\n".$seqs{$max_seq_name}."\n";
            my $nucmer_out=$nucmer_dir.$out;
            system("nucmer -coords -maxmatch -p $nucmer_out $ref_file $input_file");
            close OUT;
            my $coords_file=$nucmer_dir.$out.".coords";
            my %align;
            my %noalign;
            open(IN,$coords_file)||die("Error09: cannot read the file: $coords_file.\n");
            my $i=0;
            while(my $line=<IN>){
                $i+=1;
                if($i>5){
                    chomp $line;
                    $line=~s/^ +//g;
                    my @string=split /\s+/,$line;
                    $align{$string[12]}=1;
                }
            }
            close IN;
            my $cluster_file=$cluster_dir.$out.".cluster.fa";
            my $no_cluster_file=$nocluster_dir.$out.".nocluster.fa";
            open(OUT1,">$cluster_file")||die("Error10: cannot write the file: $cluster_file.\n");
            open(OUT2,">$no_cluster_file")||die("Error10: cannot write the file: $no_cluster_file.\n");
            while((my $key, my $value)=each %seqs){
                if(exists($align{$key})){
                    print OUT1 ">".$key."\n".$one_end_seq{$key}."\n";
                }
                else{
                    print OUT2 ">".$key."\n".$one_end_seq{$key}."\n";
                }
            }
            close OUT1;
            close OUT2;
        }
    }
}
1;
