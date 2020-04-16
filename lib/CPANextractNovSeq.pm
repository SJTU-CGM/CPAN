#!/usr/bin/perl
#Created by Duan Zhongqu at 2019/1/1.
#This script is designed to extract novel sequences from partially unaligned contigs:
#Two-end placed novel sequences and one-end placed novel sequences
package extractNovSeq;

sub extractNovSeq{
    use strict;
    use warnings;
    use Getopt::Std;
    use vars qw($opt_h $opt_l $opt_i);
    getopts("hl:i:");
    my $usage="\nUsage: cpan extractNovSeq [options] <data_directory> <output_directoty>

cpan extractNovSeq is used to extract novel sequences from partially unaligned contigs.\n.

Necessary input description:

  data_directory       <string>      This directory should contain many sub_directories
                                     named by sample names, such as Sample1, Sample2, etc.
                                     In each sub-directory, two files (*.partially.contig
                                     and *.partially.coordinate), which derived from QUAST
                                     results, should exist.

  output_directory     <string>      Both final output files and intermediate results 
                                     will be found in this directory. To avoid 
                                     overwriting of existing files. We kindly request
                                     that the output_directory should not exist. It is
                                     to say, this directory will be created by the 
                                     script itself.

Options: 

      -h                             Print this usage page.

      -i               <int>         The threshold of identity to filer alignment hits.
                                     default: 95%.

      -l               <int>         The threshold of alignemnt length to filter alignment
                                     hits. default 500 bp.
";

    die $usage if @ARGV!=2;
    die $usage if defined($opt_h);
    my ($input_dir,$out_dir)=@ARGV;
    
    #check existance of output directory
    if(-e $out_dir){
        die("Error: output directory \"$out_dir\" already exists.\nTo avoid overwriting of existing files, we kindly request that the output directory should not exist.\n");
    }
    
    #adjust directory names and create output directory
    $input_dir.="/" unless($input_dir=~/\/$/);
    $out_dir.="/" unless($out_dir=~/\/$/);
    mkdir($out_dir);
 
    #get the threshold of identity to filer alignment hits
    my $min_identity=95;
    $min_identity=$opt_i if defined $opt_i;
    
    #get the threshold of alignment length to filer alignment hits
    my $min_length=500;
    $min_length=$opt_l if defined $opt_l;
    
    #Read samples
    opendir(DATA,$input_dir)||die("Error03: cannot open data directory: $input_dir.\n");
    my @sample=readdir(DATA);
    close DATA;

    #process each sample
    foreach my $s (@sample){
        next if $s=~/^\./;
        print STDERR "Process sample $s ...\n";
        my $s_dir=$input_dir.$s."/";
        unless(-d $s_dir){
            print STDERR "Warning: cannot find $s in $input_dir. Not processing.\n";
            next;
        }
        opendir(SAM,$s_dir)||die("Warning: cannot open data directory: $s_dir.\n");
        my @files=readdir(SAM);
        close SAM;
        my $fasta_file="";
        my $coord_file="";
        foreach my $f (@files){
            next if $s=~/^\./;
            #fasta file
            if($f=~/\.partially\.contig$/){
                $fasta_file=$s_dir.$f;
            }
            #coordinate file
            if($f=~/\.partially\.coordinate$/){
                $coord_file=$s_dir.$f;
            }
        }
        if($fasta_file eq ""){
            die("Error: cannot find the fasta file of partially unaligned contigs in $s_dir. Please check the directory and rerun.\n");
        }
        if($coord_file eq ""){
            die("Error: cannot find the coordinate file of partially unaligned contigs in $s_dir. Please check the directory and rerun.\n");
        }
        #read the fasta file
        open(IN,$fasta_file)||die("Error09: cannot read the file: $fasta_file.\n");
        my %seq;
        my $name="";
        while(my $line=<IN>){
            chomp $line;
            if($line=~/^>(.+)$/){
                $name=$1;
                $seq{$name}="";
            }
            else{
                $seq{$name}.=$line;
            }
        }
        close IN;
        #read the coords file
        open(IN,$coord_file)||die("Cannot read the file: $coord_file.\n");
        readline IN;
        my %ref_pos;
        my %query_pos;
        my %ref_chr;
        while(my $line=<IN>){
            chomp $line;    
            my @string=split /\s+/,$line;
            #the alignment hit shold be meet one of the flowing two criterions:
            #aligned length >= defined_min_length, default: 500bp;
            #identity >= defined_min_identity, default: 95%;
            #alternatively, all alignment hits with identity >= 99% and length >= 100bp were also considered.
            if((($string[6]>=$min_length)&&($string[8]>=$min_identity))||(($string[6]>=100)&&($string[8]>=99))){
                if(!exists($ref_pos{$string[0]})){
                    $ref_pos{$string[0]}=$string[4].",".$string[5];
                    $ref_chr{$string[0]}=$string[1];
                    $query_pos{$string[0]}=$string[2].",".$string[3];
                }
                else{
                    $ref_pos{$string[0]}.="\t".$string[4].",".$string[5];
                    $ref_chr{$string[0]}.="\t".$string[1];
                    $query_pos{$string[0]}.="\t".$string[2].",".$string[3];
                }
            }
        }
        close IN;
        #two types of placed novel sequences: two-end placed novel sequences and one-end placed novel sequences
        my $sample_dir=$out_dir.$s."/";
        mkdir($sample_dir);
        my $one_end_file=$sample_dir.$s.".one_end.fa";
        my $two_end_file=$sample_dir.$s.".two_end.fa";
        open(ONE,">$one_end_file")||die("Error10: cannot write the file: $one_end_file.\n");
        open(TWO,">$two_end_file")||die("Error10: cannot write the file: $two_end_file.\n");
        while((my $key, my $value)=each %ref_pos){
            my @string=split "\t",$value;
            my $len=@string;
            ##one_end placed novel sequences
            if($len==1){
                my @tmp=split "_",$key;
                my $length=$tmp[1];
                my @tmp1=split ",",$ref_pos{$key};
                my $ref_start=$tmp1[0];
                my $ref_end=$tmp1[1];
                my @tmp2=split ",",$query_pos{$key};
                my $query_start=$tmp2[0];
                my $query_end=$tmp2[1];
                my $chr=$ref_chr{$key};
                my $novel_seq;
                if(($query_start<$query_end)&&(abs($query_start-1)<=100)){
                    #add 100-bp flanking sequence for each novel sequence
                    $novel_seq=substr($seq{$key},$query_end-100);
                    print ONE ">".$chr.":".$ref_end."_right"."_".$s."_".$key.":".$query_end."-".$length."\n".$novel_seq."\n";
                }
                elsif(($query_start>$query_end)&&(abs($query_end-1)<=100)){
                    $novel_seq=substr($seq{$key},$query_start-100);
                    print ONE ">".$chr.":".$ref_start."_right"."_".$s."_".$key.":".$query_start."-".$length."\n".$novel_seq."\n";
                }
                elsif((($query_start<$query_end))&&(abs($length-$query_end)<=100)){
                    $novel_seq=substr($seq{$key},0,($query_start-1+100));
                    print ONE ">".$chr.":".$ref_start."_left"."_".$s."_".$key.":1-".($query_start-1)."\n".$novel_seq."\n";
                }
                elsif((($query_start>$query_end))&&(abs($length-$query_start)<=100)){
                    $novel_seq=substr($seq{$key},0,($query_end-1+100));
                    print ONE ">".$chr.":".$ref_end."_left"."_".$s."_".$key.":1-".($query_end-1)."\n".$novel_seq."\n";
                }
            }
            ##two-end placed novel sequences
            elsif($len==2){
                my @tmp=split "_",$key;
                my $length=$tmp[1];
                my @tmp1=split "\t",$ref_pos{$key};
                my @tmp11=split ",",$tmp1[0];
                my @tmp12=split ",",$tmp1[1];
                my $ref1_start=$tmp11[0];
                my $ref1_end=$tmp11[1];
                my $ref2_start=$tmp12[0];
                my $ref2_end=$tmp12[1];
                my @tmp2=split "\t",$query_pos{$key};
                my @tmp21=split ",",$tmp2[0];
                my @tmp22=split ",",$tmp2[1];
                my $query1_start=$tmp21[0];
                my $query1_end=$tmp21[1];
                my $query2_start=$tmp22[0];
                my $query2_end=$tmp22[1];
                my @chr=split "\t",$ref_chr{$key};
                my $chr1=$chr[0];
                my $chr2=$chr[1];
                #these two anchored positions shold be on the same chromosome and the distance shold be smaller than 100 kb
                if(($chr1 eq $chr2)&&(abs($ref2_start-$ref1_end)<=100000)){
                    #the query sequence is forward mapping to the reference genome and the length of novel sequence shold be larger than 500bp
                    if(($query1_start<$query1_end)&&($query2_start<$query2_end)&&(($query2_start-$query1_end)>=500)){
                        my $novel_seq=substr($seq{$key},($query1_end-100),($query2_start-$query1_end-1+200));
                        print TWO ">".$chr1.":".$ref1_end."-".$ref2_start."_".$s."_".$key.":".($query1_end+1)."-".($query2_start-1)."\n".$novel_seq."\n";
                    }
                    #the query sequence is reverse mapping to the reference genome and the length of novel sequence shold be larger than 500bp
                    elsif(($query1_start>$query1_end)&&($query2_start>$query2_end)&&(($query2_end-$query1_start)>=500)){
                        my $novel_seq=substr($seq{$key},($query1_start-100),($query2_end-$query1_start-1+200));
                        print TWO ">".$chr1.":".$ref2_end."-".$ref1_start."_".$s."_".$key.":".($query1_start+1)."-".($query2_end-1)."\n".$novel_seq."\n";
                    }
                }
            }
            #the contig with more than two anchored positions are considered as one-end placed novel sequences
            else{
                my @tmp=split "_",$key;
                my $length=$tmp[1];
                my @tmp1=split "\t",$value;
                my @tmp2=split "\t",$query_pos{$key};
                my @tmp3=split "\t",$ref_chr{$key};
                my $num=@tmp1;
                for(my $i=1;$i<$num;$i+=1){
                    my @tmp11=split ",",$tmp1[$i-1];
                    my $ref_start=$tmp11[0];
                    my $ref_end=$tmp11[1];
                    my @tmp21=split ",",$tmp2[$i-1];
                    my $query_start=$tmp21[0];
                    my $query_end=$tmp21[1];
                    my $chr=$tmp3[$i-1];
                    my $novel_seq;
                    if(($query_start<$query_end)&&(abs($query_start-1)<=100)){
                        $novel_seq=substr($seq{$key},$query_end-100);
                        print ONE ">".$chr.":".$ref_end."_right"."_".$s."_".$key.":".$query_end."-".$length."\n".$novel_seq."\n";
                    }
                    elsif(($query_start>$query_end)&&(abs($query_end-1)<=100)){
                        $novel_seq=substr($seq{$key},$query_start-100);
                        print ONE ">".$chr.":".$ref_start."_right"."_".$s."_".$key.":".$query_start."-".$length."\n".$novel_seq."\n";
                    }
                    elsif(($query_start<$query_end)&&(abs($length-$query_end)<=100)){
                        $novel_seq=substr($seq{$key},0,($query_start-1+100));
                        print ONE ">".$chr.":".$ref_start."_left"."_".$s."_".$key.":1-".($query_start-1)."\n".$novel_seq."\n";
                    }
                    elsif(($query_start>$query_end)&&(abs($length-$query_start)<=100)){
                        $novel_seq=substr($seq{$key},0,($query_end-1+100));
                        print ONE ">".$chr.":".$ref_end."_left"."_".$s."_".$key.":1-".($query_end-1)."\n".$novel_seq."\n";
                    }
                }
            }
        }
    }
}
1;
