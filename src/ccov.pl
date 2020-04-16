#!/usr/bin/perl
use strict;
use warnings;

use Getopt::Std;
use vars qw($opt_h);
getopts("h");
my $usage="\nUsage: cpan ccov gtf_file normal_bam tumor_bam out_file normal_min_cov tumor_min_cov

cpan ccov is used to caculate the gene body coverage and CDS coverage of the protein_coding genes in .

Necessary input description:

  gtf_file             <string>      Anotation file. 
    
  normal_bam           <string>      The sorted and index bam file of normal sample.

  tumor_bam            <string>      The sorted and index bam file of tumor sample.

  out_file             <string>      Output file.

  n_minDep             <int>         Minimal depth threshold to consider position covered. 
                                     Default: 1.
      
  t_minDep             <int>         Minimal depth threshold to consider position covered. 
                                     Default: 1.

Options:
  
  -h                                 Print this Usage.

";

die $usage if @ARGV!=6;
die $usage if defined($opt_h);
my ($gff_file,$normal_bam,$tumor_bam,$out_file,$n_minDep,$t_minDep)=@ARGV;

open(IN,$gff_file)||die("Error: cannot read the file: $gff_file.\n");
my %transcripts;
my %transcripts_index;
my $k=1;
my %cdss;
while(my $line=<IN>){
    chomp $line;
    my ($chr,$source,$type,$start,$end,$score,$sym,$phase,$record)=split /\t/,$line;
    $record=~/gene_id \"([^\"]+)\"; transcript_id \"([^\"]+)\"/;
    my $gene_id=$1;
    my $tran_id=$2;
    if($type eq "transcript"){
      $transcripts{$tran_id}=$chr."_".$start."_".$end;
      $transcripts_index{$tran_id}=$k;
      $k+=1;
    }
    if($type eq "CDS"){
      if(exists($cdss{$tran_id})){
        $cdss{$tran_id}.="\t".$chr."_".$start."_".$end;
      }else{
        $cdss{$tran_id}=$chr."_".$start."_".$end;
      }
    }
}
close IN;

open(OUT,">$out_file")||die("Error: cannot read the file: $out_file.\n");
print OUT "#name\tchr\tbegin\tend\tnormal_gene_coverage\ttumor_gene_coverage\tnormal_gene_depth\ttumor_gene_depth\tnormal_transcript_coverage\ttumor_transcript_coverage\tnormal_transcript_depth\ttumor_transcript_depth\n";
foreach my $key ( sort { $transcripts_index{$a} <=> $transcripts_index{$b} } keys %transcripts_index ) {
    my $value=$transcripts{$key};
    my ($chr,$start,$end)=split /_/,$value;
    my $gene_length=$end-$start+1;
    my ($normal_gene_sum,$normal_gene_base,$tumor_gene_sum,$tumor_gene_base)=@{depth($normal_bam,$tumor_bam,$chr,$start,$end,$n_minDep,$t_minDep)};
    my $normal_gene_depth=$normal_gene_sum/$gene_length;
    my $normal_gene_coverage=$normal_gene_base/$gene_length;
    my $tumor_gene_depth=$tumor_gene_sum/$gene_length;
    my $tumor_gene_coverage=$tumor_gene_base/$gene_length;
    my $normal_transcript_sum=0;
    my $normal_transcript_base=0;
    my $tumor_transcript_sum=0;
    my $tumor_transcript_base=0;
    my $transcript_length=0;
    my $cds=$cdss{$key};
    my @string=split "\t",$cds;
    foreach my $str(@string){
        my ($cds_chr,$cds_start,$cds_end)=split /_/,$str;
        $transcript_length+=$cds_end-$cds_start+1;
        my ($normal_cds_sum,$normal_cds_base,$tumor_cds_sum,$tumor_cds_base)=@{depth($normal_bam,$tumor_bam,$cds_chr,$cds_start,$cds_end,$n_minDep,$t_minDep)};
        $normal_transcript_sum+=$normal_cds_sum;
        $normal_transcript_base+=$normal_cds_base;
	$tumor_transcript_sum+=$tumor_cds_sum;
	$tumor_transcript_base+=$tumor_cds_base;
    }
    my $normal_transcript_depth=$normal_transcript_sum/$transcript_length;
    my $normal_transcript_coverage=$normal_transcript_base/$transcript_length;
    my $tumor_transcript_depth=$tumor_transcript_sum/$transcript_length;
    my $tumor_transcript_coverage=$tumor_transcript_base/$transcript_length;
    print OUT "$key\t$chr\t$start\t$end\t$normal_gene_coverage\t$tumor_gene_coverage\t$normal_gene_depth\t$tumor_gene_depth\t$normal_transcript_coverage\t$tumor_transcript_coverage\t$normal_transcript_depth\t$tumor_transcript_depth\n";
}
close OUT;

sub depth{
    my ($normal_bam,$tumor_bam,$chr,$start,$end,$n_minDep,$t_minDep)=@_;
    my @result=readpipe("samtools depth -a -r $chr:$start-$end $normal_bam $tumor_bam");
    my $normal_sum=0;
    my $normal_base=0;
    my $tumor_sum=0;
    my $tumor_base=0;
    my @re;
    foreach my $line(@result){
        chomp $line;
        my ($chr,$pos,$normal_num,$tumor_num)=split /\t/,$line;
        $normal_sum+=$normal_num;
        if($normal_num>=$n_minDep){
            $normal_base+=1;
        }
        $tumor_sum+=$tumor_num;
        if($tumor_num>=$t_minDep){
            $tumor_base+=1;
        }
    }
    $re[0]=sprintf "%.4f",$normal_sum;
    $re[1]=sprintf "%.4f",$normal_base;
    $re[2]=sprintf "%.4f",$tumor_sum;
    $re[3]=sprintf "%.4f",$tumor_base;
    return \@re;
}
