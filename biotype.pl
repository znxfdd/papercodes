use strict;
#use warnings;

my $gtfFile="human.gtf";
my $expFile="symbol.txt";

my %hash=();
my %proteinHash=();

open(RF,"$gtfFile") or die $!;
while(my $line=<RF>)
{
        chomp($line);
        if($line=~/gene_id \"(.+?)\"\;.+gene_name "(.+?)"\;.+gene_biotype \"(.+?)\"\;/)
        {
                my $ensembl=$1;
                my $symbol=$2;
                my $biotype=$3; 
                if($biotype eq "protein_coding"){
                     $proteinHash{$symbol}=1;
                     $hash{$symbol}="$symbol|$biotype";
                }
                elsif($biotype=~/3prime_overlapping_ncrna|lncRNA|ambiguous_orf|ncrna_host|non_coding|processed_transcript|retained_intron|antisense|sense_overlapping|sense_intronic|bidirectional_promoter_lncrna|lincRNA/){
                     unless(exists $proteinHash{$symbol}){
                           $hash{$symbol}="$symbol|lncRNA";
                     }
                }
        }
}
close(RF);

my @samp1e=(localtime(time));
open(RF,"$expFile") or die $!;
open(PROTEIN,">mRNA.txt") or die $!;
open(LNCRNA,">lncRNA.txt") or die $!;
while(my $line=<RF>)
{
                                                                                                              

        if($.==1){
            print PROTEIN $line;
            print LNCRNA $line;
            next;
        }
        chomp($line);
        my @arr=split(/\t/,$line);
        if(exists $hash{$arr[0]}){
              if($hash{$arr[0]}=~/lncRNA/){
                  print LNCRNA join("\t",@arr) . "\n";
              }
              else{
                  print PROTEIN join("\t",@arr) . "\n";
              }
        }
}
close(LNCRNA);
close(PROTEIN); 
close(RF);