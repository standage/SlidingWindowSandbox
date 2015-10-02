#!/usr/bin/perl

use strict;
use warnings;
use Getopt::Long;


#### to whom it may concern: This code takes a genome along with associated annotation (gff3), and produces files with positional CpGoe (et al) values along gene frames for genes split into two classes: High CpG o/e and Low CpG o/e (those above and below the mean, respectively).  It produces a file with values for (n kb) into gene from start and stop of gene, as well as a given distance up and downstream of gene.  
#### NOTE: I wrote this code (or the base components, which remain largely unchanged) several years ago when first learning how to code.  Its messy, inefficient, and could be rewritten with a 1/4th of the variables.  If you want to rewrite it, be my guest...but don't say that I didn't warn you.

my $genomefasta = 0; #### fasta with genome (first "word" (consecutive group of non-whitespace char) must correspond to chr name in gff)
my $GFFfile0 = 0; #### gff file for annotation
my $output = "output.txt"; #### output name
my $proxDist = 1000; ### distance from gene start (and stop) to get coordinate data for (distance proximal in each dir to include along with gene frame plot)
my $geneLength = 2000; #### required length of gene to be included (eg 2000 includes genes 2kb and larger)
my $window=200;  #### length of window to use for calculating values
my $typer="gene"; ##either mRNA or gene
my $outputGeneSeqs=0; ##setting option outputs gene sequences
my $fromFront=0; ##setting option outputs front-aligned absolute gene coordinate data (eg all genes aligned at the start of gene frame w/ longest coordinate in file corresponding to last coordinate of longest gene) along with frong-back aligned
my $help =0;

GetOptions(
"out=s" => \$output, # 
"gen=s" => \$genomefasta, # 
"gff=s" => \$GFFfile0, # 
"gleng=i" => \$geneLength, # 
"prox=i" => \$proxDist, # 
"winleng=i" => \$window, # 
"type=s" => \$typer, #
"output_seqs" => \$outputGeneSeqs, #
"wfront" => \$fromFront, #
"help" => \$help); #
UsageAndExit() if $help;
UsageAndExit() unless($genomefasta && $GFFfile0);
my (%frame,%exCt,%chromTest,%chroms,%geneSeq,%frameTemp,%GFCpGs,%CtGFs,%CtGFsG,%GFGpCs)=();

############################################################################
#### Loads hashes up using input files from genome, paresed methylation ####
#### and GFF for element cacl										    ####
############################################################################

print "\nLoading genome into chromosomal hash...\n";
open(INPUT, $genomefasta) || die "$genomefasta file not found\n"; 
local $/ = "\n>"; 
for(<INPUT>){ 
	chomp;
	my $seq=$_;
	$_=~/^>*(.+)\n/;
	my $title=$1;
	$seq =~ s/^>*.+\n//;  #remove fasta header
	$seq=~ s/\n//g;
	$title=~/(^\S+)\s+.+$/;
	my $title2=$1;
	$chroms{$title2}=$seq;
}close(INPUT);


my $GFFfile = $GFFfile0."_tempcorrected.gff3";
system "perl convert_GFF.pl $GFFfile0 $typer"; #### this code


my ($frRf,$ctestRf)=readGFF($GFFfile);
%frame=%{$frRf};%chromTest=%{$ctestRf};

for my $chrom(keys %chroms){
	unless(exists $chromTest{$chrom}){
		delete $chroms{$chrom};  #limits chromosomes looped over to those with annotation information
	}
}


############################################################################
####gets sequence for introns and exons and then calcs CpGOE on them########
############################################################################
%frameTemp=%frame;
my $counting=0; 
print "\nGetting gene sequences...\n";
for my $chrom (sort keys %chroms){ ###loads up ORF seq
	for my $gene(keys %frameTemp){
		my @words=split("\t",$frameTemp{$gene});
		if ($words[0] eq $chrom){ ##if gene on scaffold in question
			my $chromSeq=$chroms{$chrom};
			my $length=($words[2]-$words[1])+2001;
			my $seqTemp=substr($chromSeq,($words[1]-1000),$length);
			if($words[3] eq "-"){
				my $tempor=REVCOM($seqTemp);  ### get reverse compliment seq if negative strand
				$geneSeq{$gene}=$tempor;
				$counting++;
			}else{
				$geneSeq{$gene}=$seqTemp;
				$counting++;
			}
			delete $frameTemp{$gene};
		}
	}
}
%frameTemp=();
print "\t$counting gene seqs read in...\n";
if($outputGeneSeqs == 1){
	print "\toutputing gene sequences...\n";
	open(OUT, ">$output".".frameSeq.fasta") or die "Cannot open the output file\n";
	for my $key (keys %geneSeq){
		my $seqr=$geneSeq{$key};
		unless($seqr=~m/N+/){
			$seqr=~s/n//g;
			$seqr=~tr/acgt/ACGT/;
			print OUT ">".$key."\n";
			print OUT $seqr."\n";
		}
	}close(OUT);
}else{ print "\tskipping output of gene frame sequences...\n"; }



print "\nCalculating CpGOE on elements...\n";
my ($CpGref,$GpCref,$GCref,$ATref) = GetCpG(\%geneSeq,$window,$geneLength);
my %CpGfr = %{$CpGref};
my %GpCfr = %{$GpCref};
my %ATs = %{$ATref};
my %GCs = %{$GCref};

my ($CpGref2,$GpCref2) = SplitCpG(\%geneSeq,$window);
my %CpGfrLow = %{$CpGref2};
my %CpGfrHigh = %{$GpCref2};



print "\nGetting front-back alignment coordinate average...\n";
for my $key (keys %CpGfr){
	my @title=split(";",$key);
	my @words=split("\t",$frame{$title[0]});
	my $length=($words[2]-$words[1])+1;
	my $class="high";
	$class="low" if exists $CpGfrLow{$title[0]};
	if($length >=$geneLength){
		my $step=($window/10);
		my $limit=((($geneLength/2)+$proxDist)/$step);
		if($title[1] <= $limit){
			if(exists $GFCpGs{$title[1]}{$class}){
				$GFCpGs{$title[1]}{$class}+=$CpGfr{$key};
				$GFGpCs{$title[1]}{$class}+=$GpCfr{$key};
				$CtGFs{$title[1]}{$class}++;
				$CtGFsG{$title[1]}{$class}++;
			}else{
				$GFCpGs{$title[1]}{$class}=$CpGfr{$key};
				$GFGpCs{$title[1]}{$class}=$GpCfr{$key};
				$CtGFs{$title[1]}{$class}=1;
				$CtGFsG{$title[1]}{$class}=1;
			}
		}elsif((($title[2]-$title[1])+1) <= $limit){
			my $posit=((($geneLength+($proxDist*2)))/$step)-($title[2]-$title[1]);
			if(exists $GFCpGs{$posit}{$class}){
				$GFCpGs{$posit}{$class}+=$CpGfr{$key};
				$GFGpCs{$posit}{$class}+=$GpCfr{$key};
				$CtGFs{$posit}{$class}++;
				$CtGFsG{$posit}{$class}++;
			}else{
				$GFCpGs{$posit}{$class}=$CpGfr{$key};
				$GFGpCs{$posit}{$class}=$GpCfr{$key};
				$CtGFs{$posit}{$class}=1;
				$CtGFsG{$posit}{$class}=1;
			}
		}
	}
}

print "\tprinting output...\n";

my @tmpclss=("low","high");

open(OUT, ">$output".".front_back_aligned.txt") or die "Cannot open the output file\n";
print OUT "type\tposition";
for my $class(@tmpclss){ print OUT "\tCpGoe_".$class."\tGpCoe_".$class."\tNumSites_".$class; }
print OUT "\n";
for my $key (sort {$a <=> $b} keys %GFCpGs){
	my $nummy=($key*($window/10));
	my $state="FRAME";
	print OUT "$state\t$nummy";
	for my $class(@tmpclss){
		my $trueCpG="NA";my $trueGpC="NA";
		if(exists $CtGFs{$key}{$class}){
			$trueCpG=$GFCpGs{$key}{$class}/$CtGFs{$key}{$class} if $CtGFs{$key}{$class} > 0;
			$trueGpC=$GFGpCs{$key}{$class}/$CtGFsG{$key}{$class} if $CtGFsG{$key}{$class} > 0;
		}
		print OUT "\t$trueCpG\t$trueGpC\t$CtGFs{$key}{$class}"
	}
	print OUT "\n";
}close(OUT);




if($fromFront == 1){
	print "\nGetting front-aligned output for associated genes...\n";
	my (%GFAT,%GFGC,%GFCpG,%CtGF,%GFGpC,%CtGFg,%allCoor)=();  ###this is a stupid number of variables for such a simple computation.  i wrote this 2 years ago when i was first learning perl...so sue me.
	for my $key (keys %CpGfr){
		my @title=split(";",$key);
		my @words=split("\t",$frame{$title[0]});
		my $length=($words[2]-$words[1])+1;
		my $class="high";
		$class="low" if exists $CpGfrLow{$title[0]};
		if($length >=$geneLength){
			$allCoor{$title[1]}=1 unless exists $allCoor{$title[1]};
			if(exists $GFCpG{$title[1]}{$class}){
				$GFCpG{$title[1]}{$class}+=$CpGfr{$key};
				$CtGF{$title[1]}{$class}++;
				$GFGpC{$title[1]}{$class}+=$GpCfr{$key};
				$CtGFg{$title[1]}{$class}++;
				$GFGC{$title[1]}{$class}+=$GCs{$key};
				$GFAT{$title[1]}{$class}+=$ATs{$key};
				
			}else{
				$GFCpG{$title[1]}{$class}=$CpGfr{$key};
				$CtGF{$title[1]}{$class}=1;
				$GFGpC{$title[1]}{$class}=$GpCfr{$key};
				$CtGFg{$title[1]}{$class}=1;
				$GFGC{$title[1]}{$class}=$GCs{$key};
				$GFAT{$title[1]}{$class}=$ATs{$key};
			}
		}
	}
	print "\tprinting output...\n";
	open(OUT, ">".$output.".front_aligned.txt") or die "Cannot open the output file\n";
	print OUT "type\tposition";
	for my $class(@tmpclss){
		print OUT "\tCpGoe_".$class."\tGpCoe_".$class."\tGCcontent_".$class."\tATcontent_".$class;
	}print OUT "\n";
	for my $key (sort {$a <=> $b} keys %allCoor){
		my $nummy=($key*($window/10))+($window/2);
		my $state="FRAME";
		print OUT "$state\t$nummy";
		for my $class(@tmpclss){
			
			my $trueCpG="NA"; my $trueGpC="NA"; my $trueGC="NA"; my $trueAT="NA";
			if(exists $CtGF{$key}{$class} && $CtGF{$key}{$class} > 0){
				$trueCpG=$GFCpG{$key}{$class}/$CtGF{$key}{$class};
				$trueGC=$GFGC{$key}{$class}/$CtGF{$key}{$class};
				$trueAT=$GFAT{$key}{$class}/$CtGF{$key}{$class};
			}
			$trueGpC=$GFGpC{$key}{$class}/$CtGFg{$key}{$class} if(exists $CtGFg{$key}{$class} && $CtGFg{$key}{$class} > 0);
			print OUT "\t$trueCpG\t$trueGpC\t$trueGC\t$trueAT";
		}print OUT "\n";
	}close(OUT);
}else{ print "skipping front-alignment averaging...\n"; }




sleep 1;

print "\nDONE\n\n";



###################################
############# SUBS ################
###################################

####													   

sub GetCpG{ #### takes hash of sequences, breaks them up into windows of length n (cmd ln) with a window every n/10 bp, then finds CpG o/e for them
	my ($seqref, $cutoff1,$lengther) = @_; 
	my %seqhash=%{$seqref};
	my (%CpGseq,%GpCseq,%ATcon,%GCcon,%geneCter,%SeqSplit)=();
	my $countTime=0;
	for my $key(sort keys %seqhash){
		my $noline=$seqhash{$key};
		
		$noline=~s/\n//g;
		my $templength=length($noline);
		if(($templength-($proxDist*2)) >=$lengther){
			my $counter=0;
			for (my $i = 0; $i < $templength ; $i+=($cutoff1/10)) {
				my $Seqer=substr($noline,$i,$cutoff1);
				$counter++;
				$SeqSplit{$key}{$counter}=$Seqer;
				$geneCter{$key}++;
			}
			$countTime++;
		}
	}
	print "\t$countTime genes pass length cutoff.\n";
	my $gnct=0;
	for my $gn (keys %SeqSplit){
		for my $ctr (keys %{$SeqSplit{$gn}}){
			my $key = $gn.";".$ctr;
			my $noline=$SeqSplit{$gn}{$ctr};
			my @letters=split(/ */,$noline);
			my $Ccount=0; my $Gcount=0; my $CGcount=0; my $GCcount=0;
			my $lastchar="X";
			my $Cfreq=0; my $Acount=0; my $Tcount=0;
			my $Gfreq=0; my $GCcontent=0; my $ATcontent=0;
			my $CGfreq=0; my $GCfreq=0; my $CpGOE=0; my $GpCOE=0;
			for(@letters){
				my $dinuc=$lastchar.$_;
				if($_ =~ m/C/i){
					$Ccount++;
				}elsif($_ =~ m/G/i){
					$Gcount++;
				}elsif($_ =~ m/A/i){
					$Acount++;
				}elsif($_ =~ m/T/i){
					$Tcount++;
				}if($dinuc =~ /CG/i){
					$CGcount++;
				}elsif($dinuc =~ /GC/i){
					$GCcount++;
				}
				$lastchar=$_;
			}
			$noline=~s/[Xx]//g;
			$noline=~s/[Nn]//g;
			my $length=length($noline);
			if ($length>10){
				$Cfreq=$Ccount/$length;
				$Gfreq=$Gcount/$length;
				$ATcontent=($Acount+$Tcount)/($Acount+$Tcount+$Ccount+$Gcount);
				$GCcontent=($Gcount+$Ccount)/($Acount+$Tcount+$Ccount+$Gcount);
				$CGfreq=$CGcount/($length);
				$GCfreq=$GCcount/($length);
				if(($Cfreq>0)&&($Gfreq>0)){
					$CpGOE=($CGfreq/($Cfreq*$Gfreq));
					$GpCOE=$GCfreq/($Cfreq*$Gfreq);
				}if($length>=($cutoff1/3)){
					my $maxer=$geneCter{$gn};
					$ATcon{$key.";".$maxer}=$ATcontent;
					$GCcon{$key.";".$maxer}=$GCcontent;
					if($CGcount > 0 &&($Cfreq>0 && $Gfreq>0)){
						$CpGseq{$key.";".$maxer}=$CpGOE;
						$GpCseq{$key.";".$maxer}=$GpCOE;
					}
				}
			}
		}
		$gnct++;
		if(($gnct % 250) == 0){
			print "\t\t".$gnct." (".(sprintf("%.2f",(($gnct/$countTime)*100)))."\%) genes processed...\n"
		}
	}
	print "\t\tdone\n";
	return (\%CpGseq,\%GpCseq,\%GCcon,\%ATcon);	
}

#### 

sub MAX{ ### returns maximum of nums
	my $current = $_[0];
	for my $num(@_){
		$current=$num if $num > $current;
	}return $current;
}
sub MIN{ ### returns min of nums
	my $current = $_[0];
	for my $num(@_){
		$current=$num if $num < $current;
	}return $current;
}

#### 

sub REVCOM {  #### gets reverse compliment of given sequence
  my ($DNA) = @_;
  my $revcomp = reverse($DNA);
  $revcomp =~ tr/ACGTacgt/TGCAtgca/;
  return $revcomp; 
}

#### 

sub SplitCpG{  #### takes hash of sequences, finds CpG o/e
	my ($seqref, $cutoff1) = @_; 
	my %seqhash=%{$seqref};
	my (%CpGseq,%CpGlow,%CpGhigh,%SeqSplit)=();
	print "\tfinding genes below and above mean CpGoe.\n";
	for my $key (keys %seqhash){
		my $noline=$seqhash{$key};
		my @letters=split(/ */,$seqhash{$key});
		my $Ccount=0; my $Gcount=0; my $CGcount=0; my $lastchar="X";
		my $Cfreq=0; my $Gfreq=0; my $CGfreq=0; my $CpGOE=0;
		for(@letters){
			my $dinuc=$lastchar.$_;
			if($_ =~ m/C/i){
				$Ccount++;
			}elsif($_ =~ m/G/i){
				$Gcount++;
			}if($dinuc =~ /CG/i){
				$CGcount++;
			}
			$lastchar=$_;
		}
		$noline=~s/N//g;
		my $length=length($noline);
		if ($length>10){
			$Cfreq=$Ccount/$length;
			$Gfreq=$Gcount/$length;
			$CGfreq=$CGcount/($length);
			if(($Cfreq>0)&&($Gfreq>0)){
				$CpGOE=($CGfreq/($Cfreq*$Gfreq));
			}if($length>=($cutoff1/3) ){
				if($CGcount > 0 &&($Cfreq>0 && $Gfreq>0)){
					$CpGseq{$key}=$CpGOE;						
				}
			}
		}
	}
	my $CpGnum=0;
	my $CpGct=0;
	for my $key (keys %CpGseq){
		$CpGnum+=$CpGseq{$key};
		$CpGct++;
	}
	my $CpGav=$CpGnum/$CpGct;
	for my $key (keys %CpGseq){
		if($CpGseq{$key} < $CpGav){
			$CpGlow{$key}=1;
		}else{
			$CpGhigh{$key}=1;
		}
	}
	return (\%CpGlow,\%CpGhigh);
}

#### 

sub readGFF{ #### reads gff file prepared by code "convert_GFF.pl" only looking at "gene" or "mRNA" entries for gene frame
	my ($geneGFF)=@_;
	my (%framer,%chromTestr,%gffhash)=();
	local $/ = "\n";
	print "\nReading in GFF file for gene coordinates...\n";
	open(INPUT, $geneGFF) || die "$geneGFF file not found\n"; ### load gff into hash in order to sort by position
	my @gffar=<INPUT>; close(INPUT);
	for my $ln(@gffar){
		$ln=~s/\s*$//g;
		my @words=split("\t", $ln);
		if($words[2] =~ /(mrna|gene)/i){
			$gffhash{$words[0]}{$words[3]}{$words[4]}=$ln;
			$chromTestr{$words[0]} = 1 unless exists $chromTestr{$words[0]};
		}
	}
	####gets coordinates from gff hash
	for my $chr(keys %gffhash){
		for my $c1(sort {$a <=> $b} keys %{$gffhash{$chr}}){
			for my $c2(sort {$a <=> $b} keys %{$gffhash{$chr}{$c1}}){
				my @words=split("\t", $gffhash{$chr}{$c1}{$c2});
				if($words[2] =~ /(mrna|gene)/i){
					$words[8]=~/ID=([^\;]+)\;/m;
					my $gene=$1;
					if(exists $framer{$gene}){
						my @words2=split("\t",$framer{$gene});
						my $start=$words2[1];my $stop=$words2[2];
						$start=$words[3] if $words[3] < $start;
						$stop=$words[4] if $words[4] > $stop;
						$framer{$gene}="$words[0]\t$start\t$stop\t$words2[3]\t$words[5]";
					}else{
						my $start=$words[3];my $stop=$words[4];
						$framer{$gene}="$words[0]\t$words[3]\t$words[4]\t$words[6]\t$words[5]";
					}
				}
			}
		}
	}
	return (\%framer,\%chromTestr);
}

####

sub UsageAndExit {
    print "\noptions:\n";
    print "\t--out:\tname used in output file(s)\n";
    print "\t--gen:\tgenome fasta file[REQUIRED]\n";
    print "\t--gff:\tannotation file in gff3 format (i think)[REQUIRED]\n";
    print "\t--gleng:\tlength of gene required to be lincluded in averaging/plotting\n";
    print "\t--winleng:\tlength of windows (step == 10th the lenght of window) used in sliding window CpG o/e\n";
    print "\t--type:\teither \"mRNA\" or \"gene\" related to which gff feature to use\n";
    print "\t--output_seqs:\tsetting this option outputs gene sequences (+prox regions) into file\n";
    print "\t--wfront:\tsetting this option results in the output of front-aligned averages as well as front-back\n";
    print "\t--prox:\tdistance (bp) from start and stop of gene to include in averaged plot output\n";
    print "\t--help:\tprint this and exit\n\n";
    print "NOTE: gff formats vary.  the code \"convert_GFF.pl\" is used by this code for parsing the genome into a standard format this code can use [REQUIRED in same directory].  The simplest way to accommodate a different annotation format is to modify that code, instead of dealing with the poorly-written parser within this code\n\n";
    print "\n";
    exit;
}
