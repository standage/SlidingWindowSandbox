#!/usr/bin/perl -w

$USAGE="\nUsage: $0 FASTAfile seqtype\n

** This script derives CpG data for a FASTA-formatted
   input file.\n";


if(scalar(@ARGV)<2){
  print STDERR "$USAGE\n";
  exit;
}

my $infname = $ARGV[0];
my $type    = $ARGV[1];
my $verbose = 0;
my $count   = 0;
my $sequence = "";

open(INFILE,"$infname");
if(! defined($_=<INFILE>)){ die $USAGE;}
print "SeqID\tSeqNbr\tSeqType\tLength\tCcount\tGcount\tCGcount\tGCcount\tCpcnt\tGpcnt\tCpGOE\tGpCOE\n";

if(! /^>/) {
  while (<INFILE>) {
    if (/^>/) {last;}
  }
}

my $header = $_;

while (<INFILE>) {
  if (/^>/) {
    &doit($header,$sequence,$type);
    $sequence = "";
    $header = $_;
    next;
  }
  chomp($sequence .= $_);
  if(/^\n$/) {next;}
}
&doit($header,$sequence,$type);



sub doit() {
  my $head = $_[0];
  my $seq  = $_[1];
  my $typ  = $_[2];
  $count++;

  if ($verbose) {&printheader($head); &printseq($seq);}
  (my $seqn) = ($head =~ m/^>(\S+)/);
  my($len,$Cct,$Gct,$CGct,$GCct) = GetCpG($seq);
  my($Cfreq,$Gfreq,$CpGOE,$GpCOE) = oeCpG($len,$Cct,$Gct,$CGct,$GCct);
  printf "$seqn\t$count\t$typ\t$len\t$Cct\t$Gct\t$CGct\t$GCct\t%5.1f\t%5.1f\t%.2f\t%.2f\n", 100*$Cfreq, 100*$Gfreq, $CpGOE, $GpCOE;
}


sub printheader() {
  my $header = $_[0];
  print $header;
}

sub printseq() {
  my $seq = $_[0];
  print $seq, "\n";
}

sub GetCpG {
	my ($seqer) = @_;

	my $noline=$seqer;
	$noline=~s/\n//g;
	my $fulllength=length($noline);
	my $elNum="NA";
	my @letters=split(/ */,$noline);
	my $Ccount=0;
	my $Gcount=0;
	my $CGcount=0;
	my $GCcount=0;
	my $lastchar="X";
	for(@letters){
		my $dinuc=$lastchar.$_;
		if(($_ eq"C")||($_ eq "c")){
			$Ccount++;
		}elsif(($_ eq"G")||($_ eq "g")){
			$Gcount++;
		}if(($dinuc eq"CG")||($dinuc eq "cg")){
			$CGcount++;
		}elsif(($dinuc eq"GC")||($dinuc eq "gc")){
			$GCcount++;
		}
		$lastchar=$_;
	}
	$noline=~s/[XNxn]//g;
	my $length=length($noline);
	return($length,$Ccount,$Gcount,$CGcount,$GCcount);
}


sub oeCpG {
	my $length = $_[0];
	my $Cct    = $_[1];
	my $Gct    = $_[2];
	my $CGct   = $_[3];
	my $GCct   = $_[4];

	my $Cfreq = 0;
	my $Gfreq = 0;
	my $CpGOE = -1;
	my $GpCOE = -1;

	if ($length > 0){
		$Cfreq=$Cct/$length;
		$Gfreq=$Gct/$length;
		if ($length > 1){
			$CGfreq=$CGct/($length-1);
			$GCfreq=$GCct/($length-1);
			if (($Cfreq*$Gfreq) > 0){
				$CpGOE =$CGfreq/($Cfreq*$Gfreq);
				$GpCOE =$GCfreq/($Cfreq*$Gfreq);
			}
		}
	}
	return($Cfreq,$Gfreq,$CpGOE,$GpCOE);
}
