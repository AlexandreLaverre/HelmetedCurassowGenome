#!/usr/bin/perl

use List::Util qw(first max maxstr min minstr reduce shuffle sum);
use strict;

##############################################################

sub readFasta{
    my $path=$_[0];
    my $reffasta=$_[1];
    my $geneids=$_[2];
   
    my @s=split("\\.",$path);
    my $ext=$s[-1];

    my $input;

    if($ext eq "gz"){
	open($input,"zcat $path |");
    }
    else{
	open($input, $path);
    }
    
    my $line=<$input>;

    while($line){
	my $b=substr $line,0,1;
	
	if($b eq ">"){
	    chomp $line;
	    my $id=substr $line,1;

	    my @s=split(" ",$id);
	    my @t=split("\\.", $s[0]);
	    $id=$t[0];

	    my $gene="NA";

	    for(my $k=1; $k<@s; $k++){
		my @t=split(":", $s[$k]);
		if($t[0] eq "gene"){
		    $gene=$t[1];
		}
	    }
	    
	    $geneids->{$id}=$gene;
	    
	    $reffasta->{$id}="";
	   

	    $line=<$input>;
	    $b=substr $line,0,1;
	    
	    while($line && !($b eq ">")){
		chomp $line;
		$reffasta->{$id}.=$line;
		$line=<$input>;
		$b=substr $line,0,1;
	    }
	}
    }

    close($input);
}

##############################################################

sub readCoordinates{
    my $pathin=$_[0];
    my $coords=$_[1];

    open(my $input, $pathin);

    my $line=<$input>;
    chomp $line;
    my %header;
    my @s=split("\t", $line);

    for(my $i=0; $i<@s; $i++){
	$header{$s[$i]}=$i;
    }

    $line=<$input>;

    while($line){
	chomp $line;

	my @s=split("\t", $line);
	my $protid=$s[$header{"Protein stable ID"}];
	my $geneid=$s[$header{"Gene stable ID"}];
	my $chr=$s[$header{"Chromosome/scaffold name"}];
	my $start=$s[$header{"Gene start (bp)"}];
	my $end=$s[$header{"Gene end (bp)"}];
	my $strand=$s[$header{"Strand"}];

	$coords->{$protid}={"chr"=>$chr, "start"=>$start, "end"=>$end, "strand"=>$strand, "gene"=>$geneid};
	
	$line=<$input>;
    }
    
    close($input);
}

##############################################################

sub readTBlastNResults{
    my $pathin=$_[0];
    my $minpcid=$_[1];
    my $maxeval=$_[2];
    my $tblastnres=$_[3];

    open(my $input, $pathin);
    
    my $line=<$input>;
 
    while($line){
	chomp $line;
	my @s=split("\t", $line);
	my @t=split("\\.",$s[0]);

	my $protein=$t[0];
	my $chr=$s[1];
	my $pcid=$s[2]+0.0;
	my $evalue=$s[10]+0.0;

	if($pcid>=$minpcid && $evalue<=$maxeval){
	    my $startquery=$s[6]+0;
	    my $endquery=$s[7]+0;
	    my $startdb=$s[8]+0;
	    my $enddb=$s[9]+0;

	    my $strand="+";

	    if($startdb>$enddb){
		$strand="-";
		$startdb=$s[9]+0;
		$enddb=$s[8]+0;
	    }
	    
	    if(exists $tblastnres->{$protein}){
		push(@{$tblastnres->{$protein}{"chr"}}, $chr);
		push(@{$tblastnres->{$protein}{"qstart"}}, $startquery);
		push(@{$tblastnres->{$protein}{"qend"}}, $endquery);
		push(@{$tblastnres->{$protein}{"dbstart"}}, $startdb);
		push(@{$tblastnres->{$protein}{"dbend"}}, $enddb);
		push(@{$tblastnres->{$protein}{"strand"}}, $strand);
	    } else{
		$tblastnres->{$protein}={"chr"=>[$chr], "qstart"=>[$startquery], "qend"=>[$endquery], "dbstart"=>[$startdb], "dbend"=>[$enddb], "strand"=>[$strand]};
	    }
	}
	
	$line=<$input>;
    }
    
    close($input);
}

##############################################################

sub computeTBlastNStatistics{
    my $tblastn=$_[0];
    my $stats=$_[1];

    foreach my $prot (keys %{$tblastn}){
	my %hashhits;

	my $nbh=@{$tblastn->{$prot}{"chr"}};

	for(my $i=0; $i<$nbh; $i++){
	    my $chr=${$tblastn->{$prot}{"chr"}}[$i];
	    my $strand=${$tblastn->{$prot}{"strand"}}[$i];

	    my $id=$chr.":".$strand;
	    
	    my $qstart=${$tblastn->{$prot}{"qstart"}}[$i];
	    my $qend=${$tblastn->{$prot}{"qend"}}[$i];
	    my $dbstart=${$tblastn->{$prot}{"dbstart"}}[$i];
	    my $dbend=${$tblastn->{$prot}{"dbend"}}[$i];

	    if(exists $hashhits{$id}){
		if(exists $hashhits{$id}{$qstart}){
		    push(@{$hashhits{$id}{$qstart}{"qend"}}, $qend);
		    push(@{$hashhits{$id}{$qstart}{"dbstart"}}, $dbstart);
		    push(@{$hashhits{$id}{$qstart}{"dbend"}}, $dbend);
		} else{
		    $hashhits{$id}{$qstart}={"qend"=>[$qend], "dbstart"=>[$dbstart], "dbend"=>[$dbend]};
		}
	    } else{
		$hashhits{$id}={$qstart=>{"qend"=>[$qend], "dbstart"=>[$dbstart], "dbend"=>[$dbend]}};
	    }
	}

	foreach my $id (keys %hashhits){
	    my @uniquestart=keys %{$hashhits{$id}};
	    my @sortedstart=sort {$a<=>$b} @uniquestart;

	    my @startblocks;
	    my @endblocks;

	    my $currentstart=$sortedstart[0];
	    my $currentend=max @{$hashhits{$id}{$currentstart}{"qend"}};

	    for(my $i=1; $i<@sortedstart; $i++){
		my $thisstart=$sortedstart[$i];
		my $thisend=max @{$hashhits{$id}{$thisstart}{"qend"}};

		if($thisstart>=$currentstart && $thisstart<=($currentend+1)){
		    if($thisend>$currentend){
			$currentend=$thisend;
		    }
		} else{
		    push(@startblocks, $currentstart);
		    push(@endblocks, $currentend);

		    $currentstart=$thisstart;
		    $currentend=$thisend;
		}
	    }
	    
	    push(@startblocks, $currentstart);
	    push(@endblocks, $currentend);
	    
	    
	    my $totlen=0;

	    for(my $i=0; $i<@startblocks; $i++){
		$totlen+=($endblocks[$i]-$startblocks[$i]+1);
	    }

	    if(exists $stats->{$prot}){
		$stats->{$prot}{$id}={"totlength"=>$totlen};
	    } else{
		$stats->{$prot}={$id=>{"totlength"=>$totlen}};
	    }
	}
    }
}

##############################################################

sub extractBestHits{
    my $blaststats=$_[0];
    my $besthits=$_[1];

    foreach my $prot (keys %{$blaststats}){
	my %hashlen;

	foreach my $id (keys %{$blaststats->{$prot}}){
	    my $len=$blaststats->{$prot}{$id};

	    if(exists $hashlen{$len}){
		push(@{$hashlen{$len}}, $id);
	    } else{
		$hashlen{$len}=[$id];
	    }
	}

	my $maxlen=max (keys %hashlen);
	my @bh=@{$hashlen{$maxlen}};

	$besthits->{$prot}=$bh[0];
    }
}

##############################################################

sub printHelp{

    my $parnames=$_[0];
    my $parvalues=$_[1];
    
    print "\n";
    print "This script checks the chromosome correspondence with a reference species, based on tblastn results.\n";
    print "\n";
    
    print "Options:\n";
    
    foreach my $par (@{$parnames}){
	print "--".$par."  [  default value: ".$parvalues->{$par}."  ]\n";
    }
    print "\n";
}

##############################################################
##############################################################

## parameters 

my %parameters;
$parameters{"pathGeneCoordinates"}="NA";
$parameters{"pathProteins"}="NA";
$parameters{"pathTBlastNResults"}="NA";
$parameters{"minPCIdentity"}="NA";
$parameters{"maxEValue"}="NA";
$parameters{"pathOutput"}="NA";

## "pathRepeatMasker", 

my @defaultpars=("pathGeneCoordinates", "pathProteins", "pathTBlastNResults", "minPCIdentity", "maxEValue", "pathOutput");
my %defaultvalues;

foreach my $par (keys %parameters){
    $defaultvalues{$par}=$parameters{$par};
}

## update arguments

my $nbargs=@ARGV;

for(my $i=0;$i<$nbargs; $i++){
    my $arg=$ARGV[$i];
    $arg=substr $arg,2;
    
    my @s=split("=",$arg);
    my $parname=$s[0];
    my $parval=$s[1];
	
    if(exists $parameters{$parname}){
	$parameters{$parname}=$parval;
    }
    else{
	print "Error: parameter ".$parname." was not recognized!!!\n";
	printHelp(\@defaultpars, \%defaultvalues);
	exit(1);
    }
}


## show parameters

print "\n";

print "Running program with the following parameters:\n";

foreach my $par (@defaultpars){
    print "--".$par."=".$parameters{$par}."\n";
}

print "\n";

##############################################################
##############################################################

print "Reading gene coordinates...\n";
my %coords;
readCoordinates($parameters{"pathGeneCoordinates"}, \%coords);
print "Done.\n";

##############################################################

print "Reading protein sequences...\n";

my %proteins;
my %geneids;

if(-e $parameters{"pathProteins"}){
    readFasta($parameters{"pathProteins"}, \%proteins, \%geneids);
} else{
    print "cannot find protein sequence file\n";
    exit(1);
}
my $nbprot=keys %proteins;

print "There are ".$nbprot." proteins in total.\n";

print "Done.\n";

##############################################################


print "Reading and analyzing tblastn results...\n";

my $minpcid=$parameters{"minPCIdentity"}+0.0;
my $maxeval=$parameters{"maxEValue"}+0.0;

print "minimum % identity: " .$minpcid."\n";
print "max e-value: ".$maxeval."\n";

my %tblastnres;

readTBlastNResults($parameters{"pathTBlastNResults"}, $minpcid, $maxeval, \%tblastnres);

my $nbhits=keys %tblastnres;

print "Found ".$nbhits." proteins with significant hits.\n";
 
print "Done.\n";

##############################################################

print "Extracting alignment statistics...\n";

my %stats;

computeTBlastNStatistics(\%tblastnres, \%stats);

print "Done.\n";

##############################################################

print "Extracting best hits...\n";

my %besthits;
extractBestHits(\%stats, \%besthits);

print "Done.\n";

##############################################################

print "Writing output...\n";

open(my $output, ">".$parameters{"pathOutput"});

my %min50prot=0;
my %min50genes=0;

print $output "ProteinID\tGeneID\tReferenceChr\tReferenceStart\tReferenceEnd\tReferenceStrand\tContig\tStrand\tTotalLength\tAlignedLength\n";

foreach my $prot (keys %besthits){
    my $gene=$coords{$prot}{"gene"};
    my $refchr=$coords{$prot}{"chr"};
    my $refstart=$coords{$prot}{"start"};
    my $refend=$coords{$prot}{"end"};
    my $refstrand=$coords{$prot}{"strand"};

    my $bh=$besthits{$prot};
    my $alnlen=$stats{$prot}{$bh}{"totlength"};
    
    my @t=split(":",$bh);
    my $tgchr=$t[0];
    my $tgstrand=$t[1];

    my $totlen=length $proteins{$prot};
    
    print $output $prot."\t".$gene."\t".$refchr."\t".$refstart."\t".$refend."\t".$refstrand."\t".$tgchr."\t".$tgstrand."\t".$totlen."\t".$alnlen."\n";

} 

close($output);

print "Done.\n";

##############################################################
##############################################################
