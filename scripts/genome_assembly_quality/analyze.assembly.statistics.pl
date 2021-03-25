#!/usr/bin/perl

use List::Util qw(first max maxstr min minstr reduce shuffle sum);
use strict;

##############################################################

sub readFasta{
    my $path=$_[0];
    my $reffasta=$_[1];
   
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
	    $id=$s[0];

	    # print "saw chromosome ".$id."\n";
	    
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

########################################################################

sub computeN50{
    my $scafseq=$_[0];

    my $totsize=0;
    my %hashsize;

    foreach my $chr (keys %{$scafseq}){
	my $size=length $scafseq->{$chr};

	if(exists $hashsize{$size}){
	    push(@{$hashsize{$size}}, $chr);
	} else{
	    $hashsize{$size}=[$chr];
	}

	$totsize+=$size;
    }

    my @uniquesizes=keys %hashsize;
    my @sortedsizes=sort {$a<=>$b} @uniquesizes;
    my @decreasingsizes=reverse @sortedsizes;

    my $halfsize=$totsize/2;

    my $n50="NA";

    my $currentsize=0;
    
    foreach my $size (@decreasingsizes){
	foreach my $chr (@{$hashsize{$size}}){
	    $currentsize+=$size;

	    if($currentsize>=$halfsize){
		$n50=$size;
		last;
	    }

	}

	if($currentsize>=$halfsize){
	    $n50=$size;
	    last;
	}
    }

    return $n50;
}

########################################################################

sub computeTotalSequenceLength{
    my $genome=$_[0];

    my $totsize=0;
    
    foreach my $chr (keys %{$genome}){
	my $size=length $genome->{$chr};
	$totsize+=$size;
    }

    return $totsize;
}

########################################################################

sub readBlatResults{
    my $pathin=$_[0];
    my $minpcid=$_[1];
    my $maxeval=$_[2];
    my $blatres=$_[3];

    open(my $input, $pathin);
    
    my $line=<$input>;
 
    while($line){
	chomp $line;
	my @s=split("\t", $line);

	my $scaffold=$s[0];
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
	    
	    if(exists $blatres->{$scaffold}){
		push(@{$blatres->{$scaffold}{"chr"}}, $chr);
		push(@{$blatres->{$scaffold}{"qstart"}}, $startquery);
		push(@{$blatres->{$scaffold}{"qend"}}, $endquery);
		push(@{$blatres->{$scaffold}{"dbstart"}}, $startdb);
		push(@{$blatres->{$scaffold}{"dbend"}}, $enddb);
		push(@{$blatres->{$scaffold}{"strand"}}, $strand);
	    } else{
		$blatres->{$scaffold}={"chr"=>[$chr], "qstart"=>[$startquery], "qend"=>[$endquery], "dbstart"=>[$startdb], "dbend"=>[$enddb], "strand"=>[$strand]};
	    }
	}
	
	$line=<$input>;
    }
    
    close($input);
}

########################################################################

sub orderBlatResultsScaffold{
    my $blatres=$_[0];
    my $sortedcoords=$_[1];
  
    foreach my $scaffold (keys %{$blatres}){
	$sortedcoords->{$scaffold}={"scaffoldstart"=>[], "scaffoldend"=>[], "chr"=>[], "chrstart"=>[], "chrend"=>[], "strand"=>[]};
	
	my %hashregions;

	my $nbhits=@{$blatres->{$scaffold}{"chr"}};

	for(my $i=0; $i<$nbhits; $i++){
	    my $chr=${$blatres->{$scaffold}{"chr"}}[$i];
	    my $qstart=${$blatres->{$scaffold}{"qstart"}}[$i];
	    my $qend=${$blatres->{$scaffold}{"qend"}}[$i];
	    my $dbstart=${$blatres->{$scaffold}{"dbstart"}}[$i];
	    my $dbend=${$blatres->{$scaffold}{"dbend"}}[$i];
	    my $strand=${$blatres->{$scaffold}{"strand"}}[$i];

	    if(exists $hashregions{$qstart}){
		push(@{$hashregions{$qstart}{"qend"}}, $qend);
		push(@{$hashregions{$qstart}{"dbstart"}}, $dbstart);
		push(@{$hashregions{$qstart}{"dbend"}}, $dbend);
		push(@{$hashregions{$qstart}{"chr"}}, $chr);
		push(@{$hashregions{$qstart}{"strand"}}, $strand);
	    } else{
		$hashregions{$qstart}={"qend"=>[$qend], "dbstart"=>[$dbstart], "dbend"=>[$dbend], "chr"=>[$chr], "strand"=>[$strand]};
	    }
	}

	my %orderedcoords;

	my @uniquestart=keys %hashregions;
	my @sortedstart=sort {$a<=>$b} @uniquestart;

	foreach my $start (@sortedstart){
	    my $nbpos=@{$hashregions{$start}{"qend"}};

	    for(my $i=0; $i<$nbpos; $i++){
		push(@{$sortedcoords->{$scaffold}{"scaffoldstart"}}, $start);
		push(@{$sortedcoords->{$scaffold}{"scaffoldend"}}, ${$hashregions{$start}{"qend"}}[$i]);
		push(@{$sortedcoords->{$scaffold}{"chr"}}, ${$hashregions{$start}{"chr"}}[$i]);
		push(@{$sortedcoords->{$scaffold}{"chrstart"}},  ${$hashregions{$start}{"dbstart"}}[$i]);
		push(@{$sortedcoords->{$scaffold}{"chrend"}}, ${$hashregions{$start}{"dbend"}}[$i]);
		push(@{$sortedcoords->{$scaffold}{"strand"}}, ${$hashregions{$start}{"strand"}}[$i]);
	    }
	}
    }
}

########################################################################

sub countBreakpointsScaffold{
    my $sortedblatcoords=$_[0];
    my $pathoutput=$_[1];
    
    my $nbbreakpoints=0;

    open(my $output, ">".$pathoutput);

    print $output "ScaffoldID\tScaffoldStart1\tScaffoldEnd1\tScaffoldStart2\tScaffoldEnd2\tChr1\tStart1\tEnd1\tStrand1\tChr2\tStart2\tEnd2\tStrand2\tBreakpointType\n";
    
    foreach my $scaffold (keys %{$sortedblatcoords}){
	my $nbhits=@{$sortedblatcoords->{$scaffold}{"scaffoldstart"}};

	if($nbhits>=2){
	    for(my $i=0; $i<($nbhits-1); $i++){
		my $thisscafstart=${$sortedblatcoords->{$scaffold}{"scaffoldstart"}}[$i];
		my $thisscafend=${$sortedblatcoords->{$scaffold}{"scaffoldend"}}[$i];

		my $nextscafstart=${$sortedblatcoords->{$scaffold}{"scaffoldstart"}}[$i+1];
		my $nextscafend=${$sortedblatcoords->{$scaffold}{"scaffoldend"}}[$i+1];
		
		my $thischr=${$sortedblatcoords->{$scaffold}{"chr"}}[$i];
		my $thisstart=${$sortedblatcoords->{$scaffold}{"chrstart"}}[$i];
		my $thisend=${$sortedblatcoords->{$scaffold}{"chrend"}}[$i];
		my $thisstrand=${$sortedblatcoords->{$scaffold}{"strand"}}[$i];

		my $nextchr=${$sortedblatcoords->{$scaffold}{"chr"}}[$i+1];
		my $nextstart=${$sortedblatcoords->{$scaffold}{"chrstart"}}[$i+1];
		my $nextend=${$sortedblatcoords->{$scaffold}{"chrend"}}[$i+1];
		my $nextstrand=${$sortedblatcoords->{$scaffold}{"strand"}}[$i+1];

		my $line=$scaffold."\t".$thisscafstart."\t".$thisscafend."\t".$nextscafstart."\t".$nextscafend."\t".$thischr."\t".$thisstart."\t".$thisend."\t".$thisstrand."\t".$nextchr."\t".$nextstart."\t".$nextend."\t".$nextstrand;
		
		if($thischr ne $nextchr){
		    print $output $line."\tDifferentChromosomes\n";
		    $nbbreakpoints++;
		} else{
		    if($thisstrand ne $nextstrand){
			print $output $line."\tDifferentStrands\n";
			$nbbreakpoints++;
		    } else{
			if(($thisstrand eq "+" && $nextstart<$thisend) || ($thisstrand eq "-" && $nextend>$thisstart)){
			    print $output $line."\tBadOrder\n";
			    $nbbreakpoints++;
			}
		    }
		}
	    }
	}
    }

    close($output);

    return $nbbreakpoints;
}

########################################################################

sub makeChrBlocks{
    my $blatres=$_[0]; ## original blat hits
    my $chrblocks=$_[1];
    
    my %hashregions;
        
    foreach my $scaffold (keys %{$blatres}){
	my $nbhits=@{$blatres->{$scaffold}{"chr"}};

	for(my $i=0; $i<$nbhits; $i++){
	    my $chr=${$blatres->{$scaffold}{"chr"}}[$i];
	    my $dbstart=${$blatres->{$scaffold}{"dbstart"}}[$i];
	    my $dbend=${$blatres->{$scaffold}{"dbend"}}[$i];
	
	    if(exists $hashregions{$chr}){
		if(exists $hashregions{$chr}{$dbstart}){
		    if($dbend>$hashregions{$chr}{$dbstart}){
			$hashregions{$chr}{$dbstart}=$dbend;
		    }
		} else{
		    $hashregions{$chr}{$dbstart}=$dbend;
		}
	    } else{
		$hashregions{$chr}={$dbstart=>$dbend};
	    }
	}
    }

    foreach my $chr (keys %hashregions){
	$chrblocks->{$chr}={"start"=>[], "end"=>[]};
	
	my @uniquestart=keys %{$hashregions{$chr}};
	my @sortedstart=sort {$a<=>$b} @uniquestart;

	my $currentstart=$sortedstart[0];
	my $currentend=$hashregions{$chr}{$currentstart};

	my $nbpos=@sortedstart;

	for(my $i=1; $i<$nbpos; $i++){
	    my $thisstart=$sortedstart[$i];
	    my $thisend=$hashregions{$chr}{$thisstart};

	    if($thisstart>=$currentstart && $thisstart<=($currentend+1)){
		if($thisend>$currentend){
		    $currentend=$thisend;
		}
	    } else{
		push(@{$chrblocks->{$chr}{"start"}}, $currentstart);
		push(@{$chrblocks->{$chr}{"end"}}, $currentend);

		$currentstart=$thisstart;
		$currentend=$thisend;
	    }
	}

	# don't forget last block
	push(@{$chrblocks->{$chr}{"start"}}, $currentstart);
	push(@{$chrblocks->{$chr}{"end"}}, $currentend);
	
    }
}

########################################################################

sub computeChromosomeCoverage{
    my $chrblocks=$_[0];
    my $chrcov=$_[1];
    
    foreach my $chr (keys %{$chrblocks}){
	my $nb=@{$chrblocks->{$chr}{"start"}};
	
	my $totsize=0;

	for(my $i=0; $i<$nb; $i++){
	    $totsize+=(${$chrblocks->{$chr}{"end"}}[$i]-${$chrblocks->{$chr}{"start"}}[$i]+1);
	}

	$chrcov->{$chr}=$totsize;
    }
}

########################################################################

sub printHelp{

    my $parnames=$_[0];
    my $parvalues=$_[1];
    
    print "\n";
    print "This script extracts genome assembly statistics, when a reference genome exists.\n";
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
$parameters{"pathAssembly"}="NA";
$parameters{"pathOutputStatistics"}="NA";

## "pathRepeatMasker", 

my @defaultpars=("pathAssembly",  "pathOutputStatistics");
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

print "Reading assembled sequences...\n";

my %assembly;
readFasta($parameters{"pathAssembly"}, \%assembly);

my $nbc=keys %assembly;
my $asslen=computeTotalSequenceLength(\%assembly);

print "There are ".$nbc." sequences and a total of ".$asslen." nucleotides.\n";

my $n50=computeN50(\%assembly);

print "N50: ".$n50."\n";
print "Done.\n";

##############################################################

print "Writing output...\n";
open(my $output, ">".$parameters{"pathOutputStatistics"});

print $output "AssemblySize:".$asslen."\n";
print $output "NbSequences:".$nbc."\n";
print $output "N50:".$n50."\n";

foreach my $chr (keys %assembly){
    my $seq=$assembly{$chr};
    my $size=length $seq;
    print $chr.":".$size."\n";
}
print "Done.\n";

##############################################################
##############################################################
