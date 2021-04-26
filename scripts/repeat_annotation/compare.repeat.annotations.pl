#!/usr/bin/perl

use List::Util qw(first max maxstr min minstr reduce shuffle sum);
use strict;

##############################################################

sub readRepeatMasker{
    my $pathin=$_[0];
    my $repeatmasker=$_[1];

    my %hashrepeats;
    my @s=split("\\.", $pathin);
    my $ext=$s[-1];
    my $input;

    if($ext eq "gz"){
	open($input,"zcat $pathin |");
    } else{
	open($input,$pathin);
    }
    
    my $line=<$input>; ## skip1
    $line=<$input>; ## skip2
    $line=<$input>; ## skip3
   
    $line=<$input>;

    while($line){
	chomp $line;
	my @s=split(" ",$line);
	
	my $chr=$s[4];
	my $start=$s[5]+0; 
	my $end=$s[6]+0;
	my $strand=$s[8];
	my $classfam=$s[10];
		
	if(exists $hashrepeats{$chr}){
	    if(exists $hashrepeats{$chr}{$start}){
		push(@{$hashrepeats{$chr}{$start}{"end"}},$end);
		push(@{$hashrepeats{$chr}{$start}{"classfam"}},$classfam);
		push(@{$hashrepeats{$chr}{$start}{"strand"}},$strand);
	    }
	    else{
		$hashrepeats{$chr}{$start}={"end"=>[$end],"classfam"=>[$classfam],"strand"=>[$strand]};
	    }
	}
	else{
	    $hashrepeats{$chr}={$start=>{"end"=>[$end],"classfam"=>[$classfam],"strand"=>[$strand]}};
	}
	
    	$line=<$input>;
    }
    
    close($input);

    ### now order repeats
    
    foreach my $chr (keys %hashrepeats){
	$repeatmasker->{$chr}={"start"=>[],"end"=>[],"classfam"=>[],"strand"=>[]};
	
	my @uniquestart=keys %{$hashrepeats{$chr}};
	my @sortedstart=sort {$a<=>$b} @uniquestart;

	foreach my $start (@sortedstart){
	    my $nb=@{$hashrepeats{$chr}{$start}{"end"}};
	    for(my $i=0; $i<$nb; $i++){

		if(${$hashrepeats{$chr}{$start}{"classfam"}}[$i] eq ""){
		    print "Found null type repeat while sorting: ".$chr." ".$start." ".${$hashrepeats{$chr}{$start}{"end"}}[$i]."\n";
		    exit(1);
		}

		push(@{$repeatmasker->{$chr}{"start"}},$start);
		push(@{$repeatmasker->{$chr}{"end"}},${$hashrepeats{$chr}{$start}{"end"}}[$i]);
		push(@{$repeatmasker->{$chr}{"classfam"}},${$hashrepeats{$chr}{$start}{"classfam"}}[$i]);
		push(@{$repeatmasker->{$chr}{"strand"}},${$hashrepeats{$chr}{$start}{"strand"}}[$i]);
	    }
	}
    }
}

####################################################################

sub extractOverlap{
    my $coords1=$_[0]; ## ordered coordinates
    my $coords2=$_[1]; ## ordered coordinates
    my $margin=$_[2];
    my $type=$_[3];
    my $overlap=$_[4];

    foreach my $chr (keys %{$coords1}){
	if(exists $coords2->{$chr}){
	    my $nbex1=@{$coords1->{$chr}{"start"}};
	    my $nbex2=@{$coords2->{$chr}{"start"}};
	    
	    my $firstj=0;
	    
	    for(my $i=0; $i<$nbex1; $i++){
		
		my $start1=${$coords1->{$chr}{"start"}}[$i]-$margin;
		my $end1=${$coords1->{$chr}{"end"}}[$i]+$margin;

		my $strand1=${$coords1->{$chr}{"strand"}}[$i];
		my $id1=$chr.":".$start1."-".$end1.":".$strand1;
			
		my $j=$firstj;
		
		while($j<$nbex2 && ${$coords2->{$chr}{"end"}}[$j]<$start1){
		    $j++;
		}
		
		$firstj=$j;
		
		while($j<$nbex2 && ${$coords2->{$chr}{"start"}}[$j]<=$end1){
		    
		    my $start2=${$coords2->{$chr}{"start"}}[$j];
		    my $end2=${$coords2->{$chr}{"end"}}[$j];
		    my $strand2=${$coords2->{$chr}{"strand"}}[$j];
		    my $classfam2=${$coords2->{$chr}{"classfam"}}[$j];
		    my $id2=$chr.":".$start2."-".$end2.":".$strand2;
		  		    
		    if(($strand1 eq $strand2 && $type eq "sense") || ($strand1 ne $strand2 && $type eq "antisense") || ($type eq "any")){
			my $M=max($start1,$start2);
			my $m=min($end1,$end2);
			
			if($M<=$m){
			    if(exists $overlap->{$id1}){
				$overlap->{$id1}{$id2}={"start"=>$M,"end"=>$m, "classfam"=>$classfam2};
			    }
			    else{
				$overlap->{$id1}={$id2=>{"start"=>$M,"end"=>$m,"classfam"=>$classfam2}};
			    }
			}
		    }
		    
		    $j++;
		}
	    }
	}
    }
}

##############################################################

sub printHelp{

    my $parnames=$_[0];
    my $parvalues=$_[1];
    
    print "\n";
    print "This script compares two repeat annotations.\n";
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
$parameters{"pathAnnotation1"}="NA";
$parameters{"pathAnnotation2"}="NA";
$parameters{"pathOutput"}="NA";

## "pathRepeatMasker", 

my @defaultpars=("pathAnnotation1", "pathAnnotation2", "pathOutput");
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

print "Reading repeat annotations...\n";

my %repeat1;
readRepeatMasker($parameters{"pathAnnotation1"}, \%repeat1);

my %repeat2;
readRepeatMasker($parameters{"pathAnnotation2"}, \%repeat2);

print "Done.\n";

##############################################################

print "Extracting overlap...\n";

my %overlap;

extractOverlap(\%repeat1, \%repeat2, 0, "sense", \%overlap);
 
print "Done.\n";

##############################################################

print "Writing output...\n";

open(my $output, ">".$parameters{"pathOutput"});

print $output "Chr\tStart\tEnd\tStrand\tType\tIntersectingElements\tIntersectingTypes\n";

foreach my $chr (keys %repeat1){
    my $nb=@{$repeat1{$chr}{"start"}};

    for(my $i=0; $i<$nb; $i++){
	my $start=${$repeat1{$chr}{"start"}}[$i];
	my $end=${$repeat1{$chr}{"end"}}[$i];
	my $strand=${$repeat1{$chr}{"strand"}}[$i];
	my $classfam=${$repeat1{$chr}{"classfam"}}[$i];

	my $id=$chr.":".$start."-".$end.":".$strand;

	my $idintersect="NA";
	my $famintersect="NA";
	
	if(exists $overlap{$id}){
	    $idintersect="";
	    $famintersect="";

	    foreach my $id2 (keys %{$overlap{$id}}){
		$idintersect.=$id2.";";
		$famintersect.=$overlap{$id}{$id2}{"classfam"}.";";
	    }

	    chop $idintersect;
	    chop $famintersect;
	}

	print $output $chr."\t".$start."\t".$end."\t".$strand."\t".$classfam."\t".$idintersect."\t".$famintersect."\n";
    }
}

close($output);

print "Done.\n";

##############################################################
##############################################################
