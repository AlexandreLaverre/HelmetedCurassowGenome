#!/usr/bin/perl

use List::Util qw(first max maxstr min minstr reduce shuffle sum);
use strict;

##############################################################

sub readChromosomeNames{
    my $pathin=$_[0];
    my $chrnames=$_[1];

    open(my $input, $pathin);
    my $line=<$input>;

    while($line){
	chomp $line;
	my @s=split("\t", $line);
	my $old=$s[0];
	my $new=$s[1];
	$chrnames->{$old}=$new;
	
	$line=<$input>;
    }
    
    close($input);
}

##############################################################

sub printHelp{

    my $parnames=$_[0];
    my $parvalues=$_[1];
    
    print "\n";
    print "This script renames chromosomes.\n";
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
$parameters{"pathChromosomeCorrespondence"}="NA";
$parameters{"pathGenomeSequence"}="NA";
$parameters{"pathOutput"}="NA";

## "pathRepeatMasker", 

my @defaultpars=("pathChromosomeCorrespondence", "pathGenomeSequence", "pathOutput");
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

print "Reading chromosome correspondence...\n";

my %chrnames; 
readChromosomeNames($parameters{"pathChromosomeCorrespondence"}, \%chrnames);
  
print "Done.\n";

##############################################################

print "Reading genome sequence and writing output...\n";

open(my $input, $parameters{"pathGenomeSequence"});
open(my $output, ">".$parameters{"pathOutput"});

my $line=<$input>;

while($line){
    chomp $line;
    my $prefix=substr $line,0,1;

    if($prefix eq ">"){
	my $oldchr=substr $line,1;

	if(exists $chrnames{$oldchr}){
	    my $newchr=$chrnames{$oldchr};
	    print $output ">".$newchr."\n";
	} else{
	    print $output $line."\n";
	}
    } else{
	print $output $line."\n";
    }
    
    $line=<$input>;
}

close($output);
close($input);

print "Done.\n";

##############################################################
