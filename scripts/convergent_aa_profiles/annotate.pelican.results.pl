use List::Util qw(first max maxstr min minstr reduce shuffle sum);
use strict;

##############################################################

sub readOrthoGroupAnnotation{
    my $pathin=$_[0];
    my $ogannot=$_[1];
      
    open(my $input, $pathin);
    
    my $line=<$input>;
    chomp $line;
    my @s=split("\t", $line);
    my %header;

    for(my $i=0; $i<@s; $i++){
	$header{$s[$i]}=$i;
    }

    $line=<$input>;

    while($line){
	chomp $line;
	my @s=split("\t", $line);

	my $hog=$s[$header{"HOG"}];
	my $og=$s[$header{"OG"}];
	
	my $ogid=$hog."_".$og;

	my $refid=$s[$header{"Reference.GeneID"}];
	my $refname=$s[$header{"Reference.GeneName"}];
	my $humid=$s[$header{"HumanOrthologue.GeneID"}];
	my $humname=$s[$header{"HumanOrthologue.GeneName"}];

	if(exists $ogannot->{$ogid}){
	    print "Weird, already saw ".$ogid."\n";
	    exit(1);
	}

	$ogannot->{$ogid}={"refgeneid"=>$refid, "refgenename"=>$refname, "humangeneid"=>$humid, "humangenename"=>$humname};
	
	$line=<$input>;
    }

    close($input);
}

##############################################################

sub readOrthoGroupCorrespondence{
    my $pathin=$_[0];
    my $ogcorr=$_[1];

    open(my $input, $pathin);
    
    my $line=<$input>;
    chomp $line;
    my @s=split("\t", $line);
    my %header;

    for(my $i=0; $i<@s; $i++){
	$header{$s[$i]}=$i;
    }

    $ogcorr->{"N1N2"}={};
    $ogcorr->{"N2N1"}={};
    
    $line=<$input>;

    while($line){
	chomp $line;
	my @s=split("\t", $line);

	my $N1=$s[$header{"N1"}];
	my $N2=$s[$header{"N2"}];

	if(exists $ogcorr->{"N1N2"}{$N1}){
	    print "Weird! alreads saw N1 id ".$N1."\n";
	    exit(1);
	}
	
	$ogcorr->{"N1N2"}{$N1}=$N2;

	if(exists $ogcorr->{"N2N1"}{$N2}){
	    print "Weird! alreads saw N2 id ".$N2."\n";
	    exit(1);
	}
	
	$ogcorr->{"N2N1"}{$N2}=$N1;
	
	$line=<$input>;
    }
    
    close($input);
}

##############################################################

sub printHelp{

    my $parnames=$_[0];
    my $parvalues=$_[1];
    
    print "\n";
    print "This script adds information to Pelican results. \n";
    print "\n";
    print "Options:\n";
    
    foreach my $par (@{$parnames}){
	print "--".$par."  [  default value: ".$parvalues->{$par}."  ]\n";
    }
    print "\n";
}

##############################################################
##############################################################

my %parameters;

$parameters{"pathPelicanResults"}="NA";
$parameters{"orthogroup"}="NA";
$parameters{"pathOrthoGroupAnnotation"}="NA";
$parameters{"pathOrthoGroupCorrespondence"}="NA";
$parameters{"pathOutput"}="NA";


my %defaultvalues;
my @defaultpars=("pathPelicanResults", "orthogroup", "pathOrthoGroupAnnotation", "pathOrthoGroupCorrespondence", "pathOutput");

my %numericpars;

foreach my $par (keys %parameters){
    $defaultvalues{$par}=$parameters{$par};
}

## check if help was asked 

foreach my $arg (@ARGV){
    if($arg eq "--help"){
	printHelp(\@defaultpars, \%defaultvalues);
	exit(0);
    }
}

## check new parameters

my $nbargs=@ARGV;

for(my $i=0; $i<$nbargs; $i++){
    my $arg=$ARGV[$i];
    $arg=substr $arg,2;
    my @s=split("=",$arg);
    my $parname=$s[0];
    my $parval=$s[1];
    
    if(exists $parameters{$parname}){
	$parameters{$parname}=$parval;
	
	if(exists $numericpars{$parname}){
	    $parameters{$parname}=$parval+0.0;
	}
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

print "Reading orthogroup annotation...\n";

my %ogannot;

readOrthoGroupAnnotation($parameters{"pathOrthoGroupAnnotation"}, \%ogannot);

my $nbog=keys %ogannot;

print "There are ".$nbog." annotated orthogroups.\n";

print "Done.\n";

##############################################################

print "Reading orthogroup correspondence...\n";

my %ogcorr;

readOrthoGroupCorrespondence($parameters{"pathOrthoGroupCorrespondence"}, \%ogcorr);

print "Done.\n";

##############################################################

print "Reading Pelican results and writing output...\n";

my $thisog=$parameters{"orthogroup"};
my $otherog;

if($thisog eq "N1"){
    $otherog="N2";
} else{
    if($thisog eq "N2"){
	$otherog="N1";
    } else{
	print "Don't know what to do for orthogroup ".$thisog."\n";
	exit(1);
    }
}


print "This orthogroup is ".$thisog.", other orthogroup is ".$otherog."\n";


my $oc=$thisog.$otherog;

open(my $input, $parameters{"pathPelicanResults"});
open(my $output, ">".$parameters{"pathOutput"});

my $line=<$input>;

chomp $line;
my @s=split("\t", $line);

my %header;
for(my $i=0; $i<@s; $i++){
    $header{$s[$i]}=$i;
}

my $lineout=$line."\tReferenceGeneID\tReferenceGeneName\tHumanGeneID\tHumanGeneName\tOtherOrthoGroupID\n";

$line=<$input>;

while($line){
    chomp $line;
    my @s=split("\t", $line);

    my $thisid=$s[$header{"alignment"}];

    my $idog="NA";

    if(exists $ogannot{$thisid}){
	$idog=$thisid;
    } else{
	print "Weird! cannot find orthogroup id in ".$line."\n";
    }

    my $geneid=$ogannot{$idog}{"refgeneid"};
    my $genename=$ogannot{$idog}{"refgenename"};
    my $humid=$ogannot{$idog}{"humangeneid"};
    my $humname=$ogannot{$idog}{"humangenename"};
    
    my $otherog="NA";

    if(exists $ogcorr{$oc}{$idog}){
	$otherog=$ogcorr{$oc}{$idog};
    }
    
    my $lineout=$line."\t".$geneid."\t".$genename."\t".$humid."\t".$humname."\t".$otherog."\n";

    print $output $lineout;
    
    $line=<$input>;
}

close($input);
close($output);

print "Done.\n";

##############################################################
