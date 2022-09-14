use List::Util qw(first max maxstr min minstr reduce shuffle sum);
use strict;

#############################################################################
#############################################################################

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

################################################################################

sub readDiamondResults{
    my $path=$_[0];
    my $minprotfr=$_[1];
    my $proteinlengths=$_[2];
    my $maxeval=$_[3];
    my $maxgapfraction=$_[4];
    my $kept=$_[5];

    open(my $input, $path);
    my $line=<$input>;

    my $nbgapped=0;
    my $nblargeeval=0;
    my $nbsmallprot=0;

    while($line){
	chomp $line;
	my @s=split("\t", $line);

	my $qid=$s[0];
	my $tid=$s[1];
	my $alnlen=$s[3]+0;
	my $qstart=$s[6]+0;
	my $qend=$s[7]+0;
	my $tstart=$s[8]+0;
	my $tend=$s[9]+0;
	my $evalue=$s[10]+0.0;
	my $totgaps=$s[12]+0; ## custom format

	if($tstart>$tend){
	    print "Weird! blastp results should not be on - strand.\n";
	    print $line."\n";
	    exit(1);
	}

	my $qlen=$qend-$qstart+1;
	my $tlen=$tend-$tstart+1;
	
	if($evalue < $maxeval){
	    my $gapfr=($totgaps+0.0)/($alnlen+0.0);
	    
	    if($gapfr>1){
		print "Weird ! gap fraction is larger than 1.\n";
		print $line."\n";
		exit(1);
	    }
	    
	    if($gapfr<$maxgapfraction){
		if(exists $proteinlengths->{$tid}){
		    my $thisprotfr=($alnlen-$totgaps+0.0)/($proteinlengths->{$tid}+0.0);
		    
		    if($thisprotfr>=$minprotfr){
			if(!exists $kept->{$qid}){
			    $kept->{$qid}=1;
			}
		    } else{
			$nbsmallprot++;	
		    }
		} else{
		    print "Weird! cannot find protein length for ".$tid."\n";
		    exit(1);
		}
	    } else{
		$nbgapped++;
	    }
	} else{
	    $nblargeeval++;
	}
	
	$line=<$input>;
    }
    
    close($input);

    print "Excluded ".$nbgapped." gapped entries (gap fraction larger than ".$maxgapfraction.").\n";
    print "Excluded ".$nblargeeval." entries with large e-value.\n";
    print "Excluded ".$nbsmallprot." entries that had small protein fractions.\n";
    
}

################################################################################

sub printHelp{

    my $parnames=$_[0];
    my $parvalues=$_[1];
    
    print "\n";
    print "This script extract statistics from diamond blastp results. \n";
    print "\n";
    print "Options:\n";
    
    foreach my $par (@{$parnames}){
	print "--".$par."  [  default value: ".$parvalues->{$par}."  ]\n";
    }
    print "\n";
}

#################################################################################
#################################################################################

my %parameters;

$parameters{"speciesList"}="NA";
$parameters{"pathsFastaProteins"}="NA";
$parameters{"pathsDiamondResults"}="NA";
$parameters{"minProteinFraction"}="NA";
$parameters{"maxEValue"}="NA";
$parameters{"maxGapFraction"}="NA";
$parameters{"pathOutput"}="NA";

my %defaultvalues;
my @defaultpars=("speciesList","pathsFastaProteins", "pathsDiamondResults", "minProteinFraction", "maxEValue", "maxGapFraction", "pathOutput");

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

###################################################################################
###################################################################################

print "Reading protein sequences...\n";

my @species=split(",", $parameters{"speciesList"});
my @pathsfasta=split(",", $parameters{"pathsFastaProteins"});
my @pathsblast=split(",", $parameters{"pathsDiamondResults"});

my $nbsp=@species;
my $nbfasta=@pathsfasta;
my $nbblast=@pathsblast;

if($nbsp!=$nbfasta || $nbsp!=$nbblast){
    print "Weird! saw ".$nbsp." species, ".$nbfasta." fasta files, ".$nbblast." blast files.\n";
    exit(1);
}

my %proteinlengths;

for(my $i=0; $i<$nbsp; $i++){
    my $sp=$species[$i];
    my $path=$pathsfasta[$i];

    my %proteins;
    readFasta($path, \%proteins);

    my $nbprot=keys %proteins;

    print "There are ".$nbprot." proteins for ".$sp."\n";

    $proteinlengths{$sp}={};
    foreach my $id (keys %proteins){
	my $length=length $proteins{$id};
	$proteinlengths{$sp}{$id}=$length;
    }
}

print "Done.\n";

###################################################################################

print "Reading and filtering tblastn results...\n";

my $minprotfr=$parameters{"minProteinFraction"}+0.0;
print "minimum protein length fraction: ".$minprotfr."\n";

my $maxeval=$parameters{"maxEValue"}+0.0;
print "maximum e-value: ".$maxeval."\n";

my $maxgapfr=$parameters{"maxGapFraction"}+0.0;
print "maximum gap fraction: ".$maxgapfr."\n";

my %kepttranscripts;

for(my $i=0; $i<$nbsp; $i++){
    my $sp=$species[$i];
    my $path=$pathsblast[$i];

    print "Reading tblastn results from ".$path." for ".$sp."\n";
    
    readDiamondResults($path, $minprotfr, $proteinlengths{$sp}, $maxeval, $maxgapfr, \%kepttranscripts);

    my $nbkept=keys %kepttranscripts;

    print "There are ".$nbkept." transcripts after reading data for ".$sp."\n";
}
    
print "Done.\n";

###################################################################################

print "Writing output...\n";

open(my $output, ">".$parameters{"pathOutput"});

print $output "TranscriptID\n";

foreach my $id (keys %kepttranscripts){
    print $output $id."\n";
}

close($output);

print "Done.\n";

###################################################################################
