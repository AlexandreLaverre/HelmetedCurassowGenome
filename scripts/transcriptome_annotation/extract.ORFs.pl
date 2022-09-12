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

sub readTBlastN{
    my $path=$_[0];
    my $minorflen=$_[1];
    my $minprotfr=$_[2];
    my $proteinlengths=$_[3];
    my $maxeval=$_[4];
    my $orf=$_[5];

    open(my $input, $path);
    my $line=<$input>;

    my $nbgapped=0;
    my $nblargeeval=0;
    my $nbmod3=0;
    my $nbsmallorf=0;
    my $nbsmallprot=0;

    while($line){
	chomp $line;
	my @s=split("\t", $line);

	my $qid=$s[0];
	my $tid=$s[1];
	my $qstart=$s[6]+0;
	my $qend=$s[7]+0;
	my $tstart=$s[8]+0;
	my $tend=$s[9]+0;
	my $evalue=$s[10]+0.0;

	my $qlen=$qend-$qstart+1;
	my $tlen=$tend-$tstart+1;
	
	if($evalue < $maxeval){
	    if($tlen%3==0){
		if($tlen==(3*$qlen)){
		    if($tlen>=$minorflen){
			if(exists $proteinlengths->{$qid}){
			    my $thisprotfr=($qlen+0.0)/($proteinlengths->{$qid}+0.0);
			    
			    if($thisprotfr>=$minprotfr){
				my $orfid=$tid.":".$tstart.":".$tend;

				if(!exists $orf->{$orfid}){
				    $orf->{$orfid}=1;
				}
			    } else{
				$nbsmallprot++;	
			    }
			} else{
			    print "Weird! cannot find protein length for ".$qid."\n";
			    exit(1);
			}
		    } else{
			$nbsmallorf++;
		    }
		} else{
		    $nbgapped++;
		}
	    } else{
		$nbmod3++;
	    }
	} else{
	    $nblargeeval++;
	}
	
	$line=<$input>;
    }
    
    close($input);

    print "Excluded ".$nbgapped." gapped entries.\n";
    print "Excluded ".$nblargeeval." entries with large e-value.\n";
    print "Excluded ".$nbmod3." entries that are not modulo 3.\n";
    print "Excluded ".$nbsmallorf." entries that had small ORFs.\n";
    print "Excluded ".$nbsmallprot." entries that had small protein fractions.\n";
    
}

################################################################################

sub extractConnectedORFs{
    my $orfs=$_[0];
    my $connected=$_[1];

    my %orfsbytx;

    foreach my $id (keys %{$orfs}){
	my @s=split(":", $id);
	my $tx=$s[0];
	my $start=$s[1]+0;
	my $end=$s[2]+0;

	if(exists $orfsbytx{$tx}){
	    push(@{$orfsbytx{$tx}{"start"}}, $start);
	    push(@{$orfsbytx{$tx}{"end"}}, $end);
	} else{
	    $orfsbytx{$tx}={"start"=>[$start], "end"=>[$end]};
	}
    }

    foreach my $tx (keys %orfsbytx){
	my $nbo=@{$orfsbytx{$tx}{"start"}};

	if($nbo>=2){
	    for(my $i=0; $i<($nbo-1); $i++){
		my $start1=${$orfsbytx{$tx}{"start"}}[$i];
		my $end1=${$orfsbytx{$tx}{"end"}}[$i];
		my $id1=$tx.":".$start1.":".$end1;
		
		for(my $j=($i+1); $j<$nbo; $j++){
		    my $start2=${$orfsbytx{$tx}{"start"}}[$j];
		    my $end2=${$orfsbytx{$tx}{"end"}}[$j];
		    my $id2=$tx.":".$start2.":".$end2;

		    my $M=max($start1, $start2);
		    my $m=min($end1, $end2);

		    if($M <= $m){
			my $diffstart=abs($start1-$start2);

			if($diffstart%3==0){
			    if(exists $connected->{$id1}){
				$connected->{$id1}{$id2}=1;
			    } else{
				$connected->{$id1}={$id1=>1, $id2=>1};
			    }

			    if(exists $connected->{$id2}){
				$connected->{$id2}{$id1}=1;
			    } else{
				$connected->{$id2}={$id1=>1, $id2=>1};
			    }
			}
		    }
		}
	    }
	}
    }
}

##############################################################

sub extractClusters{
    my $refconnected=$_[0];
    my $refclusters=$_[1];
    my $refclustid=$_[2];

    my $nbconnected=keys %{$refconnected};

    my $round=0;

    while($nbconnected>0){
	
	foreach my $key (keys %{$refconnected}){
	    addToCluster($refconnected,$refclusters,$refclustid,$key);
	}
	
	$round++;
	
	$nbconnected=keys %{$refconnected};
    }
}

##############################################################

sub addToCluster{
    my $refconnected=$_[0];
    my $refclusters=$_[1];
    my $refclustid=$_[2];
    my $key=$_[3];

    if(exists $refconnected->{$key}){
	
	## find the cluster that contains this key
	
	my $indexcluster="NA";
	
	if(exists $refclustid->{$key}){
	    $indexcluster=$refclustid->{$key};
	}
	
	## if there isn't any
	
	if($indexcluster eq "NA"){
	    my $nbclusters=keys %{$refclusters};
	    $indexcluster=$nbclusters+1;
	    $refclusters->{$indexcluster}={$key=>1};
	    $refclustid->{$key}=$indexcluster;
	}
		
	foreach my $connection (keys %{$refconnected->{$key}}){

	    ## check if this island is already in the cluster
	    
	    if(!(exists $refclusters->{$indexcluster}{$connection})){
		$refclusters->{$indexcluster}{$connection}=1;
		$refclustid->{$connection}=$indexcluster;
		addToCluster($refconnected,$refclusters,$refclustid,$connection);
	    }
	}

	## after we've checked all of its connections, remove it from the connected islands

	delete $refconnected->{$key};
    }
    
}


################################################################################

sub printHelp{

    my $parnames=$_[0];
    my $parvalues=$_[1];
    
    print "\n";
    print "This script extract ORFs from tblastn results. \n";
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
$parameters{"pathsTBlastNResults"}="NA";
$parameters{"minORFLength"}="NA";
$parameters{"minProteinFraction"}="NA";
$parameters{"maxEValue"}="NA";
$parameters{"pathOutput"}="NA";

my %defaultvalues;
my @defaultpars=("speciesList", "pathsFastaProteins", "pathsTBlastNResults", "minORFLength", "minProteinFraction", "maxEValue", "pathOutput");

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
my @pathsblast=split(",", $parameters{"pathsTBlastNResults"});

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

my $minorflen=$parameters{"minORFLength"}+0.0;
print "minimum ORF length: ".$minorflen."\n";

my $minprotfr=$parameters{"minProteinFraction"}+0.0;
print "minimum protein length fraction: ".$minprotfr."\n";

my $maxeval=$parameters{"maxEValue"}+0.0;
print "maximum e-value: ".$maxeval."\n";

my %orf;

for(my $i=0; $i<$nbsp; $i++){
    my $sp=$species[$i];
    my $path=$pathsblast[$i];

    print "Reading tblastn results from ".$path." for ".$sp."\n";
    
    readTBlastN($path, $minorflen, $minprotfr, $proteinlengths{$sp}, $maxeval, \%orf);

    my $nborf=keys %orf;

    print "There are ".$nborf." ORFs after reading data for ".$sp."\n";
}
    
print "Done.\n";

###################################################################################

print "Extracting connected ORFs...\n";

my %connectedorf;

extractConnectedORFs(\%orf, \%connectedorf);

print "Done.\n";

###################################################################################

print "Extracting ORF clusters...\n";

my %orfclust;
my %clustid;

extractClusters(\%connectedorf, \%orfclust, \%clustid);

my $nbclust=keys %orfclust;

print "There are ".$nbclust." ORF clusters.\n";

print "Done.\n";

###################################################################################
