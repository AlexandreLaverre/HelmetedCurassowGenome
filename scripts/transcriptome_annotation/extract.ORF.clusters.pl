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

sub writeSequence{
    my $sequence=$_[0];
    my $name=$_[1];
    my $output=$_[2];

    my $n=length $sequence;

    print $output ">".$name."\n";

    my $i=0;

    while($i<($n-60)){

        my $subseq=substr $sequence,$i,60;

        print $output $subseq ."\n";

        $i+=60;
    }

    if($i<$n){
        my $subseq=substr $sequence,$i;
        print $output $subseq ."\n";
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

sub readDiamondResults{
    my $path=$_[0];
    my $proteins=$_[1];
    my $minpcid=$_[2];
    my $minalnfr=$_[3];
    my $connected=$_[4];

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
	chomp $line;
	my @s=split("\t", $line);
	my $id1=$s[0];
	my $id2=$s[1];

	my $pcid=$s[2]+0.0;
	my $alnlen=$s[3]+0.0;

	if(exists $proteins->{$id1} && exists $proteins->{$id2}){
	    my $len1=length $proteins->{$id1};
	    my $len2=length $proteins->{$id2};

	    my $fraln1=($alnlen+0.0)/$len1;
	    my $fraln2=($alnlen+0.0)/$len2;

	    if($pcid>=$minpcid && ($fraln1>=$minalnfr || $fraln2>=$minalnfr)){
		if(exists $connected->{$id1}){
		    $connected->{$id1}{$id2}=1;
		} else{
		    $connected->{$id1}={$id2=>1, $id1=>1};
		}

		if(exists $connected->{$id2}){
		    $connected->{$id2}{$id1}=1;
		} else{
		    $connected->{$id2}={$id2=>1, $id1=>1};
		}
	    }

	} else{
	    print "Weird! cannot find protein sequences for ".$id1." and ".$id2."\n";
	    exit(1);
	}
	
	$line=<$input>;
    }
}

################################################################################

sub printHelp{

    my $parnames=$_[0];
    my $parvalues=$_[1];
    
    print "\n";
    print "This script extracts ORF clusters based on diamond results. \n";
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

$parameters{"pathInputCDS"}="NA";
$parameters{"pathInputProteins"}="NA";
$parameters{"pathDiamondResults"}="NA";
$parameters{"minPercentageIdentity"}="NA";
$parameters{"minAlignedLengthFraction"}="NA";
$parameters{"pathOutputClusters"}="NA";
$parameters{"pathOutputCDS"}="NA";
$parameters{"pathOutputProteins"}="NA";

my %defaultvalues;
my @defaultpars=("pathInputCDS","pathInputProteins", "pathDiamondResults", "minPercentageIdentity", "minAlignedLengthFraction", "pathOutputClusters", "pathOutputCDS", "pathOutputProteins");

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

print "Reading CDS and protein sequences...\n";

my %cds;
readFasta($parameters{"pathInputCDS"}, \%cds);

my $nb=keys %cds;
print "Found ".$nb." cds.\n";


my %proteins;
readFasta($parameters{"pathInputProteins"}, \%proteins);

my $nb=keys %proteins;
print "Found ".$nb." proteins.\n";

print "Done.\n";

###################################################################################

print "Reading diamond results and extracting connected sequences...\n";

my $minpcid = $parameters{"minPercentageIdentity"}+0.0;
my $minalnfr = $parameters{"minAlignedLengthFraction"}+0.0;

my %connected;

readDiamondResults($parameters{"pathDiamondResults"}, \%proteins, $minpcid, $minalnfr, \%connected);

print "Done.\n";

###################################################################################

print "Extracting protein clusters...\n";

my %clusters;
my %clustid;

extractClusters(\%connected, \%clusters, \%clustid);

print "Done.\n";

###################################################################################

print "Extracting representative ORFs and writing output...\n";

open(my $outputclust, ">".$parameters{"pathOutputClusters"});
open(my $outputcds, ">".$parameters{"pathOutputCDS"});
open(my $outputprotein, ">".$parameters{"pathOutputProteins"});

print $outputclust "IDCluster\tClusterMembers\tRepresentativeORF\n";

foreach my $idclust (keys %clusters){

    my %protlength;

    foreach my $idprot (keys %{$clusters{$idclust}}){
	if(exists $proteins{$idprot}){
	    my $seqprot=$proteins{$idprot};
	    my $len=length $seqprot;
	    
	    if(exists $protlength{$len}){
		push(@{$protlength{$len}}, $idprot);
	    } else{
		$protlength{$len}=[$idprot];
	    }
	}
	else{
	    print "Weird! cannot find ".$idprot." in proteins.\n";
	    exit(1);
	}
    }

    my $nbkept=keys %protlength;

    if($nbkept>0){
	my @uniquelengths=keys %protlength;
	my $maxlen = max @uniquelengths;
	
	my @idmax=@{$protlength{$maxlen}};
	my $repid=$idmax[0];
	
	print $outputclust $idclust."\t".join(",", keys %{$clusters{$idclust}})."\t".$repid."\n";
		
	if(exists $cds{$repid} && exists $proteins{$repid}){
	    my $repcds=$cds{$repid};
	    my $repprot=$proteins{$repid};
	    
	    writeSequence($repcds, $repid, $outputcds);
	    writeSequence($repprot, $repid, $outputprotein);
	} else{
	    print "Weird! cannot find ".$repid." in proteins or CDS sequences.\n";
	    exit(1);
	}
    } else{
	print "Weird, cluster ".$idclust." had no sequences.\n";
	exit(1);
    }
}
    
close($outputclust);
close($outputcds);
close($outputprotein);

print "Done.\n";

###################################################################################
###################################################################################
