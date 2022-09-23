use List::Util qw(first max maxstr min minstr reduce shuffle sum);
use strict;

#############################################################################
#############################################################################

sub readORFs{
    my $pathin=$_[0];
    my $orf=$_[1];

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

	my $txid=$s[$header{"TranscriptID"}];
	my $strand=$s[$header{"Strand"}];
	my $start=$s[$header{"Start"}]+0;
	my $end=$s[$header{"End"}]+0;
	
	my $orfid=$txid.":".$start.":".$end.":".$strand;

	if(!exists $orf->{$orfid}){
	    $orf->{$orfid}=1;
	}
	
	$line=<$input>;
    }
    
    close($input);
}

#############################################################################

sub extractConnectedORFs{
    my $orfs=$_[0];
    my $connected=$_[1];

    my %orfsbytx;

    foreach my $id (keys %{$orfs}){
	my @s=split(":", $id);
	my $tx=$s[0];
	my $start=$s[1]+0;
	my $end=$s[2]+0;
	my $strand=$s[3];
	
	if(exists $orfsbytx{$tx}){
	    push(@{$orfsbytx{$tx}{"start"}}, $start);
	    push(@{$orfsbytx{$tx}{"end"}}, $end);
	    push(@{$orfsbytx{$tx}{"strand"}}, $strand);
	} else{
	    $orfsbytx{$tx}={"start"=>[$start], "end"=>[$end], "strand"=>[$strand]};
	}
    }

    foreach my $tx (keys %orfsbytx){
	my $nbo=@{$orfsbytx{$tx}{"start"}};

	if($nbo>=2){
	    for(my $i=0; $i<($nbo-1); $i++){
		my $start1=${$orfsbytx{$tx}{"start"}}[$i];
		my $end1=${$orfsbytx{$tx}{"end"}}[$i];
		my $strand1=${$orfsbytx{$tx}{"strand"}}[$i];
		my $id1=$tx.":".$start1.":".$end1.":".$strand1;
		
		for(my $j=($i+1); $j<$nbo; $j++){
		    my $start2=${$orfsbytx{$tx}{"start"}}[$j];
		    my $end2=${$orfsbytx{$tx}{"end"}}[$j];
		    my $strand2=${$orfsbytx{$tx}{"strand"}}[$j];
		    my $id2=$tx.":".$start2.":".$end2.":".$strand2;

		    if($strand1 eq $strand2){
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

sub extractLargestORF{
    my $clusters=$_[0];
    my $largestorf=$_[1];

    foreach my $index (keys %{$clusters}){
	my %txid;
	my %strands;
	my @starts;
	my @ends;

	foreach my $orfid (keys %{$clusters->{$index}}){
	    my @s=split(":", $orfid);
	    my $tx=$s[0];
	    my $start=$s[1]+0;
	    my $end=$s[2]+0;
	    my $strand=$s[3];

	    $txid{$tx}=1;
	    $strands{$strand}=1;

	    push(@starts, $start);
	    push(@ends, $end);
	}

	my $nbtx=keys %txid;
	my $nbstrands=keys %strands;

	if($nbtx==1 && $nbstrands==1){
	    my @alltx=keys %txid;
	    my $commontxid=$alltx[0];
	    
	    my @allstrands=keys %strands;
	    my $commonstrand=$allstrands[0];
	    
	    my $minstart=min @starts;
	    my $maxend=max @ends;

	    $largestorf->{$index}={"txid"=>$commontxid, "strand"=>$commonstrand, "start"=>$minstart, "end"=>$maxend};
	    
	} else{
	    print "Weird! different transcripts and strands for cluster.\n";
	    exit(1);
	}
    }
}

################################################################################

sub printHelp{

    my $parnames=$_[0];
    my $parvalues=$_[1];
    
    print "\n";
    print "This script combines ORF coordinates. \n";
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

$parameters{"pathsORFs"}="NA";
$parameters{"pathOutput"}="NA";

my %defaultvalues;
my @defaultpars=("pathsORFs", "pathOutput");

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

print "Reading ORFs...\n";

my @paths=split(",", $parameters{"pathsORFs"});

my %orf;

foreach my $path (@paths){
    readORFs($path, \%orf);

    my $nborf=keys %orf;

    print "There are ".$nborf." ORFs after reading data from ".$path."\n";
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

print "Extracting largest ORFs for clusters...\n";

my %largestorf;

extractLargestORF(\%orfclust, \%largestorf);

print "Done.\n";

###################################################################################

print "Writing output...\n";

open(my $output, ">".$parameters{"pathOutput"});

print $output "TranscriptID\tStrand\tStart\tEnd\n";

foreach my $clustid (keys %orfclust){
    print $output $largestorf{$clustid}{"txid"}."\t".$largestorf{$clustid}{"strand"}."\t".$largestorf{$clustid}{"start"}."\t".$largestorf{$clustid}{"end"}."\n";
}

close($output);

print "Done.\n";

###################################################################################
