use List::Util qw(first max maxstr min minstr reduce shuffle sum);
use strict;

##############################################################
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
	    my $info=substr $line,1;
	    my @s=split(" ",$info);

	    my $id=$s[0]; ## if we couldn't find protein_id or gene_id, we keep CDS id
	    
	    foreach my $item (@s){
		my $prefix1=substr $item, 0, 5;

		if($prefix1 eq "gene:"){
		    my @t=split(":", $item);
		    $id=$t[1];
		    last;
		}

		my $prefix2=substr $item, 0, 12;

		if($prefix2 eq "[protein_id="){
		    my @t=split("=", $item);
		    my @u=split("]", $t[1]);
		    $id=$u[0];
		    last;
		}
	    }

	   
	    if(exists $reffasta->{$id}){
		print "Already saw ".$id."!";
		exit(1);
	    }

	     if($id eq "NA"){
		print "Weird! could not find sequence identifier in line ".$line."\n";
		exit(1);
	    }
	    
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

##################################################################################

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

############################################################################

sub readOrthogroups{
    my $path=$_[0];
    my $requiredsp=$_[1];
    my $minnbsp=$_[2];
    my $ortho=$_[3];

    ## we take only single-copy genes

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
    chomp $line;
    my @s=split("\t", $line);
    my %header;

    my @species;
    
    for(my $i=0; $i<@s; $i++){
	$header{$s[$i]}=$i;

	if($i>=3){
	    push(@species, $s[$i]);
	}
    }
    
    $line=<$input>;

    while($line){
	chomp $line;
	my @s=split("\t", $line);

	my $hog=$s[$header{"HOG"}];
	my $og=$s[$header{"OG"}];

	my $id=$hog."_".$og;

	foreach my $sp (@species){
	    my $g=$s[$header{$sp}];

	    if($g ne ""){
		my @genes=split(", ", $g);
		
		my $nbg=@genes;

		if($nbg==1){
		    if(exists $ortho->{$id}){
			$ortho->{$id}{$sp}=$genes[0];
		    } else{
			$ortho->{$id}={$sp=>$genes[0]};
		    }
		}
	    }
	}

	$line=<$input>;
    }

    close($input);

    my @allids=keys %{$ortho};
    my $nb=@allids;

    print "There were ".$nb." orthogroups before filtering.\n";

    foreach my $id (@allids){
	my $thisnb=keys %{$ortho->{$id}};

	if($thisnb < $minnbsp){
	    delete $ortho->{$id};
	} else{
	    my $hasrequired=1;

	    foreach my $sp (@{$requiredsp}){
		if(!exists $ortho->{$id}{$sp}){
		    $hasrequired=0;
		    last;
		}
	    }

	    if($hasrequired==0){
		delete $ortho->{$id};
	    }
	}
    }

    my $nbog=keys %{$ortho};

    print "Retained ".$nbog." orthogroups.\n";
}

############################################################################

sub printHelp{

    my $parnames=$_[0];
    my $parvalues=$_[1];

    print "\n";
    print "This script extracts coding sequences for predicted gene families.\n";
    print "\n";

    print "Options:\n";

    foreach my $par (@{$parnames}){
	print "--".$par."  [  default value: ".$parvalues->{$par}."  ]\n";
    }
    print "\n";
}

##########################################################################
##########################################################################

## parameters

my %parameters;
$parameters{"speciesList"}="NA";
$parameters{"pathsCDS"}="NA";
$parameters{"pathOrthogroups"}="NA";
$parameters{"requiredSpecies"}="NA";
$parameters{"minNbSpecies"}="NA";
$parameters{"dirOutput"}="NA";

my @defaultpars=("speciesList", "pathsCDS", "pathOrthogroups", "requiredSpecies", "minNbSpecies", "dirOutput");


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

print "Reading CDS sequences.\n";

my @species=split(",", $parameters{"speciesList"});
my @paths=split(",", $parameters{"pathsCDS"});

my $nbsp=@species;
my $nbp=@paths;

if($nbsp!=$nbp){
    print "Saw ".$nbsp." species and ".$nbp." paths, something is wrong.\n";
    exit(1);
}

my %cds;

for(my $i=0; $i<$nbp; $i++){
    my $sp=$species[$i];
    my $path=$paths[$i];

    print "Reading data for ".$sp." from ".$path."\n";

    $cds{$sp}={};

    readFasta($path, $cds{$sp});

    my $nb=keys %{$cds{$sp}};
    
    print "There are ".$nb." sequences for ".$sp."\n";
}

print "Done.\n";

##############################################################

print "Reading orthogroups...\n";

my @required=split(",", $parameters{"requiredSpecies"});
my $nbreq=@required;

print "There are ".$nbreq." required species: ".join(", ",@required)."\n";

my $minnb=$parameters{"minNbSpecies"}+0;

print "Minimum number of species: ".$minnb."\n";

my %orthogroups;
readOrthogroups($parameters{"pathOrthogroups"}, \@required, $minnb,\%orthogroups);
   
print "Done.\n";
    
##############################################################

print "Writing output for unaligned CDS sequences...\n";

foreach my $id (keys %orthogroups){
    my $path=$parameters{"dirOutput"}."/".$id.".unaln.fa";

    open(my $output, ">".$path);

    foreach my $sp (keys %{$orthogroups{$id}}){
	my $seqid=$orthogroups{$id}{$sp};
	
	if(exists $cds{$sp}{$seqid}){
	    my $seq=uc $cds{$sp}{$seqid};
	    writeSequence($seq, $sp, $output);
	} else{
	    print $id."\n";
	    print "Cannot find ".$seqid." for ".$sp."\n";
	    exit(1);
	}
    }

    close($output);
}

print "Done.\n";

##############################################################
