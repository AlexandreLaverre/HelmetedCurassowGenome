use List::Util qw(first max maxstr min minstr reduce shuffle sum);
use strict;

##############################################################

sub readOrthogroups{
    my $pathin=$_[0];
    my $ogmem=$_[1];
    my $oglist=$_[2];

    open(my $input, $pathin);
    my $line=<$input>;
    chomp $line;
    my @s=split("\t", $line);
    my %header;
    my @species;

    for(my $i=0; $i<@s; $i++){
	$header{$s[$i]}=$i;

	if($i>3){
	    push(@species, $s[$i]);
	}
    }
    
    $line=<$input>;

    while($line){
	chomp $line;
	my @s=split("\t", $line);

	my $hog=$s[$header{"HOG"}];
	my $og=$s[$header{"OG"}];

	my $ogid=$hog."_".$og;

	$oglist->{$ogid}=[];
	
	foreach my $sp (@species){
	    my @t=split(", ",$s[$header{$sp}]);
	    
	    foreach my $gid (@t){
		if($gid ne ""){
		    my $spid=$sp.":".$gid;
		    
		    push(@{$oglist->{$ogid}}, $spid);
		    
		    if(exists $ogmem->{$spid}){
			print "Weird, already saw ".$spid."\n";
			exit(1);
		    } else{
			$ogmem->{$spid}=$ogid;
		    }
		}
	    }
	}
	
	$line=<$input>;
    }

    close($input);
}

##############################################################

sub compareOrthoGroups{
    my $memN0=$_[0];
    my $listN0=$_[1];
    my $memN1=$_[2];
    my $listN1=$_[3];
    my $N0N1=$_[4];
    my $N1N0=$_[5];
	    
    foreach my $idN1 (keys %{$listN1}){
	my %N0;
	
	foreach my $g (@{$listN1->{$idN1}}){
	    if(exists $memN0->{$g}){
		my $idN0=$memN0->{$g};
		$N0{$idN0}=1;
	    } else{
		print "Weird! found ".$g." in N1 but not in N0.\n";
	    }
	}
	
	my @uniqueN0=keys %N0;
	my $lenN0=@uniqueN0;
	
	if($lenN0>1){
	    print "Weird!! ".$idN1." is split in N0: ".join(", ",@uniqueN0)."\n";
	    exit(1);
	}
	
	my $idN0=$uniqueN0[0];
	
	$N1N0->{$idN1}=$idN0;

	if(exists $N0N1->{$idN0}){
	    push(@{$N0N1->{$idN0}}, $idN1);
	} else{
	    $N0N1->{$idN0}=[$idN1];
	}
    }
}

##############################################################

sub printHelp{

    my $parnames=$_[0];
    my $parvalues=$_[1];
    
    print "\n";
    print "This script extracts correspondence between sets of orthogroups. \n";
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

$parameters{"pathN0"}="NA";
$parameters{"pathN1"}="NA";
$parameters{"pathN2"}="NA";
$parameters{"pathOutput"}="NA";


my %defaultvalues;
my @defaultpars=("pathN0", "pathN1",  "pathN2",  "pathOutput");

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

print "Reading orthogroups...\n";

my %memN0;
my %listN0;

readOrthogroups($parameters{"pathN0"}, \%memN0, \%listN0);

my $nbo = keys %listN0;
my $nbg = keys %memN0;

print "There are ".$nbo." orthogroups and ".$nbg." genes in N0.\n";

my %memN1;
my %listN1;

readOrthogroups($parameters{"pathN1"}, \%memN1, \%listN1);

my $nbo = keys %listN1;
my $nbg = keys %memN1;

print "There are ".$nbo." orthogroups and ".$nbg." genes in N1.\n";


my %memN2;
my %listN2;

readOrthogroups($parameters{"pathN2"}, \%memN2, \%listN2);

my $nbo = keys %listN2;
my $nbg = keys %memN2;

print "There are ".$nbo." orthogroups and ".$nbg." genes in N2.\n";
 
print "Done.\n";

##############################################################

print "Comparing orthogroups...\n";

my %N1N0;
my %N0N1;

my %N2N0;
my %N0N2;

compareOrthoGroups(\%memN0, \%listN0, \%memN1, \%listN1, \%N0N1, \%N1N0);

my $nbN0N1=keys %N0N1;

print "There are ".$nbN0N1." N0 orthogroups with correspondence in N1.\n";

compareOrthoGroups(\%memN0, \%listN0, \%memN2, \%listN2, \%N0N2, \%N2N0);


my $nbN0N2=keys %N0N2;

print "There are ".$nbN0N2." N0 orthogroups with correspondence in N2.\n";
  
print "Done.\n";

##############################################################

print "Checking correspondence and writing output...\n";

open(my $output, ">".$parameters{"pathOutput"});

print $output "N1\tN2\tN0\n";

foreach my $idN0 (keys %listN0){
    if(exists $N0N1{$idN0} && exists $N0N2{$idN0}){
	my $nbN1=@{$N0N1{$idN0}};
	my $nbN2=@{$N0N2{$idN0}};
	
	if($nbN1==1 && $nbN2==1){
	    my $idN1=${$N0N1{$idN0}}[0];
	    my $idN2=${$N0N2{$idN0}}[0];
	    
	    print $output $idN1."\t".$idN2."\t".$idN0."\n";
	}
    } 
}

close($output);

print "Done.\n";

##############################################################
