use List::Util qw(first max maxstr min minstr reduce shuffle sum);
use strict;

##############################################################
##############################################################

sub readGeneNames{
    my $pathin=$_[0];
    my $names=$_[1];

    my $input;
    my @s=split("\\.", $pathin);
    my $ext=$s[-1];

    if($ext eq "gz"){
	open($input, "zcat $pathin |");
    } else{
	open($input, $pathin);
    }

    my $line=<$input>;

    my %noname;

    while($line){
	chomp $line;
	my @s=split("\t", $line);

	my $type=$s[2];

	if($type eq "gene"){
	    my $info=$s[8];
	    my @infoarray=split(";", $info);

	    my $id=findInfo("ID=", \@infoarray);
	    my $name=findInfo("Name=", \@infoarray);

	    if($id eq "NA" || $name eq "NA"){
		if($id ne "NA"){
		    $noname{$id}=1;
		}
	    } else{
		my @u=split(":", $id);
		my $geneid=$u[1];
		$names->{$geneid}=$name;
	    }
	}

	$line=<$input>;
    }

    my $nbnotok=keys %noname;
    my $nbok=keys %{$names};

    print "Found names for ".$nbok." genes, ".$nbnotok." genes with no name\n";
}

##############################################################

sub findInfo{
    my $pattern=$_[0];
    my $array=$_[1];

    my $res="NA";

    my @grepres=grep(/${pattern}/,@{$array});

    if(@grepres==1){
	my @t=split("=", $grepres[0]);
	$res=$t[1];
    }

    return $res;
}

##############################################################

sub printHelp{

    my $parnames=$_[0];
    my $parvalues=$_[1];

    print "\n";
    print "This script extracts gene names from GFF files.\n";
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
$parameters{"pathGFF"}="NA";
$parameters{"pathOutput"}="NA";

my @defaultpars=("pathGFF", "pathOutput");


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

print "Reading gene names...\n";

my %names;
readGeneNames($parameters{"pathGFF"}, \%names);

my $nbg=keys %names;
print "Saw ".$nbg." genes with names.\n";

print "Done.\n";

##############################################################

print "Writing output...\n";

open(my $output, ">".$parameters{"pathOutput"});

foreach my $gene (keys %names){
    print $output $gene."\t".$names{$gene}."\n";
}

close($output);

print "Done.\n";

##############################################################
