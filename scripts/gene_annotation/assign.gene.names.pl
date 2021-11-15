use List::Util qw(first max maxstr min minstr reduce shuffle sum);
use strict;

###########################################################################
###########################################################################

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

    while($line){
	chomp $line;
	my @s=split("\t", $line);

	my $geneid=$s[0];
	my $name=$s[1];

	if(exists $names->{$geneid} && $names->{$geneid} ne $name){
	    print "Weird! already saw ".$geneid."\n";
	    exit(1);
	}

	$names->{$geneid}=$name;

	$line=<$input>;
    }

    close($input);
}

###########################################################################

sub readOrtho{
    my $pathin=$_[0];
    my $ref=$_[1];
    my $tg=$_[2];
    my $ortho=$_[3];

    my $input;
    my @s=split("\\.", $pathin);
    my $ext=$s[-1];

    if($ext eq "gz"){
	open($input, "zcat $pathin |");
    } else{
	open($input, $pathin);
    }

    my $line=<$input>; ## header
    chomp $line;
    my @s=split("\t", $line);
    my %header;

    for(my $i=0; $i<@s; $i++){
	$header{$s[$i]}=$i;

	print $s[$i]." ".$i."\n";
    }

    if(!exists $header{$tg}){
	print "Weird! cannot find ".$tg." column\n";
	print $pathin."\n";
	exit(1);
    }

    if(!exists $header{$ref}){
	print "Weird! cannot find ".$ref." column\n";
	print $pathin."\n";
	exit(1);
    }


    $line=<$input>;

    while($line){
	chomp $line;
	my @s=split("\t", $line);

	my $idref=$s[$header{$ref}];
	my $idtg=$s[$header{$tg}];

	my @sref=split(", ", $idref);
	my @stg=split(", ", $idtg);

	my $nbref=@sref;
	my $nbtg=@stg;

	if($nbref==1 && $nbtg==1){
	    my @tref=split("\\.", $idref);
	    my @ttg=split("\\.", $idtg);

	    my $generef=$tref[0];
	    my $genetg=$ttg[0];

	    if(exists $ortho->{$generef}){
		$ortho->{$generef}{$tg}=$genetg;
	    } else{
		$ortho->{$generef}={$tg=>$genetg};
	    }
	}

	$line=<$input>;
    }

    close($input);
}

###########################################################################

sub printHelp{

    my $parnames=$_[0];
    my $parvalues=$_[1];

    print "\n";
    print "This script assigns names to helmeted curassow genes based on 1 to 1 orthologues.\n";
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
$parameters{"refSpecies"}="NA";
$parameters{"tgSpecies"}="NA";
$parameters{"dirAnnot"}="NA";
$parameters{"dirOrtho"}="NA";
$parameters{"pathOutput"}="NA";

my @defaultpars=("refSpecies", "tgSpecies", "dirAnnot", "dirOrtho", "pathOutput");


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

print "Reading gene names for target species...\n";

my @tgsp=split(",", $parameters{"tgSpecies"});

my %genenames;

foreach my $tg (@tgsp){
    $genenames{$tg}={};

    my $path=$parameters{"dirAnnot"}."/".$tg.".genenames.txt";

    if(-e $path){
	readGeneNames($path,$genenames{$tg});

	my $nbg=keys %{$genenames{$tg}};

	print "Found ".$nbg." gene names for ".$tg."\n";
    } else{
	print "Cannot find names for ".$tg."\n";
	print $path."\n";
	exit(1);
    }
}

print "Done.\n";

##############################################################

print "Reading ortho pairs...\n";

my %ortho;

my $ref=$parameters{"refSpecies"};

foreach my $tg (@tgsp){
    my $path=$parameters{"dirOrtho"}."/".$ref."__v__".$tg.".tsv";

    if(-e $path){
	readOrtho($path,$ref, $tg, \%ortho);
    } else{
	print "Cannot find ortho for ".$tg."\n";
	print $path."\n";
	exit(1);
    }

    my $nbg=keys %ortho;

    print "There are ".$nbg." genes with 1 to 1 ortho after reading data for ".$tg.".\n";

}

print "Done.\n";

##############################################################

print "Assigning gene names and writing output...\n";

open(my $output, ">".$parameters{"pathOutput"});

my $line="GeneID";

foreach my $tg (@tgsp){
    $line.="\tGeneID.".$tg."\tGeneName.".$tg;
}

$line.="\tGeneName.".$ref;

print $output $line."\n";

my $nbnoconsensus=0;

foreach my $gene (keys %ortho){
    my %names;

    my $line=$gene;

    foreach my $tg (@tgsp){
	if(exists $ortho{$gene}{$tg}){
	    my $genetg=$ortho{$gene}{$tg};

	    print "looking at ortho ".$gene." ".$genetg." for ".$tg."\n";

	    if(exists $genenames{$tg}{$genetg}){
		my $name=$genenames{$tg}{$genetg};
		$line.="\t".$genetg."\t".$name;

		my $ucname=uc $name;

		$names{$ucname}=1;
	    } else{
		$line.="\t".$genetg."\tNA";
	    }
	} else{
	    $line.="\tNA\tNA";
	}
    }

    my $nbn=keys %names;

    if($nbn==1){
	my @names=keys %names;
	my $consensus=$names[0];

	$line.="\t".$consensus."\n";
    } else{
	$nbnoconsensus++;
	$line.="\tNA\n";
    }

    print $output $line;
}

close($output);

print "There were ".$nbnoconsensus." genes without a consensus name.\n";

print "Done.\n";

##############################################################
