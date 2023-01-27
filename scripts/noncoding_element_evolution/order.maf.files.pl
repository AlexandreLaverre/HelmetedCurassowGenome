use List::Util qw(first max maxstr min minstr reduce shuffle sum);
use strict;

#######################################################################

sub readAlignments{
    my $path=$_[0];
    my $refsp=$_[1];
    my $aln=$_[2];
    my $refcoords=$_[3];

    my @s=split("\\.",$path);
    my $ext=$s[-1];
    my $input;

    if($ext eq "gz"){
	open($input, "zcat $path |");
    } else{
	open($input, $path);
    }

    my $line=<$input>;

    my $currentindex=0;
    
    while($line){
	my $prefix=substr $line, 0, 1;

	if($prefix eq "a"){
	    $currentindex++;
	}

	if($prefix eq "s"){
	    if(exists $aln->{$currentindex}){
		push(@{$aln->{$currentindex}}, $line);
	    } else{
		$aln->{$currentindex}=[$line];
	    }

	    my @t=split("\t", $line);
	    my @u=split("\\.", $t[1]);
	    my $sp=$u[0];

	    if($sp eq $refsp){
		my $chr=$u[1];
		my $start=$t[2]+0;

		if(exists $refcoords->{$chr}){
		    if(exists $refcoords->{$chr}{$start}){
			push(@{$refcoords->{$chr}{$start}}, $currentindex);
		    } else{
			$refcoords->{$chr}{$start}=[$currentindex];
		    }
		} else{
		    $refcoords->{$chr}={$start=>[$currentindex]};
		}
	    }
	}
	
	$line=<$input>;
    }
    
    close($input);
}

#######################################################################

sub writeOrderedAlignments{
    my $aln=$_[0];
    my $refcoords=$_[1];
    my $pathoutput=$_[2];

    open(my $output, ">".$pathoutput);
    print $output "##maf version=1 scoring=N/A\n";

    foreach my $chr (keys %{$refcoords}){
	my @sortedpos=sort {$a<=>$b} (keys %{$refcoords->{$chr}});

	foreach my $pos (@sortedpos){
	    foreach my $index (@{$refcoords->{$chr}{$pos}}){
		print $output "\na\n";
		foreach my $line (@{$aln->{$index}}){
		    print $output $line;
		}
	    }
	}
    }
    
    close($output);
}


#######################################################################

sub printHelp{

    my $parnames=$_[0];
    my $parvalues=$_[1];

    print "\n";
    print "This script orders MAF files by reference species coordinates.\n";
    print "\n";

    print "Options:\n";

    foreach my $par (@{$parnames}){
	print "--".$par."  [  default value: ".$parvalues->{$par}."  ]\n";
    }
    print "\n";
}

#######################################################################
#######################################################################

my %parameters;
$parameters{"pathMAFInput"}="NA";
$parameters{"refSpecies"}="NA";
$parameters{"pathMAFOutput"}="NA";

my %defaultvalues;
my @defaultpars=("pathMAFInput", "refSpecies", "pathMAFOutput");

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

######################################################################
######################################################################

print "Reading alignment...\n";

my $refsp=$parameters{"refSpecies"};

my %aln;
my %refcoords;

readAlignments($parameters{"pathMAFInput"}, $refsp, \%aln, \%refcoords);
  
print "Done.\n";

######################################################################

print "Writing ordered alignments...\n";

writeOrderedAlignments(\%aln, \%refcoords, $parameters{"pathMAFOutput"});
 
print "Done.\n";

######################################################################
