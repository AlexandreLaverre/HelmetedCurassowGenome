use List::Util qw(first max maxstr min minstr reduce shuffle sum);
use strict;

##############################################################
##############################################################

sub readFasta{
    my $path=$_[0];
    my $reffasta=$_[1];
    my $duplicated=$_[2];
   
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
	    my $id=$s[0];
	  	   
	    if(exists $reffasta->{$id}){
		print "Already saw ".$id."!";
		$duplicated->{$id}=1;
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

sub printHelp{

    my $parnames=$_[0];
    my $parvalues=$_[1];

    print "\n";
    print "This script filters alignments.\n";
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
$parameters{"pathFastaInput"}="NA";
$parameters{"minSpecies"}="NA";
$parameters{"maxGapProportion"}="NA";
$parameters{"minUngappedLength"}="NA";
$parameters{"pathFastaOutput"}="NA";
$parameters{"pathOutputDiscarded"}="NA";

my %defaultvalues;
my @defaultpars=("pathFastaInput", "maxGapProportion", "minUngappedLength", "minSpecies", "pathFastaOutput", "pathOutputDiscarded");

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

my %aln;
my %duplicated;

readFasta($parameters{"pathFastaInput"}, \%aln, \%duplicated);
  
print "Done.\n";

######################################################################

print "Filtering alignment...\n";

my $nbdupli=keys %duplicated;
my $nbaln=keys %aln;

my $minsp=$parameters{"minSpecies"}+0;
my $maxpropgap=$parameters{"maxGapProportion"}+0.0;
my $minungaplen=$parameters{"minUngappedLength"}+0;

print "Minimum ".$minsp." species.\n";
print "Max ".$maxpropgap." gaps.\n";
print "Minimum ungapped length: ".$minungaplen."\n";

if($nbaln<$minsp){
    open(my $outerr, ">".$parameters{"pathOutputDiscarded"});
    print $outerr $nbaln." sequences.\n";
    close($outerr);
} else{
    if($nbdupli>($nbaln/3)){
	open(my $outerr, ">".$parameters{"pathOutputDiscarded"});
	print $outerr $nbdupli." duplicated sequences for ".$nbaln." species\n";
	close($outerr);
    } else{
	foreach my $id (keys %duplicated){
	    delete $aln{$id};
	}
	
	my $nbleft=keys %aln;
	
	if($nbleft>=$minsp){
	    my %toogapped;
	    
	    foreach my $sp (keys %aln){
		my $seq=$aln{$sp};
		my $totlen=length $seq;
		my $countgap= ($seq =~ tr/-//);
		my $countungap=$totlen-$countgap;
		my $propgap=($countgap+0.0)/($totlen+0.0);
		
		if($propgap>$maxpropgap && $countungap<$minungaplen){
		    $toogapped{$sp}=1;
		}
	    }
	    
	    foreach my $id (keys %toogapped){
		delete $aln{$id};
	    }
	    
	    my $nbleft2=keys %aln;
	    
	    if($nbleft2>=$minsp){
		open(my $output, ">".$parameters{"pathFastaOutput"});
		
		foreach my $sp (keys %aln){
		    my $seq=uc $aln{$sp}; ## remove masking
		    writeSequence($seq, $sp, $output);
		}
		close($output);
		
	    } else{
		open(my $outerr, ">".$parameters{"pathOutputDiscarded"});
		print $outerr $nbleft2." sequences left after removing gapped sequences.\n";
		close($outerr);
	    }
	} else{
	    open(my $outerr, ">".$parameters{"pathOutputDiscarded"});
	    print $outerr $nbleft." sequences left after removing duplicates.\n";
	    close($outerr);
	}
	
    }
}

print "Done.\n";

######################################################################
