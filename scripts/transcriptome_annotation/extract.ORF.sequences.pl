use List::Util qw(first max maxstr min minstr reduce shuffle sum);
use POSIX;
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

sub readGeneticCode{
    my $pathin=$_[0];
    my $code=$_[1];

    open(my $input, $pathin);

    my $line=<$input>;

    while($line){
	chomp $line;
	my @s=split("\t", $line);
	
	my $codon=$s[0];
	my $aa=$s[1];
	
	$code->{$codon}=$aa;

	$line=<$input>;
    }

    close($input);
}

################################################################################

sub translateCDS{
    my $cds=$_[0];
    my $geneticcode=$_[1];

    ## we allow in-frame stop codons, because of poor sequence quality
    
    my $protein="";
    # my $inframestop=0;
    my $size=length $cds;

    if($size%3!=0){
	print "Weird! CDS sequence is not multiple of 3.\n";
	exit(1);
    }

    for(my $i=0; $i<($size-2); $i+=3){
	my $codon=substr $cds, $i, 3;
	my $aa;
	    
	if(!exists $geneticcode->{$codon}){
	    print "Weird!! cannot find ".$codon." in the genetic code.\n";
	    $aa="X";
	}
	else{
	    $aa=$geneticcode->{$codon};
	}

	# if(($aa eq "*") && $i<($size-3)){
	#     $inframestop=1;
	# }

	$protein.=$aa;

    }
    
    return $protein;
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

################################################################################

sub readORFs{
    my $pathin=$_[0];
    my $hits=$_[1];

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
	
	my $id=$s[$header{"TranscriptID"}];
	my $start=$s[$header{"Start"}]+0;
	my $end=$s[$header{"End"}]+0;
	my $strand=$s[$header{"Strand"}];

	my $length=$end-$start+1;

	if($length%3!=0){
	    print "Weird!! length not multiple of 3 for ".$line.".\n";
	    exit(1);
	}

	if(exists $hits->{$id}){
	    push(@{$hits->{$id}{"start"}}, $start);
	    push(@{$hits->{$id}{"end"}}, $end);
	    push(@{$hits->{$id}{"strand"}}, $strand);
	}
	else{
	    $hits->{$id}={"start"=>[$start], "end"=>[$end], "strand"=>[$strand]};
	}

	$line=<$input>;
    }

    close($input);
}

#########################################################################################

sub reverseComplement{
    my $sequence=$_[0];

    my $rev=reverse $sequence;

    $rev=~s/A/X/g;
    $rev=~s/C/Y/g;
    $rev=~s/G/Z/g;
    $rev=~s/T/W/g;

    $rev=~s/X/T/g;
    $rev=~s/Y/G/g;
    $rev=~s/Z/C/g;
    $rev=~s/W/A/g;

    return $rev;
}

###############################################################################

sub printHelp{

    my $parnames=$_[0];
    my $parvalues=$_[1];
    
    print "\n";
    print "This script extracts ORF sequences. \n";
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

$parameters{"pathORFs"}="NA";
$parameters{"pathContigs"}="NA";
$parameters{"pathGeneticCode"}="NA";
$parameters{"pathOutputCDS"}="NA";
$parameters{"pathOutputProtein"}="NA";

my %defaultvalues;
my @defaultpars=("pathORFs","pathContigs", "pathGeneticCode",  "pathOutputCDS", "pathOutputProtein");

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

print "Reading de novo contigs...\n";

my %contigs;
readFasta($parameters{"pathContigs"}, \%contigs);

my $nb=keys %contigs;
print "Found ".$nb." contigs.\n";
print "Done.\n";

print "Reading genetic code...\n";
my %geneticcode;
readGeneticCode($parameters{"pathGeneticCode"},\%geneticcode);
print "Done.\n";

###################################################################################

print "Reading ORFs...\n";

my %orfs;
readORFs($parameters{"pathORFs"}, \%orfs);
    
print "Done.\n";

###################################################################################

print "Extracting ORFs and writing output...\n";

open(my $outputcds, ">".$parameters{"pathOutputCDS"});
open(my $outputprotein, ">".$parameters{"pathOutputProtein"});

my $nbaccepted=0;

foreach my $id (keys %orfs){

    if(exists $contigs{$id}){
	my $nborf=@{$orfs{$id}{"start"}};

	my %sequences;
       	
	for(my $i=0; $i<$nborf; $i++){
	    my $start=${$orfs{$id}{"start"}}[$i];
	    my $end=${$orfs{$id}{"end"}}[$i];
	    my $strand=${$orfs{$id}{"strand"}}[$i];

	    my $length=$end-$start+1;

	    if($length%3!=0){
		print "Weird! CDS length is not multiple of 3.\n";
		exit(1);
	    }

	    my $cds;
	    my $prot;

	    if($strand eq "+"){
		$cds=substr $contigs{$id}, ($start-1), ($end-$start+1);
		$prot=translateCDS($cds, \%geneticcode);
	    } else{
		my $seq=substr $contigs{$id}, ($start-1), ($end-$start+1);
		$cds=reverseComplement($seq);
		$prot=translateCDS($cds, \%geneticcode);
	    }
	    
	    
	    $nbaccepted++;
	    my $newid=$id.":".$start."-".$end.":".$strand;
	    
	    writeSequence($cds, $newid, $outputcds);
	    writeSequence($prot, $newid, $outputprotein);
	}

    }
    else{
	print "Weird!!! cannot find sequence for ".$id."\n";
    }
}


close($outputcds);
close($outputprotein);

print "Kept ".$nbaccepted." ORFs.\n";
print "Done.\n";

###################################################################################
###################################################################################
