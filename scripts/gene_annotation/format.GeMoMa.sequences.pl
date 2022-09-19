use List::Util qw(first max maxstr min minstr reduce shuffle sum);
use strict;

##############################################################
##############################################################

sub readGTF{
    my $pathin=$_[0];
    my $transcripts=$_[1];

    open(my $input, $pathin);

    my $line=<$input>;

    while($line){
	chomp $line;
	my @s=split("\t", $line);

	my $type=$s[2];

	if($type eq "CDS"){
	    my $chr=$s[0];
	    my $start=$s[3]+0;
	    my $end=$s[4]+0;
	    my $strand=$s[6];

	    my $info=$s[8];
	    my @infoarray=split(";", $info);

	    my $txid=findInfo("transcript_id", \@infoarray);
	    my $geneid=findInfo("gene_id", \@infoarray);

	    if($txid eq "NA" || $geneid eq "NA"){
		print "Cannot find info in line ".$line."\n";
		exit(1);
	    }

	    if(exists $transcripts->{$txid}){
		push(@{$transcripts->{$txid}{"start"}}, $start);
		push(@{$transcripts->{$txid}{"end"}}, $end);
	    } else{
		$transcripts->{$txid}={"start"=>[$start], "end"=>[$end], "chr"=>$chr, "strand"=>$strand, "gene"=>$geneid};
	    }
	}

	$line=<$input>;
    }

    close($input);
}

##############################################################

sub findInfo{
    my $pattern=$_[0];
    my $array=$_[1];

    my $res="NA";

    my @grepres=grep(/${pattern}/,@{$array});

    if(@grepres==1){
	my @t=split("\"",$grepres[0]);
	$res=$t[1];
    }

    return $res;
}

###################################################################################

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

	    # print "saw chromosome ".$id."\n";

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

##################################################################################

sub printHelp{

    my $parnames=$_[0];
    my $parvalues=$_[1];

    print "\n";
    print "This script formats sequences to display transcript and gene name.\n";
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
$parameters{"pathAnnotGTF"}="NA";
$parameters{"pathSequences"}="NA";
$parameters{"pathOutput"}="NA";

my @defaultpars=("pathAnnotGTF", "pathSequences",  "pathOutput");


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

print "Reading annotations...\n";

my %transcripts;
readGTF($parameters{"pathAnnotGTF"}, \%transcripts);

my $nbtx=keys %transcripts;

print "Found ".$nbtx." transcripts in GTF.\n";

print "Done.\n";

##############################################################

print "Reading sequences...\n";

my %sequences;
readFasta($parameters{"pathSequences"}, \%sequences);

print "Done.\n";

##############################################################

print "Writing fasta output...\n";

open(my $output, ">".$parameters{"pathOutput"});

foreach my $tx (keys %transcripts){
    my $gene=$transcripts{$tx}{"gene"};
    my $name=$tx." gene:".$gene." transcript:".$tx;
    my $seq=$sequences{$tx};

    writeSequence($seq, $name, $output);
}

close($output);

print "Done.\n";

##############################################################
