use List::Util qw(first max maxstr min minstr reduce shuffle sum);
use strict;

##############################################################
##############################################################

sub readTranscripts{
    my $pathin=$_[0];
    my $transcripts=$_[1];

    open(my $input, $pathin);

    my $line=<$input>;

    while($line){
	chomp $line;
	my @s=split("\t", $line);

	my $type=$s[2];

	if($type eq "transcript"){
	    my $info=$s[8];
	    my @infoarray=split(";", $info);

	    my $txid=findInfo("transcript_id", \@infoarray);
	    my $geneid=findInfo("gene_id", \@infoarray);

	    if($txid eq "NA" || $geneid eq "NA"){
		print "Cannot find info in line ".$line."\n";
		exit(1);
	    }

	    if(exists $transcripts->{$geneid}){
		$transcripts->{$geneid}{$txid}=1;
	    } else{
		$transcripts->{$geneid}={$txid=>1}
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

##############################################################

sub renameTranscripts{
    my $transcripts=$_[0];
    my $newnames=$_[1];

    foreach my $gene (keys %{$transcripts}){
	my @alltx=keys %{$transcripts->{$gene}};

	$newnames->{$gene}={};

	for(my $i=0; $i<@alltx; $i++){
	    my $txid=$alltx[$i];
	    my $newname=$gene."_mrna_".($i+1);

	    $newnames->{$gene}{$txid}=$newname;
	}
    }
}

##############################################################

sub printHelp{

    my $parnames=$_[0];
    my $parvalues=$_[1];

    print "\n";
    print "This script renames GeMoMa annotations.\n";
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
$parameters{"pathInputGTF"}="NA";
$parameters{"pathOutputGTF"}="NA";

my @defaultpars=("pathInputGTF", "pathOutputGTF");


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

print "Reading initial transcript names...\n";

my %transcripts;
readTranscripts($parameters{"pathInputGTF"}, \%transcripts);

my $nbg=keys %transcripts;
print "Saw ".$nbg." genes.\n";

print "Done.\n";

##############################################################

print "Renaming transcripts...\n";

my %newnames;
renameTranscripts(\%transcripts, \%newnames);

print "Done.\n";

##############################################################

print "Writing output...\n";

open(my $input, $parameters{"pathInputGTF"});
open(my $output, ">".$parameters{"pathOutputGTF"});

my $line=<$input>;

while($line){
    my $prefix=substr $line, 0,1;

    if($prefix eq "#"){
	print $output $line;
    } else{
	chomp $line;
	my @s=split("\t", $line);

	my $info=$s[8];
	my @infoarray=split(";", $info);

	my $txid=findInfo("transcript_id", \@infoarray);
	my $geneid=findInfo("gene_id", \@infoarray);

	if($txid eq "NA" || $geneid eq "NA"){
	    print "Cannot find info in line ".$line."\n";
	    exit(1);
	}

	if(exists $newnames{$geneid}{$txid}){
	    my $nn=$newnames{$geneid}{$txid};

	    print $output join("\t", @s[0..7])."\tgene_id \"".$geneid."\"; transcript_id \"".$nn."\";\n"
	} else{
	    print "Weird, cannot find new name for ".$geneid."\t".$txid."\n";
	}

    }

    $line=<$input>;
}

close($input);
close($output);

print "Done.\n";

##############################################################
