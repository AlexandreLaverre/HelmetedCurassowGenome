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

##########################################################################

sub writeGTF{
    my $transcripts=$_[0];
    my $source=$_[1];
    my $output=$_[2];

    foreach my $tx (keys %{$transcripts}){
	my $gene=$transcripts->{$tx}{"gene"};
	my $chr=$transcripts->{$tx}{"chr"};
	my $strand=$transcripts->{$tx}{"strand"};

	my $txstart=min @{$transcripts->{$tx}{"start"}};
	my $txend=max @{$transcripts->{$tx}{"end"}};

	print $output $chr."\t".$source."\ttranscript\t".$txstart."\t".$txend."\t.\t".$strand."\t.\tgene_id \"".$gene."\"; transcript_id \"".$tx."\";\n";

	my $nbex=@{$transcripts->{$tx}{"start"}};
	
	for(my $i=0; $i<$nbex; $i++){
	    my $start=${$transcripts->{$tx}{"start"}}[$i];
	    my $end=${$transcripts->{$tx}{"end"}}[$i];

	    print $output $chr."\t".$source."\tCDS\t".$start."\t".$end."\t.\t".$strand."\t.\tgene_id \"".$gene."\"; transcript_id \"".$tx."\";\n";
	}
    }
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

sub computePropStop{
    my $seq=$_[0];

    my $lastchar=substr $seq, -1, 1;

    if($lastchar eq "*"){
	chop $seq;
    }
    
    my $nbX = ($seq =~ tr/X//);
    my $nbstar = ($seq =~ tr/*//);

    my $nbtot=length $seq;
    my $propstop="NA";
    
    if($nbtot>0){
        $propstop=($nbX+$nbstar)/$nbtot;
    } else{
	print "Weird! sequence is empty: ". $seq."\n";
    }

    return $propstop;
}

#################################################################################

sub printHelp{

    my $parnames=$_[0];
    my $parvalues=$_[1];
    
    print "\n";
    print "This script filters annotations.\n";
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
$parameters{"pathProteins"}="NA";
$parameters{"minProteinLength"}="NA";
$parameters{"source"}="NA";
$parameters{"pathOutputGTF"}="NA";
$parameters{"pathOutputFasta"}="NA";

my @defaultpars=("pathAnnotGTF", "pathProteins", "minProteinLength",  "source",  "pathOutputGTF", "pathOutputFasta");


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

print "Reading proteins...\n";

my %proteins;
readFasta($parameters{"pathProteins"}, \%proteins);

print "Done.\n";

##############################################################

print "Filtering transcripts...\n";

my @alltranscripts=keys %transcripts;

my $minlen=$parameters{"minProteinLength"}+0;

my $nbtooshort=0;
my $nbwithstop=0;

foreach my $tx (@alltranscripts){
    if(!exists $proteins{$tx}){
	print "Weird! cannot find protein for ".$tx."\n";
	exit(1);
    }

    my $seq=$proteins{$tx};
    my $len=length $seq;

    if($len<$minlen){
	$nbtooshort++;
	delete $transcripts{$tx};
    } else{
	my $propstop=computePropStop($seq);

	if($propstop>0){
	    $nbwithstop++;
	    delete $transcripts{$tx};
	}
    }
}

print "Deleted ".$nbtooshort." transcripts that were too short.\n";
print "Deleted ".$nbwithstop." transcripts with stop/X codons.\n";

print "Done.\n";

##############################################################

print "Writing GTF output...\n";

open(my $output, ">".$parameters{"pathOutputGTF"});

writeGTF(\%transcripts, $parameters{"source"}, $output);
 
close($output);
print "Done.\n";

##############################################################

print "Writing fasta output...\n";

open(my $output, ">".$parameters{"pathOutputFasta"});

foreach my $tx (keys %transcripts){
    my $gene=$transcripts{$tx}{"gene"};
    my $name=$tx." gene:".$gene." transcript:".$tx;
    my $seq=$proteins{$tx};
    
    writeSequence($seq, $name, $output);
}
 
close($output);

print "Done.\n";

##############################################################
