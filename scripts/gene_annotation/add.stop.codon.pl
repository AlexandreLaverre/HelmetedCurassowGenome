use List::Util qw(first max maxstr min minstr reduce shuffle sum);
use strict;

##############################################################

sub readLastExons{
    my $pathin=$_[0];
    my $lastexon=$_[1];

    my $input;
    my @s=split("\\.",$pathin);
    my $ext=$s[-1];

    if($ext eq "gz"){
	open($input, "zcat $pathin |");
    } else{
	open($input, $pathin);
    }
    
    my $line=<$input>;
   
    
    while($line){
	$line=~s/\|/-/g;
	chomp $line;
	my @s=split("\t", $line);
	
	my $type=$s[2];

	if($type eq "CDS"){
	    my $chr=$s[0];
	    my $start=$s[3]+0; ## 1-based
	    my $end=$s[4]+0;
	    my $strand=$s[6];
	   
	    my $info=$s[8];
	    my @t=split(";", $info);
	    my $txid=findInfo("orig_transcript_id", \@t);

	    if($txid eq "NA"){
		print "Weird! cannot find transcript id in line ".$line."\n";
		exit(1);
	    }
	    
	    if(exists $lastexon->{$txid}){
		if($strand eq "+" || $strand eq "1"){
		    if($start > $lastexon->{$txid}){
			$lastexon->{$txid}=$start;
		    }
		} else{
		    if($strand eq "-" || $strand eq "-1"){
			if($start < $lastexon->{$txid}){
			    $lastexon->{$txid}=$start;
			}
		    } else{
			print "Weird strand! ".$strand."\n";
			exit(1);
		    }
		}
	    } else{
		$lastexon->{$txid}=$start;
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
   
    my $nbg=@grepres;
    my $nbreal=0;

    if(@grepres==1){
	my @t=split("=",$grepres[0]);
	$res=$t[1];
	return $res;
    } else{ 
	return "NA";
    }
}

##############################################################

sub readSequenceSizes{
    my $pathin=$_[0];
    my $sizes=$_[1];
    
    my $input;
    my @s=split("\\.",$pathin);
    my $ext=$s[-1];
    
    if($ext eq "gz"){
	open($input, "zcat $pathin |");
    } else{
	open($input, $pathin);
    }
    
    my $line=<$input>;
    
    while($line){
	$line=~s/\|/-/g;
	my $prefix=substr $line, 0, 17;
	
	if($prefix eq "##sequence-region"){
	    chomp $line;
	    my @s=split(" ", $line);
	    my $chr=$s[1];
	    my $size=$s[3]+0;

	    $sizes->{$chr}=$size;
	}
	
	$line=<$input>;
    }
    
    close($input);
}

##############################################################

sub printHelp{

    my $parnames=$_[0];
    my $parvalues=$_[1];
    
    print "\n";
    print "This script adds stop codon to CDS coordinates. \n";
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

$parameters{"pathInputGFF"}="NA";
$parameters{"pathOutputGFF"}="NA";

my %defaultvalues;
my @defaultpars=("pathInputGFF","pathOutputGFF");


my %numericpars;
my @numericpars=();

foreach my $par (@numericpars){
    $numericpars{$par}=1;
}

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

#####################################################################
#####################################################################

print "Reading CDS last exon coordinates...\n";

my %lastexon;
readLastExons($parameters{"pathInputGFF"}, \%lastexon);

my $nbtx=keys %lastexon;

print "Found last exon coordinates for ".$nbtx." transcripts.\n";

print "Done.\n";

#####################################################################

print "Reading chromosome sizes...\n";
my %chrsizes;

readSequenceSizes($parameters{"pathInputGFF"}, \%chrsizes);

print "Done.\n";

#####################################################################

print "Reading input and writing updated GFF...\n";

my $pathin=$parameters{"pathInputGFF"};
my $input;
my @s=split("\\.",$pathin);
my $ext=$s[-1];

if($ext eq "gz"){
    open($input, "zcat $pathin |");
} else{
    open($input, $pathin);
}

open(my $output, ">".$parameters{"pathOutputGFF"});

my $line=<$input>;

while($line){
    
    chomp $line;
    ## we replace "|" in transcript names
    $line=~s/\|/-/g;
    
    my @s=split("\t", $line);
    
    my $type=$s[2];
    
    if($type eq "CDS"){
	my $chr=$s[0]; ## 1-based

	if(!exists $chrsizes{$chr}){
	    print "Weird! cannot find chr size for ".$chr."\n";
	    exit(1);
	}
	
	my $size=$chrsizes{$chr};
	
	my $start=$s[3]+0; ## 1-based
	my $end=$s[4]+0; ## 1-based
	my $strand=$s[6];
	   
	my $info=$s[8];
	my @t=split(";", $info);
	my $txid=findInfo("orig_transcript_id", \@t);

	if(!exists $lastexon{$txid}){
	    print "Weird! no info for ".$txid."\n";
	    exit(1);
	}
		
	if($lastexon{$txid}==$start){
	    my $lens=@s;
	    my $newline=join("\t", @s[0..2]);

	    if($strand eq "+" || $strand eq "1"){
		## we modify the end coordinate
		if(($end+3)<=$size){
		    $newline.="\t".$start."\t".($end+3)."\t".join("\t", @s[5..($lens-1)]);
		} else{
		    $newline.="\t".$start."\t".$end."\t".join("\t", @s[5..($lens-1)]);
		}
	    } else{
		if($strand eq "-" || $strand eq "-1"){
		    ## we modify the start coordinate
		    if($start>3){
			$newline.="\t".($start-3)."\t".$end."\t".join("\t", @s[5..($lens-1)]);
		    } else{
			$newline.="\t".$start."\t".$end."\t".join("\t", @s[5..($lens-1)]);
		    }
		} else{
		    print "Weird strand for ".$txid."\n";
		    exit(1);
		}
	    }
	    print $output $newline."\n";
	    
	} else{
	    print $output $line."\n";
	}
    } else{
	print $output $line."\n";
    }
    
    $line=<$input>;
}

close($input);
close($output);

print "Done.\n";

#####################################################################
