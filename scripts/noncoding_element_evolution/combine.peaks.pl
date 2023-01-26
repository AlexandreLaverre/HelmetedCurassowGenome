use List::Util qw(first max maxstr min minstr reduce shuffle sum);
use strict;

#######################################################################

sub readSampleInfo{
    my $path=$_[0];
    my $sampleinfo=$_[1];

    open(my $input, $path);

    my %header;

    my $line=<$input>;
    chomp $line;
    my @s=split("\t", $line);
    for(my $i=0; $i<@s; $i++){
	$header{$s[$i]}=$i;
    }
    
    $line=<$input>;

    while($line){
	chomp $line;
	my @s=split("\t", $line);

	my $acc=$s[$header{"GEO_Accession (exp)"}];
	my $tissue=$s[$header{"Tissue"}];
	my $age=$s[$header{"Age"}];

	$sampleinfo->{$acc}={"tissue"=>$tissue, "age"=>$age};
	
	$line=<$input>;
    }
    
    close($input);
}

#######################################################################

## BED format: chr, start, end, other info

sub readCoordinates{
    my $paths=$_[0];
    my $samples=$_[1];
    my $coords=$_[2];

    ## paths and samples

    my @pathlist=split(",", $paths);
    my @samplelist=split(",", $samples);

    my $nbp=@pathlist;
    my $nbs=@samplelist;

    if($nbp!=$nbs){
	print "Weird! saw ".$nbs." samples and ".$nbp." paths, something is wrong.\n";
	exit(1);
    }

    my %hashcoords;

    for(my $i=0; $i<$nbp; $i++){
	my $pathin=$pathlist[$i];
	my $sampleid=$samplelist[$i];

	print "Reading data from ".$pathin." for ".$sampleid."\n";

	open(my $input, $pathin);
	my $line=<$input>;

	while($line){
	    chomp $line;
	    my @s=split("\t", $line);

	    my $chr=$s[0];
	    my $start=$s[1]+0; ## 0-based coordinates
	    my $end=$s[2]+0;
	    my $id=$sampleid.":".$chr.":".$start.":".$end;

	    if(exists $hashcoords{$chr}){
		if(exists $hashcoords{$chr}{$start}){
		    push(@{$hashcoords{$chr}{$start}{"end"}}, $end);
		    push(@{$hashcoords{$chr}{$start}{"id"}}, $id);
		} else{
		    $hashcoords{$chr}{$start}={"end"=>[$end], "id"=>[$id]};
		}
	    } else{
		$hashcoords{$chr}={$start=>{"end"=>[$end], "id"=>[$id]}};
	    }

	    $line=<$input>;
	}

	close($input);
    }

    print "ordering coordinates...\n";

    foreach my $chr (keys %hashcoords){
	$coords->{$chr}={"start"=>[], "end"=>[], "id"=>[]};

	my @allstart=keys %{$hashcoords{$chr}};
	my @orderedstart=sort {$a<=>$b} @allstart;

	foreach my $start (@orderedstart){
	    my $nb=@{$hashcoords{$chr}{$start}{"end"}};

	    for(my $i=0; $i<$nb; $i++){
		push(@{$coords->{$chr}{"start"}}, $start);
		push(@{$coords->{$chr}{"end"}}, ${$hashcoords{$chr}{$start}{"end"}}[$i]);
		push(@{$coords->{$chr}{"id"}}, ${$hashcoords{$chr}{$start}{"id"}}[$i]);
	    }
	}
    }
}

#######################################################################

sub mergeCoordinates{
    my $coords=$_[0];
    my $merged=$_[1];

    foreach my $chr (keys %{$coords}){
	$merged->{$chr}={"start"=>[], "end"=>[], "id"=>[]};

	my $nbcoords=@{$coords->{$chr}{"start"}};

	my $currentstart=${$coords->{$chr}{"start"}}[0];
	my $currentend=${$coords->{$chr}{"end"}}[0];

	my $index=0;

	for(my $i=1; $i<$nbcoords; $i++){
	    my $thisstart=${$coords->{$chr}{"start"}}[$i];
	    my $thisend=${$coords->{$chr}{"end"}}[$i];

	    if($thisstart>=$currentstart && $thisstart<=($currentend+1)){
		if($thisend>$currentend){
		    $currentend=$thisend;
		}
	    } else{
		push(@{$merged->{$chr}{"start"}}, $currentstart);
		push(@{$merged->{$chr}{"end"}}, $currentend);
		push(@{$merged->{$chr}{"id"}}, "merged_".$chr."_".$index);

		$currentstart=$thisstart;
		$currentend=$thisend;
		$index++;
	    }
	}

	## don't forget last block

	push(@{$merged->{$chr}{"start"}}, $currentstart);
	push(@{$merged->{$chr}{"end"}}, $currentend);
	push(@{$merged->{$chr}{"id"}}, "merged_".$chr."_".$index);
    }
}

#######################################################################

sub extractOverlap{
    my $coords1=$_[0];
    my $coords2=$_[1];
    my $overlap=$_[2];

    foreach my $chr (keys %{$coords1}){
	if(exists $coords2->{$chr}){
	    my $nb1=@{$coords1->{$chr}{"start"}};
	    my $nb2=@{$coords2->{$chr}{"start"}};

	    my $firstj=0;

	    for(my $i=0; $i<$nb1; $i++){
		my $start1=${$coords1->{$chr}{"start"}}[$i];
		my $end1=${$coords1->{$chr}{"end"}}[$i];
		my $id1=${$coords1->{$chr}{"id"}}[$i];

		my $j=$firstj;

		while($j<$nb2 && ${$coords2->{$chr}{"end"}}[$j]<$start1){
		    $j++;
		}

		$firstj=$j;

		while($j<$nb2 && ${$coords2->{$chr}{"start"}}[$j]<=$end1){
		    my $start2=${$coords2->{$chr}{"start"}}[$j];
		    my $end2=${$coords2->{$chr}{"end"}}[$j];
		    my $id2=${$coords2->{$chr}{"id"}}[$j];

		    my $M=max($start1, $start2);
		    my $m=min($end1, $end2);

		    if($M<=$m){
			if(exists $overlap->{$id1}){
			    push(@{$overlap->{$id1}}, $id2);
			} else{
			    $overlap->{$id1}=[$id2];
			}
		    }

		    $j++;
		}
	    }

	}
    }
}


##############################################################

sub printHelp{

    my $parnames=$_[0];
    my $parvalues=$_[1];

    print "\n";
    print "This script combines peaks across samples.\n";
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
$parameters{"pathSampleInfo"}="NA";
$parameters{"samples"}="NA";
$parameters{"pathsCoordinates"}="NA";
$parameters{"pathOutput"}="NA";

my %defaultvalues;
my @defaultpars=("pathSampleInfo", "samples", "pathsCoordinates", "pathOutput");

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

print "Reading sample info...\n";

my %sampleinfo;
readSampleInfo($parameters{"pathSampleInfo"}, \%sampleinfo);

print "Done.\n";

######################################################################

print "Reading peak coordinates...\n";

my %coords;

readCoordinates($parameters{"pathsCoordinates"}, $parameters{"samples"}, \%coords);

print "Done.\n";

######################################################################

print "Merging coordinates...\n";

my %merged;
mergeCoordinates(\%coords, \%merged);

print "Done.\n";

######################################################################

print "Computing overlap...\n";

my %overlap;

extractOverlap(\%merged, \%coords, \%overlap);

print "Done.\n";

######################################################################

print "Writing output...\n";

open(my $output, ">".$parameters{"pathOutput"});

print $output "ID\tChr\tStart\tEnd\tSamples\tTissues\tTissAge\n";

foreach my $chr (keys %merged){
    my $nb=@{$merged{$chr}{"start"}};

    for(my $i=0; $i<$nb; $i++){
	my $id=${$merged{$chr}{"id"}}[$i];
	my $start=${$merged{$chr}{"start"}}[$i];
	my $end=${$merged{$chr}{"end"}}[$i];

	if(exists $overlap{$id}){
	    my %samples;
	    my %tissues;
	    my %tissages;
	    
	    foreach my $otherid (@{$overlap{$id}}){
		my @s=split(":",$otherid);
		my $sample=$s[0];
		$samples{$sample}=1;

		if(!exists $sampleinfo{$sample}){
		    print "Weird! cannot find info for ".$sample."\n";
		    exit(1);
		} else{
		    my $tiss=$sampleinfo{$sample}{"tissue"};
		    my $age=$sampleinfo{$sample}{"age"};
		    $tissues{$tiss}=1;
		    $tissages{$tiss.":".$age}=1;
		}
	    }

	    my $samplelist=join(",", sort (keys %samples));
	    my $tisslist=join(",", sort (keys %tissues));
	    my $tissagelist=join(",", sort (keys %tissages));
	 

	    print $output $id."\t".$chr."\t".$start."\t".$end."\t".$samplelist."\t".$tisslist."\t".$tissagelist."\n";

	} else{
	    print "Weird! ".$id." doesn't overlap with anything!\n";
	    exit(1);
	}
    }
}

close($output);

print "Done.\n";

######################################################################
