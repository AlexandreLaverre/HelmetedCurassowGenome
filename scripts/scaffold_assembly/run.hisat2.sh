#!/bin/bash

export species=$1
export sample=$2
export cluster=$3
export nthreads=$4

#############################################################################

if [ ${cluster} = "pbil" ]; then
    export path=/beegfs/data/necsulea/HelmetedCurassowGenome
fi

if [ ${cluster} = "cloud" ]; then
    export path=/home/ubuntu/data/mydatalocal/HelmetedCurassowGenome
fi

export pathRNASeq=${path}/data/RNASeq/${species}
export pathDocs=${path}/docs
export pathResults=${path}/results/genome_assembly/${species}/MEGAHIT
export pathScripts=${path}/scripts/scaffold_assembly

export pathIndex=${pathResults}/final.contigs

#############################################################################

export pathR1=`grep ${species} ${pathDocs}/Samples.txt | grep RNASeq | grep $'\t'${sample}$'\t' | cut -f 7`
export pathR2=`grep ${species} ${pathDocs}/Samples.txt | grep RNASeq | grep $'\t'${sample}$'\t' | cut -f 8`

export type="NA"

if [ ${pathR2} = "." ]; then
    export type="single_end"
else
    if [ -e ${pathRNASeq}/${pathR2} ]; then
	export type="paired_end"
    else
	echo "cannot find R2"
	exit
    fi
fi

#############################################################################

export strand=""

if [ ${species} = "Basiliscus_vittatus" ]; then
    ## TruSeq Stranded
    export strand="--rna-strandness R"
fi

#############################################################################

if [ -e ${pathResults}/hisat2_${sample}/accepted_hits.sam.gz ]||[ -e ${pathResults}/hisat2_${sample}/accepted_hits.bam ]||[ -e ${pathResults}/hisat2_${sample}/accepted_hits.sam ]; then
    echo "already done (or ongoing)"
else
    
    if [ -e ${pathResults}/hisat2_${sample} ]; then
    	echo "dir output already there"
    else
    	mkdir ${pathResults}/hisat2_${sample}
    fi
      
    echo "#!/bin/bash" > ${pathScripts}/bsub_script_hisat

    if [ ${cluster} = "pbil" ]; then
	echo "#!/bin/bash " > ${pathScripts}/bsub_script_hisat
	echo "#SBATCH --job-name=hisat_${sample}" >>  ${pathScripts}/bsub_script_hisat
	echo "#SBATCH --partition=normal" >>  ${pathScripts}/bsub_script_hisat
	echo "#SBATCH --output=${pathScripts}/std_out_${sample}" >>  ${pathScripts}/bsub_script_hisat
	echo "#SBATCH --error=${pathScripts}/std_err_${sample}" >>  ${pathScripts}/bsub_script_hisat
	echo "#SBATCH --cpus-per-task=${nthreads}" >>  ${pathScripts}/bsub_script_hisat ## ${nthreads} CPU
	echo "#SBATCH --time=24:00:00" >>  ${pathScripts}/bsub_script_hisat ## 24 hours
	echo "#SBATCH --mem=4G" >>  ${pathScripts}/bsub_script_hisat ## 5g per CPU
    fi	

    if [ ${type} = "single_end" ]; then
	echo "single-end"
	echo "hisat2 --seed 19 -p ${nthreads} -x ${pathIndex} -U ${pathRNASeq}/${pathR1} -S ${pathResults}/hisat2_${sample}/accepted_hits.sam ${strand} --max-intronlen 1000000 --dta-cufflinks --no-unal --met-file ${pathResults}/hisat2_${sample}/metrics.txt  --novel-splicesite-outfile ${pathResults}/hisat2_${sample}/novel_splicesites.txt >& ${pathResults}/hisat2_${sample}/align_summary.txt">> ${pathScripts}/bsub_script_hisat
    else
	if [ ${type} = "paired_end" ]; then
	    echo "paired-end"
	    echo "hisat2 --seed 19 -p ${nthreads} -x ${pathIndex} -1 ${pathRNASeq}/${pathR1} -2 ${pathRNASeq}/${pathR2} -S ${pathResults}/hisat2_${sample}/accepted_hits.sam ${strand} --max-intronlen 1000000 --dta-cufflinks --no-unal --met-file ${pathResults}/hisat2_${sample}/metrics.txt  --novel-splicesite-outfile ${pathResults}/hisat2_${sample}/novel_splicesites.txt >& ${pathResults}/hisat2_${sample}/align_summary.txt">> ${pathScripts}/bsub_script_hisat
	fi
    fi

    echo "samtools sort -@ ${nthreads} -O bam -o ${pathResults}/hisat2_${sample}/accepted_hits.bam ${pathResults}/hisat2_${sample}/accepted_hits.sam " >> ${pathScripts}/bsub_script_hisat
    echo "gzip ${pathResults}/hisat2_${sample}/accepted_hits.sam " >> ${pathScripts}/bsub_script_hisat
    
    if [ ${cluster} = "pbil" ]; then
	sbatch ${pathScripts}/bsub_script_hisat
    fi
    
    if [ ${cluster} = "cloud" ]; then
	chmod a+x ${pathScripts}/bsub_script_hisat
	${pathScripts}/bsub_script_hisat
    fi
fi

#############################################################################
