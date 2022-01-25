#!/bin/bash

export sp=$1
export assembly=$2
export sample=$3
export cluster=$4
export nthreads=$5

#############################################################################

if [ ${cluster} = "cloud" ]; then
    export path=/ifb/data/mydatalocal/HelmetedCurassowGenome
fi

export pathRNASeq=${path}/data/RNASeq/${sp}
export pathResults=${path}/results/RNASeq_alignments/${sp}/${assembly}
export pathScripts=${path}/scripts/RNASeq_alignments

## hisat2-2.2.1

#############################################################################

if [ ${assembly} = "MEGAHIT_RAGOUT" ]; then
    export prefix="genome_sequence_renamed"
fi

#############################################################################

export pathIndex=${pathResults}/${prefix}

#############################################################################

if [ -e ${pathResults}/${sample} ]; then
    echo "output dir exists"
else
    mkdir -p ${pathResults}/${sample}
fi

#############################################################################

if [ ${sp} = "Basiliscus_vittatus" ]; then
    export strand="--rna-strandness R" ## TruSeq
else
  export strand=""
fi

#############################################################################

if [ -e ${pathResults}/${sample}/accepted_hits.sam.gz ]||[ -e ${pathResults}/${sample}/accepted_hits.bam ]||[ -e ${pathResults}/${sample}/accepted_hits.sam ]||[ -e ${pathResults}/${sample}/accepted_hits_clean.sam ]||[ -e ${pathResults}/${sample}/accepted_hits_clean.sam.gz ]; then
    echo "already done (or ongoing)"
else

    export pathR=""
    export pathR1=""
    export pathR2=""
    
    export nbR1=0
    export nbR2=0
   
    if [ -e ${pathRNASeq}/${sample}.fastq.gz ]; then
	export pathR=${pathRNASeq}/${sample}.fastq.gz,${pathR}
	export type="single_end"
    fi
    
    if [ -e ${pathRNASeq}/${sample}_R1.fastq.gz ]; then
	export pathR1=${pathRNASeq}/${sample}_R1.fastq.gz,${pathR1}
	export nbR1=$((nbR1+1))
	export type="paired_end"
    fi

     if [ -e ${pathRNASeq}/${sample}_R2.fastq.gz ]; then
	export pathR2=${pathRNASeq}/${sample}_R2.fastq.gz,${pathR2}
	export nbR2=$((nbR2+1))
	export type="paired_end"
    fi    
    
    #############################################################################
    
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
	echo "hisat2 --seed 19 -p ${nthreads} -x ${pathIndex} -U ${pathR} -S ${pathResults}/${sample}/accepted_hits.sam ${strand} --max-intronlen 1000000 --dta-cufflinks --no-unal --met-file ${pathResults}/${sample}/metrics.txt  --novel-splicesite-outfile ${pathResults}/${sample}/novel_splicesites.txt >& ${pathResults}/${sample}/align_summary.txt">> ${pathScripts}/bsub_script_hisat
    fi
        
    if [ ${type} = ${paired_end} ]; then
	echo "paired-end"
	echo "hisat2 --seed 19 -p ${nthreads} -x ${pathIndex} -1 ${pathR1} -2 ${pathR2} -S ${pathResults}/${sample}/accepted_hits.sam ${strand} --max-intronlen 1000000 --dta-cufflinks --no-unal --met-file ${pathResults}/${sample}/metrics.txt  --novel-splicesite-outfile ${pathResults}/${sample}/novel_splicesites.txt >& ${pathResults}/${sample}/align_summary.txt">> ${pathScripts}/bsub_script_hisat
    fi
    
    if [ ${cluster} = "cloud" ]; then
	chmod a+x ${pathScripts}/bsub_script_hisat
	${pathScripts}/bsub_script_hisat
    fi
    
    if [ ${cluster} = "pbil" ]; then
	sbatch ${pathScripts}/bsub_script_hisat
    fi
    
fi

#############################################################################
