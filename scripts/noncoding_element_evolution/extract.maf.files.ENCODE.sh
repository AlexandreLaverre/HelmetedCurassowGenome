#!/bin/bash

set -e 

export cluster=$1
export forceMafsInRegion=$2

#########################################################################

if [ ${cluster} = "cloud" ]; then
    export path=/ifb/data/mydatalocal/HelmetedCurassowGenome
fi

export pathHAL=${path}/results/whole_genome_alignments/avian_366
export pathResults=${path}/results/noncoding_element_evolution/ENCODE_ATAC-seq/Mouse
export pathScripts=${path}/scripts/noncoding_element_evolution

#########################################################################

export targetGenomes="Struthio_camelus,Dromaius_novaehollandiae,Gallus_gallus,Meleagris_gallopavo,Penelope_pileata,Alectura_lathami,Anas_platyrhynchos_platyrhynchos,Grus_americana,Calidris_pugnax,Upupa_epops,Rhinopomastus_cyanomelas,Strix_occidentalis,Ficedula_albicollis,Aquila_chrysaetos,Geospiza_fortis,Parus_major,Serinus_canaria,Coturnix_japonica,Casuarius_casuarius,Numida_meleagris,Anseranas_semipalmata,Anser_cygnoid,Balearica_regulorum,Bucorvus_abyssinicus,Buceros_rhinoceros,Pauxi_pauxi"

export refGenome="Gallus_gallus"

#########################################################################

if [ -e ${pathResults}/combined_peaks_galGal6_formatted.maf ]; then
    echo "already extracted MAF"
else
    docker run -v ${path}:/ifb/data/mydatalocal/HelmetedCurassowGenome --rm -t quay.io/comparative-genomics-toolkit/cactus:v1.3.0 hal2maf ${pathHAL}/366-avian.hal ${pathResults}/combined_peaks_galGal6_formatted.maf  --targetGenomes ${targetGenomes} --refGenome ${refGenome}  --noDupes --refTargets ${pathResults}/combined_peaks_galGal6_formatted.bed
fi

#########################################################################

if [ -e ${pathResults}/combined_peaks_galGal6_formatted_ordered.maf ]; then
    echo "already ordered MAF"
else
    perl ${pathScripts}/order.maf.files.pl --pathMAFInput=${pathResults}/combined_peaks_galGal6_formatted.maf --refSpecies=Gallus_gallus --pathMAFOutput=${pathResults}/combined_peaks_galGal6_formatted_ordered.maf
fi

#########################################################################

if [ -e ${pathResults}/mafs_by_element ]; then
    export nbaln=`ls ${pathResults}/mafs_by_element | grep maf | wc -l | cut -f 1 -d ' '`
    echo ${nbaln} "alignments extracted"
    export nbel=`wc -l ${pathResults}/combined_peaks_galGal6_formatted.bed | cut -f 1 -d ' ' `
    echo ${nbel} "elements originally"

    if [ ${forceMafsInRegion} = "true" ]; then
	mafsInRegion -outDir ${pathResults}/combined_peaks_galGal6_formatted.bed ${pathResults}/mafs_by_element ${pathResults}/combined_peaks_galGal6_formatted_ordered.maf
    else
	echo "not running mafsInRegion"
    fi
else
    mkdir -p ${pathResults}/mafs_by_element
    
    mafsInRegion -outDir ${pathResults}/combined_peaks_galGal6_formatted.bed ${pathResults}/mafs_by_element ${pathResults}/combined_peaks_galGal6_formatted_ordered.maf 
fi

#########################################################################

if [ -e ${pathResults}/mafs_by_element ]; then
    for file in `ls ${pathResults}/mafs_by_element | grep maf$`
    do
	export prefix=`basename ${file} .maf`
	
	if [ -e ${pathResults}/mafs_by_element/${prefix}.filtered.phy ]; then
	    echo "phylip format for "${prefix} "already done"
	else
	    if [ -e ${pathResults}/mafs_by_element/${prefix}.discarded ]; then
		echo ${prefix} "was discarded after filtering"
	    else

		export nbl=`wc -l ${pathResults}/mafs_by_element/${prefix}.maf | cut -f 1 -d ' '`

		if [ ${nbl} -gt 1 ]; then
		    echo ${prefix}
		    
		    msa_view ${pathResults}/mafs_by_element/${prefix}.maf --out-format FASTA  --missing-as-indels --clean-indels 3 --unmask > ${pathResults}/mafs_by_element/${prefix}.fa
		    
		    perl ${pathScripts}/filter.alignments.pl --minSpecies=5 --maxGapProportion=0.75 --minUngappedLength=50 --pathFastaInput=${pathResults}/mafs_by_element/${prefix}.fa --pathFastaOutput=${pathResults}/mafs_by_element/${prefix}.filtered.fa --pathOutputDiscarded=${pathResults}/mafs_by_element/${prefix}.discarded
		    
		    if [ -e ${pathResults}/mafs_by_element/${prefix}.filtered.fa ]; then
			msa_view ${pathResults}/mafs_by_element/${prefix}.filtered.fa --in-format FASTA --out-format PHYLIP > ${pathResults}/mafs_by_element/${prefix}.filtered.phy
		    fi
		else
		    echo "empty maf file" > ${pathResults}/mafs_by_element/${prefix}.discarded
		fi 
	    fi
	fi
    done
fi

#########################################################################
