#!/bin/bash

export cluster=$1

#########################################################################

if [ ${cluster} = "cloud" ]; then
    export path=/ifb/data/mydatalocal/HelmetedCurassowGenome
fi

export pathHAL=${path}/results/whole_genome_alignments/avian_366
export pathResults=${path}/results/noncoding_element_evolution/ENCODE_ATAC-seq/Mouse

#########################################################################

export targetGenomes="Struthio_camelus,Dromaius_novaehollandiae,Gallus_gallus,Meleagris_gallopavo,Penelope_pileata,Alectura_lathami,Anas_platyrhynchos_platyrhynchos,Grus_americana,Calidris_pugnax,Upupa_epops,Rhinopomastus_cyanomelas,Strix_occidentalis,Ficedula_albicollis,Aquila_chrysaetos,Geospiza_fortis,Parus_major,Serinus_canaria,Coturnix_japonica,Casuarius_casuarius,Numida_meleagris,Anseranas_semipalmata,Anser_cygnoid,Balearica_regulorum,Bucorvus_abyssinicus,Buceros_rhinoceros,Pauxi_pauxi"

export refGenome="Gallus_gallus"

#########################################################################

docker run -v ${path}:/ifb/data/mydatalocal/HelmetedCurassowGenome --rm -t quay.io/comparative-genomics-toolkit/cactus:v1.3.0 hal2maf ${pathHAL}/366-avian.hal ${pathResults}/combined_peaks_galGal6_formatted.maf  --targetGenomes ${targetGenomes} --refGenome ${refGenome}  --noDupes --refTargets ${pathResults}/combined_peaks_galGal6_formatted.bed

#########################################################################
