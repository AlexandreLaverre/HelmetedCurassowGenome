#!/bin/bash

## original script by Alexandre Laverré & Anamaria Necsulea

export dataset=$1
export refGenome=$2
export chr=$3
export start=$4
export length=$5
export cluster=$6
export nthreads=$7

#########################################################################

if [ ${cluster} = "cloud" ]; then
    export path=/ifb/data/mydatalocal/HelmetedCurassowGenome
fi

export pathHAL=${path}/results/whole_genome_alignments/avian_366
export pathResults=${path}/results/whole_genome_alignments/${dataset}

#########################################################################

if [ ${dataset} = "no_protuberance" ]; then
    export targetGenomes="Struthio_camelus,Dromaius_novaehollandiae,Gallus_gallus,Meleagris_gallopavo,Penelope_pileata,Alectura_lathami,Anas_platyrhynchos_platyrhynchos,Grus_americana,Calidris_pugnax,Upupa_epops,Rhinopomastus_cyanomelas,Strix_occidentalis,Ficedula_albicollis,Aquila_chrysaetos,Geospiza_fortis,Parus_major,Serinus_canaria,Coturnix_japonica"
fi

if [ ${dataset} = "protuberance" ]; then
    export targetGenomes="Casuarius_casuarius,Numida_meleagris,Anseranas_semipalmata,Anser_cygnoid,Balearica_regulorum,Bucorvus_abyssinicus,Buceros_rhinoceros,Pauxi_pauxi"
fi

#########################################################################

if [ -e ${pathResults}/${refGenome} ]; then
    echo "dir exists"
else
    mkdir -p ${pathResults}/${refGenome}
fi

#########################################################################

if [ ${chr} = "all" ]; then
    docker run -v ${path}:/ifb/data/mydatalocal/HelmetedCurassowGenome --rm -t quay.io/comparative-genomics-toolkit/cactus:v1.3.0 hal2mafMP.py ${pathHAL}/366-avian.hal ${pathResults}/${refGenome}/aln.maf  --targetGenomes ${targetGenomes} --refGenome ${refGenome} --numProc ${nthreads} --noDupes --splitBySequence --smallSize 100000
else
    docker run -v ${path}:/ifb/data/mydatalocal/HelmetedCurassowGenome --rm -t quay.io/comparative-genomics-toolkit/cactus:v1.3.0 hal2mafMP.py ${pathHAL}/366-avian.hal ${pathResults}/${refGenome}/aln.maf  --targetGenomes ${targetGenomes} --refGenome ${refGenome} --refSequence ${chr} --start=${start} --length=${length}  --noDupes 
fi


#########################################################################
