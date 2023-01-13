#!/usr/bin/env python3
# coding=utf-8
import sys
import os
from pathlib import Path
from Bio import AlignIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align import MultipleSeqAlignment

###################################################################################################

PathResult = sys.argv[1]
PathPHYLIP = sys.argv[2]
Prefix = sys.argv[3]
minLength = sys.argv[4]

Part = PathPHYLIP.split("/")[-1]

SpeciesList = open(PathResult + "/Genomes.txt", "r")
Species = SpeciesList.readline().strip("\n").split(",")

pathCountSeq = PathResult + "/Count_by_species/" + Part + "_FilteredPhylip.txt"
pathFilteredPhylip = PathResult + "/PHYLIPs/NoShorts/" + Part + "/"

###################################################################################################

PHYLIP_files = Path(PathPHYLIP).glob('**/*.phylip')
for phylip in PHYLIP_files:
    phylip_ID = os.path.split(phylip)[1].strip('.phylip')
    with open(phylip, "r") as InFile:
        line = InFile.readline().rstrip()
        parts = [x for x in line.split() if x]
        number_of_seqs = int(parts[0])
        length_of_seqs = int(parts[1])
        ids = []
        seqs = []

        for i in range(number_of_seqs):
            id = InFile.readline().rstrip()
            seq = InFile.readline().rstrip()

            while len(seq) < length_of_seqs:
                line = InFile.readline().strip()
                if not line:
                    break
                if line == "":
                    continue
                seq = "".join([seq, line.strip().replace(" ", "")])
                if len(seq) > length_of_seqs:
                    raise ValueError("Found a record of length %i, ""should be %i" % (len(seq), length_of_seqs))

            if seq.count("-")/length_of_seqs < float(minLength):
                ids.append(id)
                seqs.append(seq)

        #print("Number of original seq:", number_of_seqs, "; Number of filtered seqs:", len(ids))
        records = (SeqRecord(Seq(s), id=i) for (i, s) in zip(ids, seqs))

    if len(ids) > 1:
        with open(pathFilteredPhylip + phylip_ID + ".rphylip", "w") as OutFile:
            AlignIO.write(MultipleSeqAlignment(records), OutFile, "phylip-relaxed")

    #### Count species
    dic_sp = {}
    for sp in Species:
        dic_sp[sp] = str(ids.count(sp))

    # Writing output
    with open(pathCountSeq, "a") as countSeq:
        if os.stat(pathCountSeq).st_size == 0:
            countSeq.write("element.ID" + "\t" + '\t'.join(dic_sp.keys()) + "\n")

        countSeq.write(phylip_ID + "\t" + '\t'.join(dic_sp.values()) + "\n")

###################################################################################################
