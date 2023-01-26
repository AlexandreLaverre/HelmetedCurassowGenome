#!/usr/bin/env python3
# coding=utf-8
import sys
from Bio import AlignIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align import MultipleSeqAlignment

###################################################################################################

PathAlign = sys.argv[1]
IDList = sys.argv[2]
PathResults = sys.argv[3]
minLength = sys.argv[4]

Part = IDList.split("/")[-1].strip(".txt")
pathFiltered = PathResults + "/FilteredAlign/" + Part + "/"

###################################################################################################
error_align = open(pathFiltered + "noalignment.txt", "a")

with open(IDList, 'r') as f:
    for ID in f.readlines():
        ID = ID.strip("\n")
        ids = []
        seqs = []

        align = AlignIO.read(PathAlign + "/" + ID + ".aln.best_hmm.fasta.fsa", "fasta")
        align_length = align.get_alignment_length()

        if align_length > 0:
            # Keep sequences with low number of GAP
            for seq in align:
                if seq.seq.count("-")/align_length < float(minLength):
                    ids.append(seq.id)
                    seqs.append(seq.seq)

            records = (SeqRecord(s, id=i) for (i, s) in zip(ids, seqs))

            if len(ids) > 1:
                with open(pathFiltered + ID + ".rphylip", "w") as OutFile:
                    AlignIO.write(MultipleSeqAlignment(records), OutFile, "phylip-relaxed")
        else:
            error_align.write(ID)

error_align.close()
###################################################################################################
