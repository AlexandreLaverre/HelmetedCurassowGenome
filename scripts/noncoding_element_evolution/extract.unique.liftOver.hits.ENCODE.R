######################################################################

peaks=read.table("../../results/noncoding_element_evolution/ENCODE_ATAC-seq/Mouse/combined_peaks_galGal6_merged100bp.bed",h=F, stringsAsFactors=F)

print(paste("all start < end:", all(peaks$V2<peaks$V3)))

######################################################################

dupli=peaks$V4[which(duplicated(peaks$V4))]
peaks=peaks[which(!peaks$V4%in%dupli),]

######################################################################

write.table(peaks[,1:4], file="../../results/noncoding_element_evolution/ENCODE_ATAC-seq/Mouse/combined_peaks_galGal6_unique.bed", row.names=F, col.names=F, quote=F, sep="\t")

######################################################################
