######################################################################

peaks=read.table("../../results/noncoding_element_evolution/ENCODE_ATAC-seq/Mouse/combined_peaks_galGal6_merged100bp.bed",h=F, stringsAsFactors=F)

print(paste("all start < end:", all(peaks$V2<peaks$V3)))

######################################################################

dupli=peaks$V4[which(duplicated(peaks$V4))]
peaks=peaks[which(!peaks$V4%in%dupli),]

######################################################################

write.table(peaks[,1:4], file="../../results/noncoding_element_evolution/ENCODE_ATAC-seq/Mouse/combined_peaks_galGal6_unique.bed", row.names=F, col.names=F, quote=F, sep="\t")

######################################################################

accepted.chr=paste("chr", c(as.character(1:33), "Z", "W"), sep="")

peaks=peaks[which(peaks$V1%in%accepted.chr),]
peaks$V1=paste("Gallus_gallus.",unlist(lapply(peaks$V1, function(x) substr(x, 4, nchar(x)))), sep="")

write.table(peaks,file="../../results/noncoding_element_evolution/ENCODE_ATAC-seq/Mouse/combined_peaks_galGal6_formatted.bed",row.names=F, col.names=F, quote=F, sep="\t")

######################################################################
