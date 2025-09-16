# HelmetedCurassowGenome

De novo genome assembly and evolutionary analyses of the Helmeted Curassow (*Pauxi pauxi*).  
This project aimed to generate a high-quality reference genome for this endangered bird species and to investigate gene family evolution, protein evolution, conserved elements, and molecular convergence across birds and squamates.

## Project Overview

1. **Genome Assembly**
   - Raw Illumina reads assembled with MEGAHIT.
   - Scaffolding performed with Ragout.

2. **Genome Annotation**
   - Gene annotation with GeMoMa using homology-based predictions.
   - Repeat annotation with RepeatModeler and RepeatMasker.

3. **Comparative Genomics**
   - Reference taxa: *Gallus gallus* (chicken) and *Anolis carolinensis* (green anole). 
   - Whole-genome alignments with Progressive Cactus.
   - Orthogroup inference with OrthoFinder.   
   - Detection of conserved elements with phastCons and phyloP.  
   - Copy number variation analysis using LUMPY and CNVnator.  
   - Synteny analysis across reference bird genomes.  
   - Analysis of sequence composition (GC content, repeat density, etc.).

5. **Evolutionary Analyses**
   - **Gene evolution:** phylogenetic reconstruction with IQ-TREE and PhyML.   
   - Gene family evolution with CAFE5.  
   - Selection analyses with PAML.  
   - Protein evolutionary rate with RERconverge.  
   - Convergent amino acid profiles with Pelican.
  
   - **Conserved element evolution:** phylogenetic reconstruction with IQ-TREE and PhyML; rate tests with RERconverge. 
