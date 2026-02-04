# ArchR-Transcription-Factor-Motif-Analysis
This repository provides code and resources for performing Transcription Factor (TF) motif analysis on single-cell data using ArchR. The analysis focuses on identifying and visualizing transcription factor binding motifs to uncover key regulatory elements driving gene expression in different cell types or conditions.

## 1. Key Features:
****************************************************************************************************
Key Features:

Motif Discovery: Identifies enriched transcription factor motifs within accessible chromatin regions, linking TFs to regulatory mechanisms.

Motif Activity in Single Cells: Analyzes the activity of TF motifs across single cells, providing insights into cell-type specific transcriptional regulation.

Integration with scRNA-seq and scATAC-seq: Leverages multi-omic data to correlate gene expression with chromatin accessibility at the single-cell level.

Visualization: Generates heatmaps, motif enrichment plots, and motif activity scores to visualize transcription factor dynamics across cells.

Differential Motif Analysis: Identifies differences in TF motif accessibility between conditions, cell types, or developmental stages.

****************************************************************************************************

                                   SCRIPT USAGE  

                                  
library(dplyr)

library(Seurat)

library(scales)

library(cowplot)

library(ggplot2)

library(RColorBrewer)

library(gplots)

library(sctransform)

library(SingleCellExperiment)

library(cluster)

#library(factoextra)

#library(intrinsicDimension)

#library(DoubletFinder)

library(ArchR)

library("monocle3")

#library(BSgenome.Hsapiens.UCSC.hg38)

addArchRThreads(threads = 50)

#addArchRGenome("hg38")

options(future.seed = TRUE)

options(future.globals.maxSize = 800 * 1024^3)


# 1️⃣ Load your updated ArchR project
projEM <- loadArchRProject("/rsrch3/home/genomic_med/ssatpati/Tcell_Exaustion/3.Analysis-Mito0/Common_Cells/ArchR-Save-after-Motif")


# 4️⃣ Differential peaks: LMT vs PT
markers_C4vsC7 <- getMarkerFeatures(
  ArchRProj = projEM,
  useMatrix = "PeakMatrix",
  groupBy = "Clusters",
  maxCells = 200000,
  testMethod = "wilcoxon",
  useGroups = "C7",    # test group
  bgdGroups = "C4"   # background
)

# 6️⃣ Motif enrichment for LMT vs PT
enrich_C4vsC7_Up <- peakAnnoEnrichment(
  seMarker = markers_C4vsC7,
  ArchRProj = projEM,
  peakAnnotation = "Motif",
  cutOff = "FDR <= 0.01 & Log2FC >= 0.5"
)
  
# 7️⃣ Motif enrichment for NLMT vs PT
enrich_C4vsC7_down <- peakAnnoEnrichment(
  seMarker = markers_C4vsC7,
  ArchRProj = projEM,
  peakAnnotation = "Motif",
  cutOff = "FDR <= 0.01 & Log2FC <= -0.5"
)

du <- data.frame(TF = rownames(enrich_C4vsC7_Up), mlog10Padj = assay(enrich_C4vsC7_Up)[,1])
du <- du[order(du$mlog10Padj, decreasing = TRUE),]
du$rank <- seq_len(nrow(du))
head(du)

write.csv(du,file = "TFmotifs_C4_vs_C7_up.csv", row.names = FALSE, quote = FALSE)

ggUp <- ggplot(du, aes(rank, mlog10Padj, color = mlog10Padj)) +
  geom_point(size = 1) +
  ggrepel::geom_label_repel(
        data = du[rev(seq_len(30)), ], aes(x = rank, y = mlog10Padj, label = TF),
        size = 1.5,
        nudge_x = 2,
        color = "black"
  ) + theme_ArchR() +
  ylab("-log10(P-adj) Motif Enrichment") +
  xlab("Rank Sorted TFs Up Enriched") +
  scale_color_gradientn(colors = paletteContinuous(set = "comet"))

ggUp
## Warning: ggrepel: 23 unlabeled data points (too many overlaps). Consider
## increasing max.overlaps

dn <- data.frame(TF = rownames(enrich_C4vsC7_down), mlog10Padj = assay(enrich_C4vsC7_down)[,1])
dn <- dn[order(dn$mlog10Padj, decreasing = TRUE),]
dn$rank <- seq_len(nrow(dn))


head(dn)

write.csv(dn,file = "TFmotifs_C4_vs_C7_down.csv", row.names = FALSE, quote = FALSE)

ggDo <- ggplot(dn, aes(rank, mlog10Padj, color = mlog10Padj)) +
  geom_point(size = 1) +
  ggrepel::geom_label_repel(
        data = dn[rev(seq_len(30)), ], aes(x = rank, y = mlog10Padj, label = TF),
        size = 1.5,
        nudge_x = 2,
        color = "black"
  ) + theme_ArchR() +
  ylab("-log10(FDR) Motif Enrichment") +
  xlab("Rank Sorted TFs Dn Enriched") +
  scale_color_gradientn(colors = paletteContinuous(set = "comet"))

ggDo

plotPDF(ggUp, ggDo, name = "TFmotifs_C4_vs_C7-Markers-Motifs-Enriched", width = 5, height = 5, ArchRProj = projEM, addDOC = FALSE)



####Second##
# 4️⃣ Differential peaks: LMT vs PT
markers_C4vsC8 <- getMarkerFeatures(
  ArchRProj = projEM,
  useMatrix = "PeakMatrix",
  groupBy = "Clusters",
  maxCells = 200000,
  testMethod = "wilcoxon",
  useGroups = "C8",    # test group
  bgdGroups = "C4"   # background
)

# 6️⃣ Motif enrichment for LMT vs PT
enrich_C4vsC8_Up <- peakAnnoEnrichment(
  seMarker = markers_C4vsC8,
  ArchRProj = projEM,
  peakAnnotation = "Motif",
  cutOff = "FDR <= 0.01 & Log2FC >= 0.5"
)

# 7️⃣ Motif enrichment for NLMT vs PT
enrich_C4vsC8_down <- peakAnnoEnrichment(
  seMarker = markers_C4vsC8,
  ArchRProj = projEM,
  peakAnnotation = "Motif",
  cutOff = "FDR <= 0.01 & Log2FC <= -0.5"
)

du <- data.frame(TF = rownames(enrich_C4vsC8_Up), mlog10Padj = assay(enrich_C4vsC8_Up)[,1])
du <- du[order(du$mlog10Padj, decreasing = TRUE),]
du$rank <- seq_len(nrow(du))
head(du)

write.csv(du,file = "TFmotifs_C4_vs_C8_up.csv", row.names = FALSE, quote = FALSE)

ggUp <- ggplot(du, aes(rank, mlog10Padj, color = mlog10Padj)) +
  geom_point(size = 1) +
  ggrepel::geom_label_repel(
        data = du[rev(seq_len(30)), ], aes(x = rank, y = mlog10Padj, label = TF),
        size = 1.5,
        nudge_x = 2,
        color = "black"
  ) + theme_ArchR() +
  ylab("-log10(P-adj) Motif Enrichment") +
  xlab("Rank Sorted TFs Up Enriched") +
  scale_color_gradientn(colors = paletteContinuous(set = "comet"))

ggUp
## Warning: ggrepel: 23 unlabeled data points (too many overlaps). Consider
## increasing max.overlaps
dn <- data.frame(TF = rownames(enrich_C4vsC8_down), mlog10Padj = assay(enrich_C4vsC8_down)[,1])
dn <- dn[order(dn$mlog10Padj, decreasing = TRUE),]
dn$rank <- seq_len(nrow(dn))


head(dn)

write.csv(dn,file = "TFmotifs_C4_vs_C8_down.csv", row.names = FALSE, quote = FALSE)

ggDo <- ggplot(dn, aes(rank, mlog10Padj, color = mlog10Padj)) +
  geom_point(size = 1) +
  ggrepel::geom_label_repel(
        data = dn[rev(seq_len(30)), ], aes(x = rank, y = mlog10Padj, label = TF),
        size = 1.5,
        nudge_x = 2,
        color = "black"
  ) + theme_ArchR() +
  ylab("-log10(FDR) Motif Enrichment") +
  xlab("Rank Sorted TFs Dn Enriched") +
  scale_color_gradientn(colors = paletteContinuous(set = "comet"))

ggDo

plotPDF(ggUp, ggDo, name = "TFmotifs_C4_vs_C8-Markers-Motifs-Enriched", width = 5, height = 5, ArchRProj = projEM, addDOC = FALSE)



#Third

# 4️⃣ Differential peaks: LMT vs PT
markers_C3vsC8 <- getMarkerFeatures(
  ArchRProj = projEM,
  useMatrix = "PeakMatrix",
  groupBy = "Clusters",
  maxCells = 200000,
  testMethod = "wilcoxon",
  useGroups = "C8",    # test group
  bgdGroups = "C3"   # background
)

# 6️⃣ Motif enrichment for LMT vs PT
enrich_C3vsC8_Up <- peakAnnoEnrichment(
  seMarker = markers_C3vsC8,
  ArchRProj = projEM,
  peakAnnotation = "Motif",
  cutOff = "FDR <= 0.01 & Log2FC >= 0.5"
)

# 7️⃣ Motif enrichment for NLMT vs PT
enrich_C3vsC8_down <- peakAnnoEnrichment(
  seMarker = markers_C3vsC8,
  ArchRProj = projEM,
  peakAnnotation = "Motif",
  cutOff = "FDR <= 0.01 & Log2FC <= -0.5"
) 

du <- data.frame(TF = rownames(enrich_C3vsC8_Up), mlog10Padj = assay(enrich_C3vsC8_Up)[,1])
du <- du[order(du$mlog10Padj, decreasing = TRUE),]
du$rank <- seq_len(nrow(du))
head(du)
  
write.csv(du,file = "TFmotifs_C3_vs_C8_up.csv", row.names = FALSE, quote = FALSE)
  
ggUp <- ggplot(du, aes(rank, mlog10Padj, color = mlog10Padj)) +
  geom_point(size = 1) +
  ggrepel::geom_label_repel(
        data = du[rev(seq_len(30)), ], aes(x = rank, y = mlog10Padj, label = TF),
        size = 1.5,
        nudge_x = 2,
        color = "black"
  ) + theme_ArchR() +
  ylab("-log10(P-adj) Motif Enrichment") +
  xlab("Rank Sorted TFs Up Enriched") +
  scale_color_gradientn(colors = paletteContinuous(set = "comet"))

ggUp
## Warning: ggrepel: 23 unlabeled data points (too many overlaps). Consider
## increasing max.overlaps
dn <- data.frame(TF = rownames(enrich_C3vsC8_down), mlog10Padj = assay(enrich_C3vsC8_down)[,1])
dn <- dn[order(dn$mlog10Padj, decreasing = TRUE),]
dn$rank <- seq_len(nrow(dn)) 
  
  
head(dn)
        
write.csv(dn,file = "TFmotifs_C3_vs_C8_down.csv", row.names = FALSE, quote = FALSE)

ggDo <- ggplot(dn, aes(rank, mlog10Padj, color = mlog10Padj)) +
  geom_point(size = 1) +
  ggrepel::geom_label_repel(
        data = dn[rev(seq_len(30)), ], aes(x = rank, y = mlog10Padj, label = TF),
        size = 1.5,
        nudge_x = 2,
        color = "black"
  ) + theme_ArchR() +
  ylab("-log10(FDR) Motif Enrichment") +
  xlab("Rank Sorted TFs Dn Enriched") +
  scale_color_gradientn(colors = paletteContinuous(set = "comet"))

ggDo

plotPDF(ggUp, ggDo, name = "TFmotifs_C3_vs_C8-Markers-Motifs-Enriched", width = 5, height = 5, ArchRProj = projEM, addDOC = FALSE)


#Fourth
# 4️⃣ Differential peaks: LMT vs PT
markers_C3vsC7 <- getMarkerFeatures(
  ArchRProj = projEM,
  useMatrix = "PeakMatrix",
  groupBy = "Clusters",
  maxCells = 200000,
  testMethod = "wilcoxon",
  useGroups = "C7",    # test group
  bgdGroups = "C3"   # background
)

# 6️⃣ Motif enrichment for LMT vs PT
enrich_C3vsC7_Up <- peakAnnoEnrichment(
  seMarker = markers_C3vsC7,
  ArchRProj = projEM,
  peakAnnotation = "Motif",
  cutOff = "FDR <= 0.01 & Log2FC >= 0.5"
)

# 7️⃣ Motif enrichment for NLMT vs PT
enrich_C3vsC7_down <- peakAnnoEnrichment(
  seMarker = markers_C3vsC7,
  ArchRProj = projEM,
  peakAnnotation = "Motif",
  cutOff = "FDR <= 0.01 & Log2FC <= -0.5"
)

du <- data.frame(TF = rownames(enrich_C3vsC7_Up), mlog10Padj = assay(enrich_C3vsC7_Up)[,1])
du <- du[order(du$mlog10Padj, decreasing = TRUE),]
du$rank <- seq_len(nrow(du))
head(du)

write.csv(du,file = "TFmotifs_C3_vs_C7_up.csv", row.names = FALSE, quote = FALSE)

ggUp <- ggplot(du, aes(rank, mlog10Padj, color = mlog10Padj)) +
  geom_point(size = 1) +
  ggrepel::geom_label_repel(
        data = du[rev(seq_len(30)), ], aes(x = rank, y = mlog10Padj, label = TF),
        size = 1.5,
        nudge_x = 2,
        color = "black"
  ) + theme_ArchR() +
  ylab("-log10(P-adj) Motif Enrichment") +
  xlab("Rank Sorted TFs Up Enriched") +
  scale_color_gradientn(colors = paletteContinuous(set = "comet"))

ggUp
## Warning: ggrepel: 23 unlabeled data points (too many overlaps). Consider
## increasing max.overlaps
dn <- data.frame(TF = rownames(enrich_C3vsC7_down), mlog10Padj = assay(enrich_C3vsC7_down)[,1])
dn <- dn[order(dn$mlog10Padj, decreasing = TRUE),]
dn$rank <- seq_len(nrow(dn))


head(dn)

write.csv(dn,file = "TFmotifs_C3_vs_C7_down.csv", row.names = FALSE, quote = FALSE)

ggDo <- ggplot(dn, aes(rank, mlog10Padj, color = mlog10Padj)) +
  geom_point(size = 1) +
  ggrepel::geom_label_repel(
        data = dn[rev(seq_len(30)), ], aes(x = rank, y = mlog10Padj, label = TF),
        size = 1.5,
        nudge_x = 2,
        color = "black"
  ) + theme_ArchR() +
  ylab("-log10(FDR) Motif Enrichment") +
  xlab("Rank Sorted TFs Dn Enriched") +
  scale_color_gradientn(colors = paletteContinuous(set = "comet"))

ggDo

plotPDF(ggUp, ggDo, name = "TFmotifs_C3_vs_C7-Markers-Motifs-Enriched", width = 5, height = 5, ArchRProj = projEM, addDOC = FALSE)


#New_one
# 4️⃣ Differential peaks: LMT vs PT
markers_C3vsC9 <- getMarkerFeatures(
  ArchRProj = projEM,
  useMatrix = "PeakMatrix",
  groupBy = "Clusters",
  maxCells = 200000,
  testMethod = "wilcoxon",
  useGroups = "C9",    # test group
  bgdGroups = "C3"   # background
)

# 6️⃣ Motif enrichment for LMT vs PT
enrich_C3vsC9_Up <- peakAnnoEnrichment(
  seMarker = markers_C3vsC9,
  ArchRProj = projEM,
  peakAnnotation = "Motif",
  cutOff = "FDR <= 0.01 & Log2FC >= 0.5"
)

# 7️⃣ Motif enrichment for NLMT vs PT
enrich_C3vsC9_down <- peakAnnoEnrichment(
  seMarker = markers_C3vsC9,
  ArchRProj = projEM,
  peakAnnotation = "Motif",
  cutOff = "FDR <= 0.01 & Log2FC <= -0.5"
)

du <- data.frame(TF = rownames(enrich_C3vsC9_Up), mlog10Padj = assay(enrich_C3vsC9_Up)[,1])
du <- du[order(du$mlog10Padj, decreasing = TRUE),]
du$rank <- seq_len(nrow(du))
head(du)

write.csv(du,file = "TFmotifs_C3_vs_C9_up.csv", row.names = FALSE, quote = FALSE)

ggUp <- ggplot(du, aes(rank, mlog10Padj, color = mlog10Padj)) +
  geom_point(size = 1) +
  ggrepel::geom_label_repel(
        data = du[rev(seq_len(30)), ], aes(x = rank, y = mlog10Padj, label = TF),
        size = 1.5,
        nudge_x = 2,
        color = "black"
  ) + theme_ArchR() +
  ylab("-log10(P-adj) Motif Enrichment") +
  xlab("Rank Sorted TFs Up Enriched") +
  scale_color_gradientn(colors = paletteContinuous(set = "comet"))

ggUp
## Warning: ggrepel: 23 unlabeled data points (too many overlaps). Consider
## increasing max.overlaps
dn <- data.frame(TF = rownames(enrich_C3vsC9_down), mlog10Padj = assay(enrich_C3vsC9_down)[,1])
dn <- dn[order(dn$mlog10Padj, decreasing = TRUE),]
dn$rank <- seq_len(nrow(dn))


head(dn)

write.csv(dn,file = "TFmotifs_C3_vs_C9_down.csv", row.names = FALSE, quote = FALSE)


ggDo <- ggplot(dn, aes(rank, mlog10Padj, color = mlog10Padj)) +
  geom_point(size = 1) +
  ggrepel::geom_label_repel(
        data = dn[rev(seq_len(30)), ], aes(x = rank, y = mlog10Padj, label = TF),
        size = 1.5,
        nudge_x = 2,
        color = "black"
  ) + theme_ArchR() +
  ylab("-log10(FDR) Motif Enrichment") +
  xlab("Rank Sorted TFs Dn Enriched") +
  scale_color_gradientn(colors = paletteContinuous(set = "comet"))

ggDo

plotPDF(ggUp, ggDo, name = "TFmotifs_C3_vs_C9-Markers-Motifs-Enriched", width = 5, height = 5, ArchRProj = projEM, addDOC = FALSE)



#Newtwo
# 4️⃣ Differential peaks: LMT vs PT
markers_C3vsC10 <- getMarkerFeatures(
  ArchRProj = projEM,
  useMatrix = "PeakMatrix",
  groupBy = "Clusters",
  maxCells = 200000,
  testMethod = "wilcoxon",
  useGroups = "C10",    # test group
  bgdGroups = "C3"   # background
)

# 6️⃣ Motif enrichment for LMT vs PT
enrich_C3vsC10_Up <- peakAnnoEnrichment(
  seMarker = markers_C3vsC10,
  ArchRProj = projEM,
  peakAnnotation = "Motif",
  cutOff = "FDR <= 0.01 & Log2FC >= 0.5"
)

# 7️⃣ Motif enrichment for NLMT vs PT
enrich_C3vsC10_down <- peakAnnoEnrichment(
  seMarker = markers_C3vsC10,
  ArchRProj = projEM,
  peakAnnotation = "Motif",
  cutOff = "FDR <= 0.01 & Log2FC <= -0.5"
)

du <- data.frame(TF = rownames(enrich_C3vsC10_Up), mlog10Padj = assay(enrich_C3vsC10_Up)[,1])
du <- du[order(du$mlog10Padj, decreasing = TRUE),]
du$rank <- seq_len(nrow(du))
head(du)

write.csv(du,file = "TFmotifs_C3_vs_C10_up.csv", row.names = FALSE, quote = FALSE)

ggUp <- ggplot(du, aes(rank, mlog10Padj, color = mlog10Padj)) +
  geom_point(size = 1) +
  ggrepel::geom_label_repel(
        data = du[rev(seq_len(30)), ], aes(x = rank, y = mlog10Padj, label = TF),
        size = 1.5,
        nudge_x = 2,
        color = "black"
  ) + theme_ArchR() +
  ylab("-log10(P-adj) Motif Enrichment") +
  xlab("Rank Sorted TFs Up Enriched") +
  scale_color_gradientn(colors = paletteContinuous(set = "comet"))

ggUp
## Warning: ggrepel: 23 unlabeled data points (too many overlaps). Consider
## increasing max.overlaps
dn <- data.frame(TF = rownames(enrich_C3vsC10_down), mlog10Padj = assay(enrich_C3vsC10_down)[,1])
dn <- dn[order(dn$mlog10Padj, decreasing = TRUE),]
dn$rank <- seq_len(nrow(dn))


head(dn)

write.csv(dn,file = "TFmotifs_C3_vs_C10_down.csv", row.names = FALSE, quote = FALSE)

ggDo <- ggplot(dn, aes(rank, mlog10Padj, color = mlog10Padj)) +
  geom_point(size = 1) +
  ggrepel::geom_label_repel(
        data = dn[rev(seq_len(30)), ], aes(x = rank, y = mlog10Padj, label = TF),
        size = 1.5,
        nudge_x = 2,
        color = "black"
  ) + theme_ArchR() +
  ylab("-log10(FDR) Motif Enrichment") +
  xlab("Rank Sorted TFs Dn Enriched") +
  scale_color_gradientn(colors = paletteContinuous(set = "comet"))

ggDo

plotPDF(ggUp, ggDo, name = "TFmotifs_C3_vs_C10-Markers-Motifs-Enriched", width = 5, height = 5, ArchRProj = projEM, addDOC = FALSE)

#other_one
# 4️⃣ Differential peaks: LMT vs PT
markers_C4vsC10 <- getMarkerFeatures(
  ArchRProj = projEM,
  useMatrix = "PeakMatrix",
  groupBy = "Clusters",
  maxCells = 200000,
  testMethod = "wilcoxon",
  useGroups = "C10",    # test group
  bgdGroups = "C4"   # background
)

# 6️⃣ Motif enrichment for LMT vs PT
enrich_C4vsC10_Up <- peakAnnoEnrichment(
  seMarker = markers_C4vsC10,
  ArchRProj = projEM,
  peakAnnotation = "Motif",
  cutOff = "FDR <= 0.01 & Log2FC >= 0.5"
)

# 7️⃣ Motif enrichment for NLMT vs PT
enrich_C4vsC10_down <- peakAnnoEnrichment(
  seMarker = markers_C4vsC10,
  ArchRProj = projEM,
  peakAnnotation = "Motif",
  cutOff = "FDR <= 0.01 & Log2FC <= -0.5"
)

du <- data.frame(TF = rownames(enrich_C4vsC10_Up), mlog10Padj = assay(enrich_C4vsC10_Up)[,1])
du <- du[order(du$mlog10Padj, decreasing = TRUE),]
du$rank <- seq_len(nrow(du))
head(du)

write.csv(du,file = "TFmotifs_C4_vs_C10_up.csv", row.names = FALSE, quote = FALSE)

ggUp <- ggplot(du, aes(rank, mlog10Padj, color = mlog10Padj)) +
  geom_point(size = 1) +
  ggrepel::geom_label_repel(
        data = du[rev(seq_len(30)), ], aes(x = rank, y = mlog10Padj, label = TF),
        size = 1.5,
        nudge_x = 2,
        color = "black"
  ) + theme_ArchR() +
  ylab("-log10(P-adj) Motif Enrichment") +
  xlab("Rank Sorted TFs Up Enriched") +
  scale_color_gradientn(colors = paletteContinuous(set = "comet"))

ggUp
## Warning: ggrepel: 23 unlabeled data points (too many overlaps). Consider
## increasing max.overlaps
dn <- data.frame(TF = rownames(enrich_C4vsC10_down), mlog10Padj = assay(enrich_C4vsC10_down)[,1])
dn <- dn[order(dn$mlog10Padj, decreasing = TRUE),]
dn$rank <- seq_len(nrow(dn))


head(dn)

write.csv(dn,file = "TFmotifs_C4_vs_C10_down.csv", row.names = FALSE, quote = FALSE)

ggDo <- ggplot(dn, aes(rank, mlog10Padj, color = mlog10Padj)) +
  geom_point(size = 1) +
  ggrepel::geom_label_repel(
        data = dn[rev(seq_len(30)), ], aes(x = rank, y = mlog10Padj, label = TF),
        size = 1.5,
        nudge_x = 2,
        color = "black"
  ) + theme_ArchR() +
  ylab("-log10(FDR) Motif Enrichment") +
  xlab("Rank Sorted TFs Dn Enriched") +
  scale_color_gradientn(colors = paletteContinuous(set = "comet"))

ggDo

plotPDF(ggUp, ggDo, name = "TFmotifs_C4_vs_C10-Markers-Motifs-Enriched", width = 5, height = 5, ArchRProj = projEM, addDOC = FALSE)


#othertwo
# 4️⃣ Differential peaks: LMT vs PT
markers_C4vsC9 <- getMarkerFeatures(
  ArchRProj = projEM,
  useMatrix = "PeakMatrix",
  groupBy = "Clusters",
  maxCells = 200000,
  testMethod = "wilcoxon",
  useGroups = "C9",    # test group
  bgdGroups = "C4"   # background
)

# 6️⃣ Motif enrichment for LMT vs PT
enrich_C4vsC9_Up <- peakAnnoEnrichment(
  seMarker = markers_C4vsC9,
  ArchRProj = projEM,
  peakAnnotation = "Motif",
  cutOff = "FDR <= 0.01 & Log2FC >= 0.5"
)

# 7️⃣ Motif enrichment for NLMT vs PT
enrich_C4vsC9_down <- peakAnnoEnrichment(
  seMarker = markers_C4vsC9,
  ArchRProj = projEM,
  peakAnnotation = "Motif",
  cutOff = "FDR <= 0.01 & Log2FC <= -0.5"
)

du <- data.frame(TF = rownames(enrich_C4vsC9_Up), mlog10Padj = assay(enrich_C4vsC9_Up)[,1])
du <- du[order(du$mlog10Padj, decreasing = TRUE),]
du$rank <- seq_len(nrow(du))
head(du)

write.csv(du,file = "TFmotifs_C4_vs_C9_up.csv", row.names = FALSE, quote = FALSE)

ggUp <- ggplot(du, aes(rank, mlog10Padj, color = mlog10Padj)) +
  geom_point(size = 1) +
  ggrepel::geom_label_repel(
        data = du[rev(seq_len(30)), ], aes(x = rank, y = mlog10Padj, label = TF),
        size = 1.5,
        nudge_x = 2,
        color = "black"
  ) + theme_ArchR() +
  ylab("-log10(P-adj) Motif Enrichment") +
  xlab("Rank Sorted TFs Up Enriched") +
  scale_color_gradientn(colors = paletteContinuous(set = "comet"))

ggUp
## Warning: ggrepel: 23 unlabeled data points (too many overlaps). Consider
## increasing max.overlaps
dn <- data.frame(TF = rownames(enrich_C4vsC9_down), mlog10Padj = assay(enrich_C4vsC9_down)[,1])
dn <- dn[order(dn$mlog10Padj, decreasing = TRUE),]
dn$rank <- seq_len(nrow(dn))


head(dn)

write.csv(dn,file = "TFmotifs_C4_vs_C9_down.csv", row.names = FALSE, quote = FALSE)

ggDo <- ggplot(dn, aes(rank, mlog10Padj, color = mlog10Padj)) +
  geom_point(size = 1) +
  ggrepel::geom_label_repel(
        data = dn[rev(seq_len(30)), ], aes(x = rank, y = mlog10Padj, label = TF),
        size = 1.5,
        nudge_x = 2,
        color = "black"
  ) + theme_ArchR() +
  ylab("-log10(FDR) Motif Enrichment") +
  xlab("Rank Sorted TFs Dn Enriched") +
  scale_color_gradientn(colors = paletteContinuous(set = "comet"))

ggDo

plotPDF(ggUp, ggDo, name = "TFmotifs_C4_vs_C9-Markers-Motifs-Enriched", width = 5, height = 5, ArchRProj = projEM, addDOC = FALSE)

#C3 vs C4

# 4️⃣ Differential peaks: LMT vs PT
markers_C4vsC3 <- getMarkerFeatures(
  ArchRProj = projEM,
  useMatrix = "PeakMatrix",
  groupBy = "Clusters",
  maxCells = 200000,
  testMethod = "wilcoxon",
  useGroups = "C3",    # test group
  bgdGroups = "C4"   # background
)

# 6️⃣ Motif enrichment for LMT vs PT
enrich_C4vsC3_Up <- peakAnnoEnrichment(
  seMarker = markers_C4vsC3,
  ArchRProj = projEM,
  peakAnnotation = "Motif",
  cutOff = "FDR <= 0.01 & Log2FC >= 0.5"
)

# 7️⃣ Motif enrichment for NLMT vs PT
enrich_C4vsC3_down <- peakAnnoEnrichment(
  seMarker = markers_C4vsC3,
  ArchRProj = projEM,
  peakAnnotation = "Motif",
  cutOff = "FDR <= 0.01 & Log2FC <= -0.5"
)

du <- data.frame(TF = rownames(enrich_C4vsC3_Up), mlog10Padj = assay(enrich_C4vsC3_Up)[,1])
du <- du[order(du$mlog10Padj, decreasing = TRUE),]
du$rank <- seq_len(nrow(du))
head(du)

write.csv(du,file = "TFmotifs_C4_vs_C3_up.csv", row.names = FALSE, quote = FALSE)

ggUp <- ggplot(du, aes(rank, mlog10Padj, color = mlog10Padj)) +
  geom_point(size = 1) +
  ggrepel::geom_label_repel(
        data = du[rev(seq_len(30)), ], aes(x = rank, y = mlog10Padj, label = TF),
        size = 1.5,
        nudge_x = 2,
        color = "black"
  ) + theme_ArchR() +
  ylab("-log10(P-adj) Motif Enrichment") +
  xlab("Rank Sorted TFs Up Enriched") +
  scale_color_gradientn(colors = paletteContinuous(set = "comet"))

ggUp
## Warning: ggrepel: 23 unlabeled data points (too many overlaps). Consider
## increasing max.overlaps
dn <- data.frame(TF = rownames(enrich_C4vsC3_down), mlog10Padj = assay(enrich_C4vsC3_down)[,1])
dn <- dn[order(dn$mlog10Padj, decreasing = TRUE),]
dn$rank <- seq_len(nrow(dn))


head(dn)

write.csv(dn,file = "TFmotifs_C4_vs_C3_down.csv", row.names = FALSE, quote = FALSE)

ggDo <- ggplot(dn, aes(rank, mlog10Padj, color = mlog10Padj)) +
  geom_point(size = 1) +
  ggrepel::geom_label_repel(
        data = dn[rev(seq_len(30)), ], aes(x = rank, y = mlog10Padj, label = TF),
        size = 1.5,
        nudge_x = 2,
        color = "black"
  ) + theme_ArchR() +
  ylab("-log10(FDR) Motif Enrichment") +
  xlab("Rank Sorted TFs Dn Enriched") +
  scale_color_gradientn(colors = paletteContinuous(set = "comet"))

ggDo

plotPDF(ggUp, ggDo, name = "TFmotifs_C4_vs_C3-Markers-Motifs-Enriched", width = 5, height = 5, ArchRProj = projEM, addDOC = FALSE)


saveArchRProject(
  ArchRProj = projEM,
  outputDirectory = "Save-ProjEM-TF-analyses",
  load = FALSE
)


#Monocle3_Trajectory

activated_to_exhaustion_rna <- getMonocleTrajectories(
    ArchRProj = projEM,
    name = "Activated_to_Exhaustion_Pathway_RNA",
    useGroups = c("C4","C3","C5","C6","C10","C9","C7", "C8"),
    principalGroup = "C4",
    groupBy = "Clusters",
    embedding = "UMAP_RNA",
    clusterParams = list(k = 50),
    seed = 1
)

projEM <- addMonocleTrajectory(
    ArchRProj = projEM,
    name = "activated_to_exhaustion_rna",
    useGroups = c("C4","C3","C5","C6","C10","C9","C7", "C8"),
    groupBy = "Clusters",
    monocleCDS = activated_to_exhaustion_rna,
    force = TRUE
)

p_exhaustion <- plotTrajectory(
    ArchRProj = projEM,
    trajectory = "activated_to_exhaustion_rna",
    colorBy = "cellColData",
    name = "activated_to_exhaustion_rna",
    embedding = "UMAP_RNA",
    addArrow = TRUE
)

png("Activated_to_Exhaustion_Trajectory_RNA.png", width = 2000, height = 2000, res = 300)
print(p_exhaustion[[1]])
dev.off()


activated_to_exhaustion_atac <- getMonocleTrajectories(
    ArchRProj = projEM,
    name = "Activated_to_Exhaustion_Pathway_ATAC",
    useGroups = c("C4","C3","C5","C6","C10","C9","C7", "C8"),
    principalGroup = "C4",
    groupBy = "Clusters",
    embedding = "UMAP_ATAC",
    clusterParams = list(k = 50),
    seed = 1
)

projEM <- addMonocleTrajectory(
    ArchRProj = projEM,
    name = "activated_to_exhaustion_atac",
    useGroups = c("C4","C3","C5","C6","C10","C9","C7", "C8"),
    groupBy = "Clusters",
    monocleCDS = activated_to_exhaustion_atac,
    force = TRUE
)

p_exhaustion1 <- plotTrajectory(
    ArchRProj = projEM,
    trajectory = "activated_to_exhaustion_atac",
    colorBy = "cellColData",
    name = "activated_to_exhaustion_atac",
    embedding = "UMAP_ATAC",
    addArrow = TRUE
)

png("Activated_to_Exhaustion_Trajectory_ATAC.png", width = 2000, height = 2000, res = 300)
print(p_exhaustion1[[1]])
dev.off()

saveArchRProject(
  ArchRProj = projEM,
  outputDirectory = "Save-ProjEM-monocle-analyses",
  load = FALSE
)



