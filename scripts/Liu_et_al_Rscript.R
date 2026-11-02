# Set input directories
setwd("/Users/liumin/Desktop/R_Input/")

# =================================================================
# Alpha diversity analysis
# Indices: Observed richness,Shannon-Wiener index and Pielou index
# =================================================================
install.packages("vegan")
library(vegan)
# Import ASV abundance table
asv<-read.delim("d1_ASV_table.txt",header = TRUE,row.names = 1,sep = "\t",check.names = FALSE,stringsAsFactors = FALSE)
asv<-t(asv)
#Calculate alpha diversity indices
Observed_richness <- vegan::specnumber(asv)
Shannon_Wiener <- vegan::diversity(asv,index = "shannon",base = exp(1))
Pielou <- Shannon_Wiener / log(Observed_richness)
# Combine results
alpha_diversity <- data.frame(SampleID = rownames(asv),`Observed richness` = as.numeric(Observed_richness),`Shannon-Wiener index` = as.numeric(Shannon_Wiener),`Pielou index` = as.numeric(Pielou), row.names = NULL,check.names = FALSE)
# Export results
write.csv(alpha_diversity, "d1_alpha_diversity.csv", row.names = FALSE)

# ======================================================================================================================================
# Paired t-test of alpha diversity and Figure 3 visualization
# A total of 25 pairs of young (YL) and old (OL) leaf samples (YL1–OL1 to YL25–OL25; 50 samples in total) were included in the analysis.
# Figure 3:
# Bar plots showing mean ± standard error (SE).
# Different letters indicate significant differences (P < 0.05).
# ======================================================================================================================================
install.packages(c("ggplot2", "patchwork"))
library(ggplot2)
library(patchwork)
# Import data and match paired samples
alpha_diversity <- read.csv("d1_alpha_diversity.csv",header = TRUE,check.names = FALSE,stringsAsFactors = FALSE)
alpha_diversity$Group <- sub("[0-9]+$", "", alpha_diversity$SampleID)
alpha_diversity$PairID <- as.integer(sub("^[A-Za-z]+", "", alpha_diversity$SampleID))
YL <- alpha_diversity[alpha_diversity$Group == "YL", ]
OL <- alpha_diversity[alpha_diversity$Group == "OL", ]
YL <- YL[order(YL$PairID), ]
OL <- OL[order(OL$PairID), ]
indices <- c("Observed richness","Shannon-Wiener index","Pielou index")
# Shapiro-Wilk normality tests
normality_results <- do.call(rbind,lapply(indices,function(i){
test <- shapiro.test(YL[[i]]-OL[[i]])
data.frame(Index=i,W=unname(test$statistic),P_value=test$p.value)
}))
# Paired t-tests
ttest_results <- do.call(rbind,lapply(indices,function(i){
test <- t.test(YL[[i]],OL[[i]],paired=TRUE)
data.frame(Index=i,t=unname(test$statistic),P_value=test$p.value)
}))
# Benjamini-Hochberg correction
ttest_results$P_adjusted <- p.adjust(ttest_results$P_value,method="BH")
# Calculate mean and standard error
summary_results <- do.call(rbind,lapply(indices,function(i){
data.frame(Index=i,Group=c("YL","OL"),
Mean=c(mean(YL[[i]]),mean(OL[[i]])),
SE=c(sd(YL[[i]]),sd(OL[[i]]))/sqrt(nrow(YL)))
}))
# Prepare significance letters
plot_data <- merge(summary_results,ttest_results[,c("Index","P_value")],by="Index")
plot_data$Letter <- "a"
for(i in indices){
rows <- plot_data$Index==i
p <- ttest_results$P_adjusted[ttest_results$Index==i]
if(p<0.05){
m <- plot_data$Mean[rows]
plot_data$Letter[rows] <- ifelse(m==max(m),"a","b")
  }
}
# Create bar plots
plots <- lapply(indices, function(index) {
df <- plot_data[plot_data$Index == index, ]
df$Group <- factor(df$Group, levels = c("YL", "OL"))
ymax <- max(df$Mean + df$SE)
upper <- ymax * 1.08
ggplot(df, aes(x = Group, y = Mean, fill = Group)) + geom_col(width = 0.8, color = "grey70", linewidth = 0.3) + geom_errorbar(aes(ymin = Mean - SE, ymax = Mean + SE), width = 0.12, linewidth = 0.6, color = "black") + geom_text(aes(y = Mean + SE, label = Letter), vjust = -0.4, size = 6, family = "serif") + scale_fill_manual(values = c("YL" = "blue", "OL" = "red")) + scale_y_continuous(limits = c(0, upper), expand = expansion(mult = c(0, 0))) + labs(x = NULL, y = index) + theme_classic(base_size = 16, base_family = "serif") + theme(legend.position = "none", panel.border = element_rect(color = "black", fill = NA, linewidth = 0.6), axis.line = element_blank(), axis.ticks = element_line(color = "black", linewidth = 0.6), axis.ticks.length = unit(0.20, "cm"), axis.text = element_text(color = "black", size = 14), axis.text.y = element_text(margin = margin(r = 2)), axis.text.x = element_text(margin = margin(t = 3)), axis.title.y = element_text(color = "black", size = 16), plot.margin = margin(5, 5, 5, 5))
})
# Combine three plots
Figure3 <- wrap_plots(plots,ncol = 3)
# Export Figure 3
ggsave("Figure 3.pdf", plot = Figure3,width = 12, height = 4.5,units = "in", device = cairo_pdf)

# ======================================================================================================================================
# Indicator ASV analysis and bipartite network visualization
# A total of 25 pairs of young (YL) and old (OL) leaf samples (YL1–OL1 to YL25–OL25; 50 samples in total) were included in the analysis.
# Permutations were restricted to the YL and OL samples within each individual tree.
# ======================================================================================================================================
install.packages(c("indicspecies", "permute", "igraph"))
library(indicspecies)
library(permute)
library(igraph)
# Data input and preparation
asv <- read.delim("d1_ASV_table.txt", header=TRUE, row.names=1, check.names=FALSE, stringsAsFactors=FALSE)
taxonomy <- read.delim("d2_taxonomy_table.txt", header=TRUE, check.names=FALSE, stringsAsFactors=FALSE)
metadata <- read.delim("d3_metadata_table.txt", header=TRUE, check.names=FALSE, stringsAsFactors=FALSE)
asv <- t(as.matrix(asv))
storage.mode(asv) <- "numeric"
taxonomy$Class <- trimws(sub("^c__","",taxonomy$Class))
taxonomy$Class[is.na(taxonomy$Class) | taxonomy$Class==""] <- "Unclassified"
# Indicator species analysis
set.seed(123)
perm <- how(nperm=999, within=Within(type="free"), blocks=factor(metadata$PairID))
ind <- multipatt(asv, factor(metadata$Group,levels=c("YL","OL")), func="IndVal.g", max.order=2, control=perm)
# Extract significant indicator ASVs
results <- data.frame(ASV=rownames(ind$sign), ind$sign, check.names=FALSE)
sig <- results[!is.na(results$p.value) & results$p.value < 0.05,,drop=FALSE]
if(nrow(sig)==0) stop("No significant indicator ASVs detected.")
sig$Group <- apply(sig[,c("s.YL","s.OL"),drop=FALSE], 1, function(x) paste(c("YL","OL")[x==1],collapse="+"))
# Assign fungal classes to significant indicator ASVs
sig$Class <- taxonomy$Class[match(sig$ASV,taxonomy[["ASV ID"]])]
sig$Class[is.na(sig$Class) | sig$Class==""] <- "Unclassified"
# Bipartite network construction
edges <- do.call(rbind, lapply(c("YL","OL"), function(group){
x <- sig[sig[[paste0("s.",group)]]==1,,drop=FALSE]
data.frame(from=rep(group,nrow(x)), to=x$ASV, r=x$stat, stringsAsFactors=FALSE)
}))
if(nrow(edges)==0) stop("No network edges detected.")
g <- graph_from_data_frame(edges, directed=FALSE)
V(g)$type <- V(g)$name %in% c("YL","OL")
base_cols <- c(Agaricomycetes="#950095", Agaricostilbomycetes="#3435A9", Dothideomycetes="#FF0800", Eurotiomycetes="#00F000", Exobasidiomycetes="#00530C", Lecanoromycetes="#A5FF8B", Leotiomycetes="#4B9D9E", Microbotryomycetes="#FF1682", Sordariomycetes="#F000D9", Taphrinomycetes="#F2FF00", Tremellomycetes="#00E9E9", Unclassified="#888888")
classes <- sort(unique(sig$Class))
extra <- setdiff(classes,names(base_cols))
extra_cols <- if(length(extra)>0) setNames(hcl.colors(length(extra),"Dark 3"),extra) else character(0)
cols <- c(base_cols,extra_cols)[classes]
V(g)$Class <- sig$Class[match(V(g)$name,sig$ASV)]
V(g)$color <- unname(cols[V(g)$Class])
V(g)$color[V(g)$type] <- "white"
V(g)$size <- ifelse(V(g)$type,20,9)
V(g)$label <- ifelse(V(g)$type,V(g)$name,NA_character_)
V(g)$frame.color <- ifelse(V(g)$type,"black",NA)
stopifnot(!anyNA(V(g)$color))
set.seed(8046)
lay <- layout_with_fr(g,niter=999,weights=NA)
lay[,1] <- lay[,1]-mean(range(lay[,1]))
lay[,2] <- lay[,2]-mean(range(lay[,2]))
lay[,1] <- lay[,1]/max(abs(lay[,1]))*3.2
lay[,2] <- lay[,2]/max(abs(lay[,2]))*4.2
pdf("Figure2_FR.pdf",width=12,height=8,family="Times")
layout(matrix(c(1,2),nrow=1),widths=c(1.8,1.3))
par(mar=c(0.5,0,0.5,0))
plot(g, layout=lay, rescale=FALSE, vertex.size=V(g)$size*3.00, vertex.color=V(g)$color, vertex.frame.color=V(g)$frame.color, vertex.label=V(g)$label, vertex.label.family="Times", vertex.label.font=1, vertex.label.cex=1.1, vertex.label.color="black", edge.color="grey75", edge.width=1.2, edge.curved=0, asp=0, xlim=c(-4.5,4.5), ylim=c(-5.5,5.5), main="")
par(mar=c(0,0,0,0))
plot.new()
legend(x=0.15, y=0.5, xjust=0, yjust=0.5, legend=names(cols), col=cols, pch=16, pt.cex=2.023, cex=1.367, bty="n", y.intersp=1.517)
dev.off()
# Export analysis results
write.csv(results,"Figure2_IndVal_results.csv",row.names=FALSE)
write.csv(sig,"Figure2_significant_ASVs.csv",row.names=FALSE)

# =======================================================
# NMDS, PERMANOVA, and Bray-Curtis dissimilarity analysis
# =======================================================
install.packages(c("vegan", "permute", "ggplot2", "patchwork"))
library(vegan)
library(permute)
library(ggplot2)
library(patchwork)
# Import data
asv <- read.delim("d1_ASV_table.txt", row.names = 1, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE)
asv <- data.frame(t(asv))
metadata <- read.delim("d3_metadata_table.txt", sep = "\t", stringsAsFactors = FALSE, check.names = FALSE)
# Calculate Bray-Curtis dissimilarity
bray_dis <- vegdist(asv, method = "bray")
# NMDS analysis
set.seed(123)
nmds_dis <- metaMDS(bray_dis, k = 2, trymax = 200, trace = FALSE)
# Extract NMDS site scores
nmds_dis$stress
nmds_dis_site <- data.frame(scores(nmds_dis, display = "sites"))
nmds_dis_site$Samples <- rownames(nmds_dis_site)
nmds_dis_site$Group <- metadata$Group[match(nmds_dis_site$Samples, metadata$Samples)]
nmds_dis_site$PairID <- metadata$PairID[match(nmds_dis_site$Samples, metadata$Samples)]
# NMDS visualization
group_colors <- c(YL = "#0000FF", OL = "#FF0000")
Figure4A <- ggplot(nmds_dis_site, aes(NMDS1, NMDS2)) + stat_ellipse(aes(fill = Group, color = Group), geom = "polygon", type = "t", level = 0.95, alpha = 0.10, linewidth = 0.5, show.legend = FALSE) + geom_point(aes(color = Group), size = 2.5) + scale_color_manual(values = group_colors, breaks = c("YL", "OL")) + scale_fill_manual(values = group_colors) + annotate("text", x = -Inf, y = Inf, label = sprintf("Stress = %.2f", nmds_dis$stress), hjust = -0.08, vjust = 1.5, size = 4.5, family = "serif") + labs(x = "NMDS1", y = "NMDS2", color = NULL) + theme_classic(base_size = 14, base_family = "serif") + theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 0.75), axis.line = element_blank(), axis.text = element_text(color = "black", size = 12, family = "serif"), axis.title = element_text(color = "black", size = 14, family = "serif"), axis.ticks = element_line(color = "black", linewidth = 0.5), axis.ticks.length = grid::unit(1.35, "mm")) + theme(legend.position = c(0.98, 0.98), legend.justification = c(1, 1), legend.title = element_blank(), legend.background = element_blank(), legend.key = element_blank(), legend.key.size = grid::unit(3, "mm"), legend.spacing.y = grid::unit(0, "mm"), legend.text = element_text(size = 12, family = "serif")) + guides(color = guide_legend(override.aes = list(size = 2.5)))
# PERMANOVA analysis
set.seed(123)
perm <- how(nperm = 999, within = Within(type = "free"), blocks = metadata$PairID)
adonis_result <- adonis2(bray_dis ~ Group, data = metadata, permutations = perm)
adonis_result
# Calculate Bray-Curtis dissimilarities
dis <- as.matrix(bray_dis)
YL <- as.character(metadata$Samples[metadata$Group == "YL"])
OL <- as.character(metadata$Samples[metadata$Group == "OL"])
dis_YL <- as.vector(as.dist(dis[YL, YL]))
dis_OL <- as.vector(as.dist(dis[OL, OL]))
dat <- data.frame(dis = c(dis_YL, dis_OL), Group = factor(c(rep("YL", length(dis_YL)), rep("OL", length(dis_OL))), levels = c("YL", "OL")))
# Wilcoxon rank-sum test
wilcox_result <- wilcox.test(dis ~ Group, data = dat, paired = FALSE, exact = FALSE)
p_value <- wilcox_result$p.value
# Assign significance letters
letters <- if (p_value < 0.001) c("a", "b") else c("a", "a")
label_data <- data.frame(Group = factor(c("YL", "OL"), levels = c("YL", "OL")), Letter = letters, y = c(max(dis_YL), max(dis_OL)) + 0.04)
# Bray-Curtis dissimilarity boxplot
Figure4B <- ggplot(dat, aes(Group, dis)) + geom_boxplot(aes(fill = Group), width = 0.6, linewidth = 0.5) + scale_fill_manual(values = group_colors) + geom_text(data = label_data, aes(x = Group, y = y, label = Letter), inherit.aes = FALSE, size = 5, family = "serif") + coord_cartesian(ylim = c(0, max(1, label_data$y + 0.04))) + labs(x = NULL, y = "Bray-Curtis dissimilarity") + theme_classic(base_size = 14, base_family = "serif") + theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 0.75), axis.line = element_blank(), axis.text = element_text(color = "black", size = 12, family = "serif"), axis.title = element_text(color = "black", size = 14, family = "serif"), axis.ticks = element_line(color = "black", linewidth = 0.5), axis.ticks.length = grid::unit(1.35, "mm")) + theme(legend.position = "none")
# Combine NMDS and Bray-Curtis boxplot
Figure4_combined <- Figure4A + Figure4B + plot_layout(ncol = 2, widths = c(1.15, 1)) + plot_annotation(tag_levels = "A") & theme(plot.tag = element_text(size = 16, face = "bold", family = "serif"))
ggsave("Figure4.pdf", Figure4_combined, width = 11, height = 5, units = "in")

# ================================================================================
# Neutral community model (Figure 5A) and modified stochasticity ratio (Figure 5B)
# ================================================================================
install.packages(c("Hmisc", "minpack.lm", "ggplot2", "patchwork", "NST"))
library(Hmisc)
library(minpack.lm)
library(ggplot2)
library(patchwork)
library(NST)
# Data input and preparation
asv <- read.delim("d1_ASV_table.txt", row.names = 1, header = TRUE, sep = "\t", check.names = FALSE)
YL_cols <- grep("^YL[0-9]+$", colnames(asv))
OL_cols <- grep("^OL[0-9]+$", colnames(asv))
YL <- t(as.matrix(asv[, YL_cols, drop = FALSE]))
OL <- t(as.matrix(asv[, OL_cols, drop = FALSE]))
group <- read.delim("d3_metadata_table.txt", row.names = 1, check.names = FALSE)
asv_nst <- t(as.matrix(asv))
# Fit Sloan's neutral community model
fit_ncm <- function(spp) {
stopifnot(!anyNA(spp), all(is.finite(spp)), all(spp >= 0),all(spp == floor(spp)), all(rowSums(spp) > 0))
spp <- spp[, colSums(spp) > 0, drop = FALSE]
N <- mean(rowSums(spp))
p <- colMeans(spp) / N
freq <- colMeans(spp > 0)
d <- 1 / N
model <- minpack.lm::nlsLM(freq ~ pbeta(d, N * m * p, N * m * (1 - p), lower.tail = FALSE), start = list(m = 0.01), lower = c(m = 1e-10), upper = c(m = 1))
m <- unname(coef(model)["m"])
predicted <- pbeta(d, N * m * p, N * m * (1 - p), lower.tail = FALSE)
ci <- Hmisc::binconf(predicted * nrow(spp), nrow(spp), alpha = 0.05, method = "wilson")
R2 <- 1 - sum((freq - predicted)^2) / sum((freq - mean(freq))^2)
dat <- data.frame(ASV = names(p), p = as.numeric(p), freq = as.numeric(freq), predicted = as.numeric(predicted), Lower = ci[, "Lower"], Upper = ci[, "Upper"])
dat$Category <- ifelse(dat$freq > dat$Upper, "Above", ifelse(dat$freq < dat$Lower, "Below", "Neutral"))
dat <- dat[order(dat$p), ]
list(data = dat, m = m, R2 = R2, N = N, model = model)
}
ncm_YL <- fit_ncm(YL)
ncm_OL <- fit_ncm(OL)
# Plot Figure 5A
plot_ncm <- function(result, group) {
dat <- result$data
dat$logp <- log10(dat$p)
ggplot(dat, aes(x = logp)) + geom_line(aes(y = predicted), color = "blue", linewidth = 0.7) + geom_line(aes(y = Lower), color = "blue", linewidth = 0.6, linetype = "dashed") + geom_line(aes(y = Upper), color = "blue", linewidth = 0.6, linetype = "dashed") + geom_point(aes(y = freq, color = Category), size = 1.5) + scale_color_manual(values = c(Above = "blue", Below = "red", Neutral = "black")) + scale_y_continuous(breaks = seq(0, 1, 0.2)) + coord_cartesian(ylim = c(0, 1.02)) + annotate("text", x = -Inf, y = Inf, label = sprintf("R² = %.4f\nm = %.4f", result$R2, result$m), hjust = -0.05, vjust = 1.1, size = 4.5, family = "serif") + labs(title = group, x = "Mean relative abundance (log10)", y = "Occurrence frequency of ASVs") + theme_classic(base_size = 14, base_family = "serif") + theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 0.75), axis.line = element_blank(), axis.text = element_text(color = "black", size = 12), axis.ticks = element_line(color = "black", linewidth = 0.5), axis.ticks.length = grid::unit(1.35, "mm"), axis.title = element_text(size = 14, color = "black"), plot.title = element_text(hjust = 0.5, size = 16), legend.position = "none")
}
p_YL <- plot_ncm(ncm_YL, "YL") + labs(tag = "A") + theme(plot.tag = element_text(size = 16, face = "bold", family = "serif"))
p_OL <- plot_ncm(ncm_OL, "OL")
Figure5A <- (p_YL | p_OL) + plot_layout(ncol = 2, widths = c(1, 1), axis_titles = "collect")
# Calculate modified stochasticity ratio using abundance-based Jaccard (Ruzicka) distance
set.seed(123)
tnst <- tNST(comm = asv_nst, group = group, dist.method = "jaccard", null.model = "PF", rand = 1000, nworker = 1)
# Extract the values of modified stochasticity ratio
mst <- as.data.frame(tnst$index.pair.grp)
stopifnot(all(c("group", "name1", "name2","MST.ij.ruzicka") %in% names(mst)))
mst_data <- data.frame(Group = factor(mst$group, levels = c("YL", "OL")), Sample1 = mst$name1, Sample2 = mst$name2, MST = mst$MST.ij.ruzicka)
# Wilcoxon signed-rank test
sample_mst <- rbind(data.frame(Sample = mst_data$Sample1, MST = mst_data$MST), data.frame(Sample = mst_data$Sample2, MST = mst_data$MST))
sample_mst <- aggregate(MST ~ Sample, sample_mst, mean)
YL_mst <- sample_mst$MST[match(paste0("YL", 1:25), sample_mst$Sample)]
OL_mst <- sample_mst$MST[match(paste0("OL", 1:25), sample_mst$Sample)]
p_value <- wilcox.test(YL_mst, OL_mst, paired = TRUE, exact = FALSE)$p.value
# Assign significance letters
letters <- if (p_value < 0.05) {if (mean(YL_mst) > mean(OL_mst)) c("a", "b") else c("b", "a")} else {c("a", "a")}
label_data <- data.frame(Group = factor(c("YL", "OL"), levels = c("YL", "OL")), Letter = letters, y = c(max(mst_data$MST[mst_data$Group == "YL"]) + 0.05, max(mst_data$MST[mst_data$Group == "OL"]) + 0.05))
ylim_top <- max(1, ceiling(max(label_data$y + 0.05) / 0.25) * 0.25)
# Plot Figure 5B
mst_data$X <- ifelse(mst_data$Group == "YL", 1, 2)
label_data$X <- ifelse(label_data$Group == "YL", 1, 2)
Figure5B <- ggplot(mst_data, aes(x = X, y = MST, color = Group)) + geom_hline(yintercept = 0.5, linetype = "dashed", linewidth = 0.6) + geom_boxplot(aes(group = Group), width = 0.65, fill = NA, linewidth = 0.75, outlier.shape = NA) + geom_jitter(width = 0.12, height = 0, size = 1.5, alpha = 0.7) + geom_text(data = label_data, aes(x = X, y = y, label = Letter), inherit.aes = FALSE, family = "serif", size = 5, color = "black") + scale_color_manual(values = c(YL = "blue", OL = "red")) + scale_x_continuous(breaks = c(1, 2), labels = c("YL", "OL"), limits = c(0.25, 2.75), expand = expansion(mult = 0)) + scale_y_continuous(breaks = seq(0, ylim_top, 0.25), labels = function(x) sprintf("%.2f", x)) + coord_cartesian(ylim = c(0, ylim_top), expand = FALSE) + labs(tag = "B", x = NULL, y = "Modified stochasticity ratio (MST)") + theme_classic(base_size = 14, base_family = "serif") + theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 0.75), axis.line = element_blank(), axis.text = element_text(color = "black", size = 12), axis.ticks = element_line(color = "black", linewidth = 0.5), axis.ticks.length = grid::unit(1.35, "mm"), axis.title.y = element_text(size = 14), plot.tag = element_text(size = 16, face = "bold", family = "serif"), legend.position = "none")
ggsave("Figure5A.pdf", Figure5A, width = 10, height = 5, units = "in")
ggsave("Figure5B.pdf", Figure5B, width = 5, height = 5, units = "in")

# ===================================================================
# Co-occurrence networks and network topological properties(Figure 6)
# ===================================================================
install.packages(c("igraph", "psych", "ggplot2", "ggraph", "patchwork", "scales", "cowplot"))
library(igraph)
library(psych)
library(ggplot2)
library(ggraph)
library(patchwork)
library(scales)
library(cowplot)
library(grid)
# Data input and preparation
asv <- read.delim("d1_ASV_table.txt", header = TRUE, row.names = 1, sep = "\t", check.names = FALSE)
taxonomy <- read.delim("d2_taxonomy_table.txt", header = TRUE, row.names = 1, sep = "\t", check.names = FALSE)
YL_samples <- grep("^YL", colnames(asv), value = TRUE)
OL_samples <- grep("^OL", colnames(asv), value = TRUE)
asv_YL <- asv[, YL_samples, drop = FALSE]
asv_OL <- asv[, OL_samples, drop = FALSE]
# Construct co-occurrence networks
build_network <- function(otu, tax) {
otu <- as.matrix(otu)
storage.mode(otu) <- "numeric"
otu <- otu[, colSums(otu) > 0, drop = FALSE]
otu <- otu[rowSums(otu >= 2) >= 3, , drop = FALSE]
otu_rel <- sweep(otu, 2, colSums(otu), "/")
otu_rel <- otu_rel[rowMeans(otu_rel) > 0.0001, , drop = FALSE]
otu_rel <- otu_rel[apply(otu_rel, 1, sd) > 0, , drop = FALSE]
occor <- psych::corr.test(t(otu_rel), method = "spearman", adjust = "none")
r <- as.matrix(occor$r)
p_raw <- as.matrix(occor$p)
upper <- upper.tri(p_raw)
p_adj <- matrix(1, nrow(p_raw), ncol(p_raw), dimnames = dimnames(p_raw))
p_adj[upper] <- p.adjust(p_raw[upper], method = "BH")
p_adj[lower.tri(p_adj)] <- t(p_adj)[lower.tri(p_adj)]
r[p_adj >= 0.01 | abs(r) < 0.6] <- 0
r <- (r + t(r)) / 2
diag(r) <- 0
g <- igraph::graph_from_adjacency_matrix(r, mode = "undirected", weighted = TRUE, diag = FALSE)
g <- igraph::delete_vertices(g, igraph::V(g)[igraph::degree(g) == 0])
idx <- match(igraph::V(g)$name, rownames(tax))
igraph::V(g)$Class <- as.character(tax$Class[idx])
igraph::V(g)$Class[is.na(igraph::V(g)$Class) | igraph::V(g)$Class == ""] <- "Unclassified"
igraph::E(g)$correlation <- igraph::E(g)$weight
igraph::E(g)$sign <- ifelse(igraph::E(g)$correlation > 0, "Positive", "Negative")
igraph::V(g)$degree <- igraph::degree(g)
return(g)
}
g_YL <- build_network(asv_YL, taxonomy)
g_OL <- build_network(asv_OL, taxonomy)
# Plot co-occurrence networks
class_cols <- c(Agaricomycetes = "#950095", Agaricostilbomycetes = "#3435A9", Dothideomycetes = "#FF0800", Eurotiomycetes = "#00F000", Microbotryomycetes = "#FF1682", Pucciniomycetes = "#008080", Sordariomycetes = "#F000D9", Taphrinomycetes = "#F2FF00", Tremellomycetes = "#00E9E9", Ustilaginomycotina_Incertae_sedis = "#A640C5",Others = "#BDBDBD")
class_order <- c(sort(setdiff(names(class_cols), "Others")), "Others")
plot_network <- function(g, group_name) {
g_plot <- g
cls <- trimws(as.character(igraph::V(g_plot)$Class))
cls <- sub("^c__", "", cls)
cls[is.na(cls) | cls == "" | !cls %in% names(class_cols)] <- "Others"
igraph::V(g_plot)$Class_plot <- factor(cls, levels = class_order)
igraph::E(g_plot)$weight <- abs(igraph::E(g_plot)$correlation)
set.seed(123)
ggraph(g_plot, layout = "fr") + geom_edge_link(aes(edge_colour = sign, edge_width = abs(correlation)), alpha = 0.7, show.legend = FALSE) + scale_edge_width(range = c(0.5, 2), guide = "none") + scale_edge_colour_manual(values = c(Positive = "#FF0000", Negative = "#0000FF"), guide = "none") + geom_node_point(aes(colour = Class_plot, size = degree), alpha = 0.9) + scale_colour_manual(values = class_cols, limits = class_order, breaks = class_order, drop = FALSE, name = NULL, na.translate = FALSE) + scale_size_continuous(range = c(1.2, 5), guide = "none") + ggtitle(group_name) + theme_void(base_family = "Times New Roman") + theme(plot.title = element_text(family = "Times New Roman", size = 16, face = "plain", hjust = 0.5, colour = "black"),legend.position = "none", plot.margin = margin(5, 5, 0, 5))
}
p_YL <- plot_network(g_YL, "YL")
p_OL <- plot_network(g_OL, "OL")
legend_df <- data.frame(Class = factor(class_order, levels = class_order), x = 1, y = seq_along(class_order))
legend_plot <- ggplot(legend_df, aes(x = x, y = y, colour = Class)) + geom_point(shape = 16, size = 3.5, alpha = 1) + scale_colour_manual(values = class_cols, limits = class_order, breaks = class_order, drop = FALSE, name = NULL) + guides(colour = guide_legend(nrow = 3, byrow = TRUE, override.aes = list(shape = 16, size = 3.5, alpha = 1))) + theme_void(base_family = "Times New Roman") + theme(legend.position = "bottom", legend.direction = "horizontal", legend.text = element_text(family = "Times New Roman", size = 10, face = "plain", colour = "black"), legend.title = element_blank(), legend.key = element_blank(), legend.background = element_blank(), legend.box = "horizontal", legend.spacing.x = unit(0.25, "cm"))
legend_class <- cowplot::get_legend(legend_plot)
network_combined <- cowplot::plot_grid(p_YL, p_OL, ncol = 2, align = "hv", labels = NULL)
network_with_A <- cowplot::ggdraw(network_combined) + cowplot::draw_label( "A", x = 0.015, y = 0.98, hjust = 0, vjust = 1, fontfamily = "Times New Roman", fontface = "bold", size = 14, colour = "black")
Figure6A <- cowplot::plot_grid(network_with_A, legend_class, ncol = 1, rel_heights = c(1, 0.22))
# Calculate network topological metrics
get_metrics <- function(g, group_name) {
node_number <- igraph::vcount(g)
edge_number <- igraph::ecount(g)
bt <- igraph::betweenness(g, directed = FALSE, weights = NA, normalized = FALSE)
mean_betweenness <- if (length(bt) > 0) {mean(bt)} else {NA_real_}
assortativity <- if (node_number > 1 && edge_number > 0) {igraph::assortativity_degree(g, directed = FALSE)} else {NA_real_}
positive_n <- sum(igraph::E(g)$correlation > 0, na.rm = TRUE)
negative_n <- sum(igraph::E(g)$correlation < 0, na.rm = TRUE)
negative_positive_ratio <- if (positive_n > 0) {negative_n / positive_n} else {NA_real_}
data.frame(Group = group_name, Nodes = node_number, Edges = edge_number, Betweenness = mean_betweenness, Assortativity = assortativity, Positive_edges = positive_n, Negative_edges = negative_n, Negative_Positive_Ratio = negative_positive_ratio)
}
metrics <- rbind(get_metrics(g_YL, "YL"), get_metrics(g_OL, "OL"))
# Define a common plotting function
bar_cols <- c(YL = "#0000FF", OL = "#FF0000")
plot_metric <- function(metric, y_label, y_max, y_step, digits, data = metrics, error = FALSE, letter_data = NULL, error_type = "SD") {
df <- data.frame(Group = factor(data$Group, levels = c("YL", "OL")), Value = data[[metric]])
p <- ggplot(df, aes(x = Group, y = Value, fill = Group)) + geom_col(width = 0.65, colour = "black", linewidth = 0.5, na.rm = TRUE) + scale_fill_manual(values = bar_cols, guide = "none") + scale_y_continuous(breaks = seq(0, y_max, by = y_step), labels = function(x) {sprintf(paste0("%.", digits, "f"), x)}, expand = expansion(mult = c(0, 0))) + coord_cartesian(ylim = c(0, y_max), clip = "off") + labs(x = NULL, y = y_label) + theme_classic(base_family = "Times New Roman", base_size = 11) + theme(panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5), axis.line = element_blank(), axis.ticks = element_line(colour = "black", linewidth = 0.5), axis.ticks.length = grid::unit(3, "pt"), axis.text = element_text(family = "Times New Roman", size = 10, colour = "black"), axis.title = element_text(family = "Times New Roman", size = 11, colour = "black"), axis.title.y = element_text(margin = margin(r = 5)), legend.position = "none", plot.margin = margin(5, 6, 5, 6))
if (error) {error_df <- data.frame(Group = factor(data$Group, levels = c("YL", "OL")), Lower = pmax(0, data$Mean - data[[error_type]]), Upper = data$Mean + data[[error_type]])
p <- p + geom_errorbar(data = error_df, aes(x = Group, ymin = Lower, ymax = Upper), inherit.aes = FALSE, width = 0.16, linewidth = 0.5, colour = "black", na.rm = TRUE)}
if (!is.null(letter_data)) {p <- p + geom_text(data = letter_data, aes(x = Group, y = y, label = Letter), inherit.aes = FALSE, family = "Times New Roman", size = 4, colour = "black")}
return(p)
}
# Plot network topological metrics
p1 <- plot_metric(metric = "Nodes", y_label = "Node number", y_max = 100, y_step = 20, digits = 0)
p2 <- plot_metric(metric = "Edges", y_label = "Edge number", y_max = 100, y_step = 20, digits = 0)
p3 <- plot_metric(metric = "Betweenness", y_label = "Betweenness", y_max = 1.0, y_step = 0.2, digits = 1)
p4 <- plot_metric(metric = "Assortativity", y_label = "Assortativity", y_max = 1.0, y_step = 0.2, digits = 1)
p5 <- plot_metric(metric = "Negative_Positive_Ratio", y_label = paste0("Ratio of negative correlation and\n", "positive correlation"), y_max = 0.1, y_step = 0.02, digits = 2)
# Calculate network robustness
set.seed(123)
make_unweighted_matrix <- function(g) {
m <- as.matrix(igraph::as_adjacency_matrix(g, attr = NULL, sparse = FALSE))
m[m != 0] <- 1
diag(m) <- 0
return(m)
}
network_YL <- make_unweighted_matrix(g_YL)
network_OL <- make_unweighted_matrix(g_OL)
rand.remov.once <- function(netRaw, rm.percent = 0.5) {
n <- nrow(netRaw)
id.rm <- sample(seq_len(n), size = round(n * rm.percent), replace = FALSE)
net.Raw <- netRaw
net.Raw[id.rm, ] <- 0
net.Raw[, id.rm] <- 0
sum(rowSums(net.Raw) > 0) / n
}
nperm <- 1000
robustness_YL <- replicate(nperm, rand.remov.once(network_YL, 0.5))
robustness_OL <- replicate(nperm, rand.remov.once(network_OL, 0.5))
# Compare robustness using t-test
test_robustness <- t.test(robustness_YL, robustness_OL, alternative = "two.sided", paired = FALSE, var.equal = FALSE)
p_value_robustness <- test_robustness$p.value
if (p_value_robustness < 0.05) {
letter_YL <- "b"
letter_OL <- "a"
} else {
letter_YL <- "a"
letter_OL <- "a"
}
robustness_summary <- data.frame(Group = factor(c("YL", "OL"), levels = c("YL", "OL")), Mean = c(mean(robustness_YL), mean(robustness_OL)), SD = c(sd(robustness_YL), sd(robustness_OL)), SE = c(sd(robustness_YL) / sqrt(nperm), sd(robustness_OL) / sqrt(nperm)), Letter = c(letter_YL, letter_OL))
letter_df <- robustness_summary
# Plot network robustness
letter_df$y <- pmin(0.385, letter_df$Mean + letter_df$SD + 0.015)
p6 <- plot_metric(metric = "Mean", y_label = "Robustness", y_max = 0.4, y_step = 0.1, digits = 1, data = robustness_summary, error = TRUE, letter_data = letter_df, error_type = "SD")
# Calculate network vulnerability
global_efficiency <- function(g) {
n <- igraph::vcount(g)
if (n < 2) {return(0)}
d <- igraph::distances(g, mode = "all", weights = NA)
inv_d <- 1 / d
inv_d[!is.finite(inv_d)] <- 0
diag(inv_d) <- 0
sum(inv_d) / (n * (n - 1))
}
calculate_vulnerability <- function(g) {
g <- igraph::as.undirected(g, mode = "collapse")
g <- igraph::delete_vertices(g, igraph::V(g)[igraph::degree(g) == 0])
n <- igraph::vcount(g)
E_original <- global_efficiency(g)
node_vulnerability <- vapply(seq_len(n),
function(i) {
g_removed <- igraph::delete_vertices(g, i)
(E_original - global_efficiency(g_removed)) / E_original}, numeric(1))
names(node_vulnerability) <- igraph::V(g)$name
return(node_vulnerability)
}
vulnerability_YL <- calculate_vulnerability(g_YL)
vulnerability_OL <- calculate_vulnerability(g_OL)
vulnerability_summary <- data.frame(Group = factor(c("YL", "OL"), levels = c("YL", "OL")), Maximum_vulnerability = c(max(vulnerability_YL), max(vulnerability_OL)))
# Plot network vulnerability
p7 <- plot_metric(metric = "Maximum_vulnerability", y_label = "Vulnerability", y_max = 0.09, y_step = 0.02, digits = 2, data = vulnerability_summary)
#  Assemble Figure 6B and Figure 6C
Figure6B_raw <- (p1 | p2 | p3 | p4) + plot_layout(ncol = 4)
Figure6B <- cowplot::ggdraw(Figure6B_raw) + cowplot::draw_label("B", x = 0.005, y = 0.99, hjust = 0, vjust = 1, fontfamily = "Times New Roman", fontface = "bold", size = 14, colour = "black")
Figure6C_raw <- (p6 | p7 | p5 | plot_spacer()) + plot_layout(ncol = 4)
Figure6C <- cowplot::ggdraw(Figure6C_raw) + cowplot::draw_label("C", x = 0.005, y = 0.99, hjust = 0, vjust = 1, fontfamily = "Times New Roman", fontface = "bold", size = 14, colour = "black")

# ============================================================
# Random forest and Spearman correlation analysis(Figure 7)
# ============================================================
install.packages(c("randomForest", "ggplot2", "patchwork"))
library(randomForest)
library(ggplot2)
library(patchwork)
# Data input and preparation
asv <- read.delim("d1_ASV_table.txt", header = TRUE, row.names = 1, check.names = FALSE)
taxonomy <- read.delim("d2_taxonomy_table.txt", header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)
traits <- read.delim("d4_traits_table.txt", header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)
asv <- asv[, traits$Samples, drop = FALSE]
rel <- sweep(asv, 2, colSums(asv), "/") * 100
# Random forest model construction
trees <- unique(as.character(traits$PairID))
set.seed(123)
train_trees <- sample(trees, floor(0.8 * length(trees)))
train_idx <- as.character(traits$PairID) %in% train_trees
rel_rf <- rel[apply(rel[, train_idx, drop = FALSE], 1, sd) > 0,, drop = FALSE]
rf_data <- data.frame(Group = factor(traits$Group, levels = c("YL", "OL")), t(rel_rf), check.names = FALSE)
rf_train <- rf_data[train_idx, , drop = FALSE]
set.seed(123)
rf_model <- randomForest(Group ~ ., data = rf_train, ntree = 1000, importance = TRUE)
# Extract and annotate Top 20 important ASVs
imp <- importance(rf_model, type = 1, scale = FALSE)
imp_df <- data.frame(ASV = rownames(imp), importance = imp[, "MeanDecreaseAccuracy"])
top20 <- head(imp_df[order(imp_df$importance, decreasing = TRUE), ], 20)
top20$Class_original <- sub("^c__", "", as.character(taxonomy$Class[match(top20$ASV, taxonomy[["ASV ID"]])]))
top20$Class_original[is.na(top20$Class_original) | top20$Class_original == ""] <- "Unclassified"
class_cols <- c(Agaricomycetes = "#950095", Agaricostilbomycetes = "#3435A9", Dothideomycetes = "#FF0800", Eurotiomycetes = "#00F000", Microbotryomycetes = "#FF1682", Pucciniomycetes = "#008080", Sordariomycetes = "#F000D9", Taphrinomycetes = "#F2FF00", Tremellomycetes = "#00E9E9", Ustilaginomycotina_Incertae_sedis = "#A640C5", Others = "#BDBDBD")
top20$Class <- ifelse(top20$Class_original %in% names(class_cols), top20$Class_original, "Others")
top20$label <- paste0(top20$ASV, ": ", top20$Class_original)
top_ids <- as.character(top20$ASV)
asv_order <- rev(top_ids)
top20$ASV <- factor(top20$ASV, levels = asv_order)
label_map <- setNames(top20$label, as.character(top20$ASV))
# Visualization of random forest variable importance (Figure7A)
Figure7A <- ggplot(top20, aes(x = importance, y = ASV, fill = Class)) + geom_col(width = 0.55) + scale_fill_manual(values = class_cols) + scale_y_discrete(limits = asv_order, labels = label_map) + labs(x = "Mean decrease in accuracy", y = NULL, title = "A") + theme_bw(base_size = 10, base_family = "Times New Roman") + theme(panel.grid = element_blank(), panel.border = element_rect(color = "black", fill = NA, linewidth = 0.9), axis.line = element_blank(), legend.position = "none", axis.text = element_text(size = 10, color = "black"), axis.title.x = element_text(size = 10, color = "black"), plot.title = element_text(face = "bold", size = 12, hjust = 0), plot.margin = margin(5, 5, 5, 5))
# Spearman correlation analysis
trait_labels <- c("LCC", "LNC", "LPC", "Chl", "SOD", "CAT", "PRO", "MDA", "SPC")
trait_data <- data.frame(lapply(traits[, 2:10], as.numeric), check.names = FALSE)
names(trait_data) <- trait_labels
cor_results <- expand.grid(ASV = top_ids, Trait = trait_labels, stringsAsFactors = FALSE)
cor_results$r <- NA_real_
cor_results$p <- NA_real_
for (i in seq_len(nrow(cor_results))) {
x <- as.numeric(rel[cor_results$ASV[i], ])
y <- trait_data[[cor_results$Trait[i]]]
if (sd(x) > 0 && sd(y) > 0) {
test <- suppressWarnings(cor.test(x, y, method = "spearman", exact = FALSE))
cor_results$r[i] <- unname(test$estimate)
cor_results$p[i] <- test$p.value}
}
cor_results$p_adj <- p.adjust(cor_results$p, method = "BH")
cor_results$stars <- ifelse(is.na(cor_results$p_adj), "", ifelse(cor_results$p_adj < 0.001, "***", ifelse(cor_results$p_adj < 0.01, "**", ifelse(cor_results$p_adj < 0.05, "*", ""))))
cor_results$ASV <- factor(cor_results$ASV, levels = asv_order)
cor_results$Trait <- factor(cor_results$Trait, levels = trait_labels)
 # Visualization of Spearman correlations(Figure7B)
Figure7B <- ggplot(cor_results, aes(x = Trait, y = ASV, fill = r)) + geom_tile(color = "grey80", linewidth = 0.35) + geom_text(aes(label = stars), size = 10.5 / ggplot2::.pt, family = "Times New Roman") + scale_fill_gradient2(low = "#9C7A35", mid = "white", high = "#24873B", midpoint = 0, limits = c(-1, 1), name = "Spearman r", na.value = "grey90") + scale_y_discrete(limits = asv_order, labels = label_map) + scale_x_discrete(limits = trait_labels) + labs(x = NULL, y = NULL, title = "B") + theme_minimal(base_size = 10.5, base_family = "Times New Roman") + theme(panel.grid = element_blank(), axis.text.x = element_text(size = 10.5, angle = 90, vjust = 0.5, hjust = 1, color = "black"), axis.text.y = element_text(size = 10.5, color = "black"), plot.title = element_text(face = "bold", size = 12, hjust = 0), plot.margin = margin(5, 5, 5, 5))
Figure7 <- Figure7A / Figure7B + patchwork::plot_layout(ncol = 1, heights = c(1, 1.1))
  
# ============================================
# Piecewise Structural Equation Modeling (SEM)
# ============================================  
install.packages(c("lmerTest", "piecewiseSEM", "performance", "ggplot2", "patchwork"))
library(lmerTest)
library(piecewiseSEM)
library(performance)
library(ggplot2)
library(patchwork)
# Data input and preparation  
dat <- read.delim("d5_sem_table.txt", check.names = FALSE, stringsAsFactors = FALSE)
traits <- c("Carboncontent", "Chlorophyllcontent", "SODactivity", "MDAcontent")
responses <- c("Diversity", "Community")
vars <- c("Leafage", traits, responses)
dat$Leafage <- match(dat$Leafage, c("YL", "OL")) - 1
dat$PairID <- factor(dat$PairID)
dat[vars] <- as.data.frame(scale(dat[vars]))
# Piecewise SEM
models <- lapply(c(traits, responses), function(y) {
predictors <- if (y %in% traits) "Leafage" else c("Leafage", traits)
lmer(as.formula(paste(y, "~", paste(predictors, collapse = " + "), "+ (1 | PairID)")), data = dat, REML = FALSE)
})
names(models) <- c(traits, responses)
sem <- do.call(psem, unname(models))
sem_result <- summary(sem)
sem_result$AIC <- AIC(sem, AIC.type = "dsep")
print(sem_result)
paths <- do.call(rbind, lapply(names(models), function(y) {
x <- summary(models[[y]])$coefficients
data.frame(Predictor = rownames(x)[-1], Response = y, Estimate = x[-1, "Estimate"], SE = x[-1, "Std. Error"], P = x[-1, "Pr(>|t|)"], row.names = NULL)
}))
r2 <- sapply(models, function(m) r2_nakagawa(m, tolerance = 0)$R2_marginal)
# Direct, indirect and total effects plots (Figure8B-C)
coef_path <- function(x, y) paths$Estimate[paths$Predictor == x & paths$Response == y]
effects <- do.call(rbind, lapply(responses, function(y) {
p <- c("Leafage", traits)
direct <- sapply(p, coef_path, y = y)
indirect <- c(sum(sapply(traits, function(t) coef_path("Leafage", t) * coef_path(t, y))), rep(0, length(traits)))
data.frame(Response = y, Predictor = p, Direct = unname(direct), Indirect = indirect, Total = unname(direct + indirect))
}))
cols <- c("Direct effect" = "#0000FF", "Indirect effect" = "#008000", "Total effect" = "#FF0000")
plot_effect <- function(response, tag) {
d <- subset(effects, Response == response)
d <- d[match(c("Leafage", traits), d$Predictor), ]
d$Predictor <- factor(d$Predictor, levels = c("Leafage", traits), labels = c("Leaf age", "LCC", "Chl", "SOD", "MDA"))
d <- reshape(d, varying = c("Direct", "Indirect", "Total"), v.names = "Effect", timevar = "Type", times = names(cols), direction = "long")
d$Type <- factor(d$Type, levels = names(cols))
ggplot(d, aes(Predictor, Effect, fill = Type)) +
geom_hline(yintercept = 0, linetype = "dotted") + geom_col(position = position_dodge(0.8), width = 0.75) + scale_fill_manual(values = cols) + guides(fill = guide_legend(title = NULL, keywidth = grid::unit(0.3, "cm"), keyheight = grid::unit(0.3, "cm"))) + scale_y_continuous(breaks = seq(-0.5, 1, 0.5), expand = expansion(mult = 0)) + coord_cartesian(ylim = c(-0.5, 1)) + labs(x = NULL, y = "Standardized effect size", tag = tag) + theme_bw(base_size = 11, base_family = "Times New Roman") + theme(panel.grid = element_blank(), panel.border = element_rect(color = "black", fill = NA), axis.text = element_text(size = 10, color = "black"), axis.title.y = element_text(size = 12), axis.ticks = element_line(color = "black"), axis.ticks.length = grid::unit(3, "pt"), plot.tag = element_text(face = "bold", size = 13), plot.tag.position = c(0.02, 0.98), legend.position = c(0.98, 0.98), legend.justification = c(1, 1), legend.background = element_blank(), legend.key = element_blank(), legend.text = element_text(size = 10), legend.spacing.y = grid::unit(0.05, "cm"), plot.margin = margin(5, 10, 5, 5))
}
Figure8B <- plot_effect("Diversity", "B")
Figure8C <- plot_effect("Community", "C")
Figure8BC <- Figure8B + Figure8C + plot_layout(ncol = 2)


# ===============================================================================================================================================
# Statistical analysis of leaf traits (Table S2)
# A total of 25 pairs of young (YL) and old (OL) leaf samples (YL1–OL1 to YL25–OL25; 50 samples in total) were included in the analysis.
# Paired t-tests for Carbon content,Nitrogen content,Phosphorus content,Chlorophyll content,SOD activity,CAT activity,Proline content,MDA content
# Wilcoxon signed-rank test for Soluble protein content
# ===============================================================================================================================================
# Import data
traits <- read.delim("d4_traits_table.txt", header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)
# Extract paired samples
YL <- traits[traits$Group == "YL", ]
OL <- traits[traits$Group == "OL", ]
# Define leaf traits
indices <- c("Carbon content", "Nitrogen content", "Phosphorus content",
             "Chlorophyll content", "SOD activity", "CAT activity",
             "Proline content", "MDA content", "Soluble protein content")
# Shapiro-Wilk normality tests
normality_results <- data.frame(Index=indices,P_value=sapply(indices,function(i)
shapiro.test(YL[[i]]-OL[[i]])$p.value))
# Paired t-tests and Wilcoxon signed-rank test
test_results <- data.frame(Index=indices,P_value=sapply(indices,function(i){
if(i=="Soluble protein content")
wilcox.test(YL[[i]],OL[[i]],paired=TRUE,exact=FALSE)$p.value
else t.test(YL[[i]],OL[[i]],paired=TRUE)$p.value
}))
# Benjamini-Hochberg correction
test_results$P_adjusted <- p.adjust(test_results$P_value,method="BH")
# Calculate mean, SE and significance letters
summary_results <- do.call(rbind,lapply(indices,function(i){
m <- c(mean(YL[[i]]),mean(OL[[i]]))
se <- c(sd(YL[[i]]),sd(OL[[i]]))/sqrt(nrow(YL))
p <- test_results$P_value[test_results$Index==i]
letter <- if(p>=0.05) c("a","a") else if(m[1]>m[2]) c("a","b") else c("b","a")
data.frame(Index=i,Group=c("YL","OL"),Mean=m,SE=se,Letter=letter)
}))
# Export results
write.csv(summary_results, "d4_traits_mean_SE_results.csv", row.names = FALSE)