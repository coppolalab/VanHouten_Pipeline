#For WGCNA
library(WGCNA)
library(flashClust)
ALLOW_WGCNA_THREADS = 8
allowWGCNAThreads() #For baseline processing
library(limma)
library(sva)
library(peer)

#Functional programming
library(magrittr)
library(purrr)
library(functional)
library(lambda.r)

#Data arrangement
library(reshape2)
library(plyr)
library(dplyr)
library(tidyr)

#String operations
library(stringr)

#Plotting
library(ggplot2)
library(extrafont)
library(Cairo)

#Reading and writing tables
library(readr)
library(openxlsx)

gen.pcaplot <- function(filename, dataset, facet.bool, size.height, size.width)
{
    colnames(dataset)[2] <- "Module"
    dataset$Module %<>% str_replace("ME", "") 
    p <- ggplot(dataset, aes(x = as.numeric(x), y = value, fill = Module, color = Module)) + geom_point()
    p <- p + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
    p <- p + theme(axis.ticks.x = element_blank()) + ylab("First Principal Component")
    p <- p + theme(axis.title.x = element_blank()) + scale_x_continuous(as.numeric(unique(dataset$x)))
    p <- p + scale_color_manual(values = sort(unique(dataset$Module)))
    if (facet.bool == TRUE) {
        p <- p + facet_wrap(~ Module)
        p <- p + theme(legend.position = "none")
    } 
    CairoPDF(filename, height = size.height, width = size.width)
    print(p)
    dev.off()
}

gen.heatmap <- function(dataset, ME.genes)
{
    color <- as.character(unique(dataset$module.colors))
    dataset %<>% select(-module.colors) %>% scale
    max.dataset <- max(abs(dataset))
    CairoPDF(paste("./modules/", color, sep = ""), width = 21, height = 12)
    par(mar = c(3.5,3,2,3))
    par(oma = c(4,0,2,0))
    plotMat(dataset, zlim = c(-max.dataset, max.dataset), main = paste(color, " (", nrow(dataset), ")", sep = ""))

    ME.genes.plot <- select(ME.genes, Sample.ID, matches(color))
    p <- ggplot(ME.genes.plot, aes_string(x = "Sample.ID", y = color))
    p <- p + geom_bar(stat = "identity") + xlab("Eigengene Expression")#+ ylim(c(-6, 16)) 
    p <- p + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
    p <- p + theme(axis.text.x = element_text(angle = 90, size = 2))  
    p <- p + theme(axis.title.x = element_blank(), axis.ticks.x = element_blank())
    print(p)
    dev.off()
}

gen.cor <- function(dataset, trait.df)
{
    dataset %<>% mutate(SampleID = rownames(dataset))
    sample.key <- paste(trait.df$SampleID, collapse = "|")
    dataset.reduce <- filter(dataset, grepl(sample.key, SampleID))
    module.trait.cor <- cor(select(dataset.reduce, -SampleID), select(trait.df, -SampleID), use = "p")
    module.trait.cor.pval <- corPvalueStudent(module.trait.cor, nrow(dataset.reduce))
    colnames(module.trait.cor.pval) %<>% paste(".p.value", sep = "") 
    return(cbind(module.trait.cor, module.trait.cor.pval))
}

gen.text.heatmap <- function(cor.dataset, text.matrix, x.names, y.names, maintitle, filename)
{
    width.dynamic = 3 + (1 * ncol(text.matrix))
    CairoPDF(filename, width = width.dynamic, height = 10)
    par(mar = c(8, 8, 3, 3))
    labeledHeatmap(Matrix = cor.dataset, xLabels = x.names, yLabels = y.names, ySymbols = y.names, yColorLabels = TRUE, colors = greenWhiteRed(50), textMatrix = text.matrix, setStdMargins = F, cex.text = 0.5, zlim = c(-1,1), main = maintitle)
    dev.off()
}

#Peer analysis
gen.peer <- function(num.factors, intensities, use.covariates, covariates)
{
    model = PEER()
    PEER_setNk(model,num.factors)
    PEER_setPhenoMean(model, as.matrix(t(intensities)))
    PEER_setAdd_mean(model, TRUE)
    if (use.covariates == TRUE)
    {
        PEER_setCovariates(model, as.matrix(covariates))
    }
    PEER_setNmax_iterations(model, 1000)
    PEER_update(model)
    residuals.PEER = t(PEER_getResiduals(model))
    rownames(residuals.PEER) = rownames(intensities)
    colnames(residuals.PEER) = colnames(intensities)

    write_csv(data.frame(residuals.PEER), path = paste("residuals_", num.factors, sep = "", ".csv"))
    write_csv(data.frame(PEER_getX(model)), path = paste("factor_", num.factors, sep = "", ".csv"))
    write_csv(data.frame(PEER_getW(model)), path = paste("weight_", num.factors, sep = "", ".csv"))
    write_csv(data.frame(PEER_getAlpha(model)), path = paste("precision_", num.factors, sep = "", ".csv"))

    CairoPDF(file = paste("model", num.factors, ".pdf", sep = ""), width = 10, height = 10)
    PEER_plotModel(model)
    dev.off()

    CairoPDF(file = paste("precision_", num.factors, ".pdf", sep = ""), width = 10, height = 10)
    plot(PEER_getAlpha(model), col = "red", lwd = 4, main = paste("precision", num.factors, "factor", sep = " "))
    dev.off()
}

load(file = "../baseline/save/intensities.avg.norm.rda")
load(file = "../baseline/save/intensities1.norm.rda")
load(file = "../baseline/save/targets.final.rda")

targets.gaa <- filter(targets.final, Status == "Patient")
intensities1.gaa <- select(data.frame(intensities1.norm), one_of(targets.gaa$SampleID)) %>% as.matrix %>% normalizeBetweenArrays

model.sex.full <- model.matrix(~ 0 + factor(targets.gaa$Sex))[,-1]
targets.gaa$Age %<>% as.numeric 
draw.age.cov <- log2(targets.gaa$Age)
targets.gaa$Onset %<>% as.numeric

covariates.full <- cbind(Sex = model.sex.full, Age = draw.age.cov, GAA1 = log2(targets.gaa$GAA1), Onset = log2(targets.gaa$Onset))
gen.peer(5, intensities1.gaa, TRUE, covariates.full)

PEER.weights <- read_csv("./weight_5.csv") %>% select(-(X1:X5))
PEER.weights.sums <- colSums(abs(PEER.weights)) %>% data.frame
PEER.weights.sums$Factor <- 1:nrow(PEER.weights.sums)
colnames(PEER.weights.sums)[1] <- "Weight"

p <- ggplot(PEER.weights.sums, aes(x = factor(Factor), y = as.numeric(Weight), group = 1)) + geom_line(color = "blue") 
p <- p + theme_bw() + xlab("Factor") + ylab("Weigth")
CairoPDF("./PEER_weights", height = 4, width = 6)
print(p)
dev.off()

intensities.PEER <- t(intensities1.gaa)
save(intensities.PEER, file = "./save/intensities.peer.rda")
intensities1.gaa %<>% t

#Calculate scale free topology measures for different values of power adjacency function
powers <- c(c(1:10), seq(from = 12, to = 39, by = 2))

sft.PEER <- pickSoftThreshold(intensities1.gaa, powerVector = powers, verbose = 5, networkType = "signed")
sft.PEER.df <- sft.PEER$fitIndices
save(sft.PEER, file = "./save/sft.peer.rda")

sft.PEER.df$multiplied <- sft.PEER.df$SFT.R.sq * -sign(sft.PEER.df$slope)
p <- ggplot(sft.PEER.df, aes(x = Power,  y = multiplied, label = Power)) + geom_point() + geom_text(vjust = -0.6, size = 4, col = "red")
p <- p + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + geom_hline(aes(yintercept = 0.9))
p <- p + xlab("Soft Threshold") + ylab("Scale Free Topology Model Fit, signed R^2") + ggtitle("Scale Independence")
CairoPDF(file = "./scaleindependence", width = 6, height = 6)
print(p)
dev.off()

p <- ggplot(sft.PEER.df, aes(x = Power,  y = mean.k.)) + geom_point() 
p <- p + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
p <- p + xlab("Soft Threshold") + ylab("Mean Connectivity") + ggtitle("Mean Connectivity")
CairoPDF(file = "./meanconnectivity", width = 6, height = 6)
print(p)
dev.off()

softPower <- 20
adjacency.PEER <- adjacency(intensities.PEER, power = softPower, type = "signed")
save(adjacency.PEER, file = "./save/adjacency.PEER.rda")

TOM.PEER <- TOMsimilarity(adjacency.PEER)
dissimilarity.TOM <- 1 - TOM.PEER
save(TOM.PEER, file = "./save/tom.PEER.rda")
save(dissimilarity.TOM, file = "./save/dissimilarity.TOM.rda")

geneTree = flashClust(as.dist(dissimilarity.TOM), method = "average")
save(geneTree, file = "./save/gene.tree.rda")

CairoPDF(file = "./1-genecluster", height = 10, width = 15)
plot(geneTree, xlab = "", sub = "", main = "Gene clustering on TOM-based dissimilarity", labels = FALSE, hang = 0.04)
dev.off()

min.module.size <- 50

#Identify modules using dynamic tree cutting with hybrid clustering
dynamic.modules <- cutreeDynamic(dendro = geneTree, method = "hybrid", distM = dissimilarity.TOM, cutHeight = 0.995, deepSplit = 2, pamRespectsDendro = FALSE, minClusterSize = min.module.size, verbose = 2)
dynamic.colors <- labels2colors(dynamic.modules)
save(dynamic.colors, file = "./save/dynamic.colors.rda")

CairoPDF(file = "./1-gene_dendrogram_and_module_colors_min50", height = 10, width = 15)
plotDendroAndColors(geneTree, dynamic.colors, "Dynamic Tree Cut", dendroLabels = FALSE, hang = 0.03, addGuide = TRUE, guideHang = 0.05)
dev.off()

#Calculate module eigengenes
ME.list <- moduleEigengenes(intensities.PEER, colors = dynamic.colors, softPower = softPower, nPC = 1)
ME.genes <- ME.list$eigengenes
MEDiss <- 1 - cor(ME.genes)
METree <- flashClust(as.dist(MEDiss), method = "average")
save(METree, file = "./save/me.tree.rda")

CairoPDF(file = "./3-module_eigengene_clustering_min50", height = 10, width = 15)
plot(METree, xlab = "", sub = "", main = "")
dev.off()

#Check if any modules are too similar and merge them.  Possibly not working.
ME.dissimilarity.threshold <- 0.20
merge.all <- mergeCloseModules(intensities.PEER, dynamic.colors, cutHeight = ME.dissimilarity.threshold, verbose = 3) #PC analysis may be failing because of Intel MKL Lapack routine bug.  Test with openBLAS in R compiled with gcc.
merged.colors <- merge.all$colors
merged.genes <- merge.all$newMEs

CairoPDF("4-module_eigengene_clustering_min50", height = 10, width = 15)
plotDendroAndColors(geneTree, cbind(dynamic.colors, merged.colors), c("Dynamic Tree Cut", "Merged dynamic"), dendroLabels = FALSE, hang = 0.03, addGuide = TRUE, guideHang = 0.05, main = "")
dev.off()

#Use merged eigengenes 
module.colors <- merged.colors
save(module.colors, file = "./save/module.colors.rda")
color.order <- c("grey", standardColors(50))
modules.labels <- match(module.colors, color.order)
save(modules.labels, file = "./save/modules.labels.rda")
ME.genes <- merged.genes
save(ME.genes, file = "./save/me.genes.rda")

all.degrees <- intramodularConnectivity(adjacency.PEER, module.colors)
expr.data.t <- mutate(data.frame(t(intensities.PEER)), ProbeName = colnames(intensities.PEER)) 
gene.info.join <- join(expr.data.t, intensities.avg.norm$genes)
gene.info <- select(gene.info.join, ProbeName:Description) %>% mutate(module.colors = module.colors, mean.count = apply(intensities.PEER, 2, mean)) %>% data.frame(all.degrees)

CairoPDF("5-eigengenes", height = 10, width = 18)
par(cex = 0.7)
plotEigengeneNetworks(ME.genes, "", marDendro = c(0,4,1,2), marHeatmap = c(3,4,1,2), cex.adjacency = 0.3, cex.preservation = 0.3, plotPreservation = "standard")
dev.off()

write_csv(data.frame(table(module.colors)), path = "./final_eigengenes.csv") 
gene.info$kscaled <- by(gene.info, gene.info$module.colors, select, kWithin) %>% llply(function(x) { x / max (x) }) %>% reduce(c)
save(gene.info, file = "./save/gene.info.rda")

gene.module.membership <- as.data.frame(cor(intensities.PEER, ME.genes, use = "p"))
module.membership.pvalue <- as.data.frame(corPvalueStudent(as.matrix(gene.module.membership), nrow(intensities.PEER)))
names(gene.module.membership) <- str_replace(names(ME.genes),"ME", "MM.")
colnames(module.membership.pvalue) <- str_replace(names(ME.genes),"ME", "MM.pvalue.")

module.membership <- cbind(select(gene.info, ProbeName:Description), gene.module.membership, module.membership.pvalue)
write_csv(module.membership, "module_membership.csv")

colnames(gene.module.membership) %<>% str_replace("MM.", "")
colnames(module.membership.pvalue) %<>% str_replace("MM.pvalue.", "")
gene.module.membership$ProbeName <- rownames(gene.module.membership)
module.membership.pvalue$ProbeName <- rownames(module.membership.pvalue)

gene.module.membership.long <- gather(gene.module.membership, module.comparison, correlation, ivory:grey)
module.membership.pvalue.long <- gather(module.membership.pvalue, module.comparison, p.value, ivory:grey)
membership.join <- join(gene.module.membership.long, module.membership.pvalue.long)
eigengene.connectivity <- join(membership.join, gene.info) %>% select(ProbeName, GeneName:kscaled, module.comparison:p.value)
write_csv(eigengene.connectivity, "eigengene_connectivity.csv")

all.smooth <- apply(ME.genes, 2, smooth.spline, spar = 0.4) %>% llply(`[`, "y") 
smooth.df <- data.frame(all.smooth)
colnames(smooth.df) <- names(all.smooth)
smooth.df$x <- as.factor(1:nrow(smooth.df))
smooth.plot <- melt(smooth.df, id.vars = "x")

gen.pcaplot("all_principal_components", smooth.plot, FALSE, 10, 15)
gen.pcaplot("facet_principal_components", smooth.plot, TRUE, 13, 25)

sample.ids <- factor(rownames(intensities.PEER), levels = rownames(intensities.PEER))
colnames(ME.genes) %<>% str_replace("ME", "")
ME.genes.plot <- mutate(data.frame(ME.genes), Sample.ID = sample.ids)
expr.data.plot <- data.frame(t(intensities.PEER), module.colors)
by(expr.data.plot, expr.data.plot$module.colors, gen.heatmap, ME.genes.plot)

annot.reduce <- select(intensities.avg.norm$genes, ProbeName, SystematicName, GeneName, Description)
modules.out <- select(expr.data.plot, module.colors)
modules.out$ProbeName <- rownames(modules.out)
modules.out %<>% join(annot.reduce)
write.xlsx(modules.out, "modules_out.xlsx")

targets.age <- select(targets.gaa, SampleID, Age) %>% filter(!is.na(Age))
targets.age$Age %<>% as.numeric
cor.age <- gen.cor(ME.genes, targets.age)

targets.sex <- filter(targets.gaa, Sex != "UNKNOWN")
targets.sex.m <- model.matrix( ~ 0 + factor(targets.sex$Sex) )[,-1] %>% data.frame 
colnames(targets.sex.m) <- "Sex"
targets.sex.m %<>% mutate(SampleID = targets.gaa$SampleID)
cor.sex <- gen.cor(ME.genes, targets.sex.m)

targets.gaa.cor <- select(targets.gaa, SampleID, GAA1) %>% filter(!is.na(GAA1)) 
cor.gaa <- gen.cor(ME.genes, targets.gaa.cor)

targets.onset <- select(targets.gaa, SampleID, Onset) %>% filter(!is.na(Onset))
cor.onset <- gen.cor(ME.genes, targets.onset)

module.traits.all <- cbind(cor.age, cor.sex, cor.gaa, cor.onset) %>% data.frame
module.traits.pval <- select(module.traits.all, contains("p.value")) %>% as.matrix
module.traits.cor <- select(module.traits.all, -contains("p.value")) %>% as.matrix

module.trait.out <- data.frame(Module = rownames(module.traits.cor), module.traits.cor, module.traits.pval)
write_csv(module.trait.out, "module_trait_cor.csv")
text.matrix.traits <- paste(signif(module.traits.cor, 2), '\n(', signif(module.traits.pval, 1), ')', sep = '')
dim(text.matrix.traits) = dim(module.traits.cor)
gen.text.heatmap(module.traits.cor, text.matrix.traits, colnames(module.traits.cor), colnames(ME.genes), "", "module-trait relationships")

#text.matrix.PEER <- paste(signif(module.PEER.cor, 2), '\n(', signif(module.PEER.pval, 1), ')', sep = '')
#text.matrix.PEERcor <- paste(signif(cor.PEER.cor, 2), '\n(', signif(cor.PEER.pval, 1), ')', sep = '')

#dim(text.matrix.PEER) = dim(module.PEER.cor)
#dim(text.matrix.PEERcor) = dim(cor.PEER.cor)

#write_csv(module.PEER.out, "module_PEER_cor.csv")
#write_csv(cor.PEER.out, "cor_PEER_cor.csv")

#module.PEER.out <- data.frame(Factor = rownames(module.PEER.cor), module.PEER.cor, module.PEER.pval)
#cor.PEER.out <- data.frame(Module = rownames(cor.PEER), cor.PEER.cor, cor.PEER.pval)

#module.PEER.all <- cor.gaa.PEER %>% data.frame
#module.PEER.pval <- select(module.PEER.all, contains("p.value")) %>% as.matrix
#module.PEER.cor <- select(module.PEER.all, -contains("p.value")) %>% as.matrix

#cor.PEER.cor <- select(cor.PEER, -contains("p.value")) %>% as.matrix
#cor.PEER.pval <- select(cor.PEER, contains("p.value")) %>% as.matrix

#gen.text.heatmap(module.PEER.cor, text.matrix.PEER, colnames(module.PEER.cor), rownames(module.PEER.cor), "", "PEER factor-trait relationships")
#gen.text.heatmap(cor.PEER.cor, text.matrix.PEERcor, colnames(cor.PEER.cor), colnames(ME.genes), "", "module-PEER factor relationships")
#module.traits.all <- cbind(cor.status, cor.age, cor.sex) %>% data.frame
#module.traits.pval <- select(module.traits.all, contains("p.value")) %>% as.matrix
#module.traits.cor <- select(module.traits.all, -contains("p.value")) %>% as.matrix

#module.PEER.all <- cbind(cor.gaa.PEER, cor.onset.PEER) %>% data.frame
#module.PEER.pval <- select(module.PEER.all, contains("p.value")) %>% as.matrix
#module.PEER.cor <- select(module.PEER.all, -contains("p.value")) %>% as.matrix

#cor.PEER.cor <- select(cor.PEER, -contains("p.value")) %>% as.matrix
#cor.PEER.pval <- select(cor.PEER, contains("p.value")) %>% as.matrix

#module.trait.out <- data.frame(Module = rownames(module.traits.cor), module.traits.cor, module.traits.pval)
#module.PEER.out <- data.frame(Factor = rownames(module.PEER.cor), module.PEER.cor, module.PEER.pval)
#cor.PEER.out <- data.frame(Module = rownames(cor.PEER), cor.PEER.cor, cor.PEER.pval)

#write_csv(module.trait.out, "module_trait_cor.csv")
#write_csv(module.PEER.out, "module_PEER_cor.csv")
#write_csv(cor.PEER.out, "cor_PEER_cor.csv")

#text.matrix.traits <- paste(signif(module.traits.cor, 2), '\n(', signif(module.traits.pval, 1), ')', sep = '')
#text.matrix.PEER <- paste(signif(module.PEER.cor, 2), '\n(', signif(module.PEER.pval, 1), ')', sep = '')
#text.matrix.PEERcor <- paste(signif(cor.PEER.cor, 2), '\n(', signif(cor.PEER.pval, 1), ')', sep = '')
#dim(text.matrix.traits) = dim(module.traits.cor)
#dim(text.matrix.PEER) = dim(module.PEER.cor)
#dim(text.matrix.PEERcor) = dim(cor.PEER.cor)

#gen.text.heatmap(module.traits.cor, text.matrix.traits, colnames(module.traits.cor), colnames(ME.genes), "", "module-trait relationships")
#gen.text.heatmap(module.PEER.cor, text.matrix.PEER, colnames(module.PEER.cor), rownames(module.PEER.cor), "", "PEER factor-trait relationships")
#gen.text.heatmap(cor.PEER.cor, text.matrix.PEERcor, colnames(cor.PEER.cor), colnames(ME.genes), "", "module-PEER factor relationships")

#targets.final$SampleID %<>% str_replace(" ", "")
#rownames(ME.genes) <- rownames(intensities.PEER)

#targets.age <- select(targets.final, SampleID, Age) %>% filter(!is.na(Age))
#targets.age$Age %<>% as.numeric
#cor.age <- gen.cor(ME.genes, targets.age)
##cor.age.PEER <- gen.cor(PEER.factors, targets.age)

#targets.sex <- filter(targets.final, Sex != "UNKNOWN")
#targets.sex.m <- model.matrix( ~ 0 + factor(targets.sex$Sex) )[,-1] %>% data.frame #%>% mutate(SampleID = targets.sex$SampleID)
#colnames(targets.sex.m) <- c("Sex")
#targets.sex.m %<>% mutate(SampleID = targets.final$SampleID)
#cor.sex <- gen.cor(ME.genes, targets.sex.m)
##cor.sex.PEER <- gen.cor(PEER.factors, targets.sex.m)

#targets.gaa <- select(targets.final, SampleID, GAA1 ) %>% filter(!is.na(GAA1))# %>% filter(!is.na(GAA2))
##cor.gaa <- gen.cor(ME.genes, targets.gaa)
#cor.gaa.PEER <- gen.cor(PEER.factors, targets.gaa)

#targets.onset <- select(targets.final, SampleID, Onset) %>% filter(!is.na(Onset))
#cor.onset.PEER <- gen.cor(PEER.factors, targets.onset)

#targets.diff <- select(targets.final, SampleID, Diff.Baseline)
#cor.diff <- gen.cor(ME.genes, targets.diff)

#PEER.factors.df <- mutate(PEER.factors, SampleID = targets.final$SampleID)
#cor.PEER <- gen.cor(ME.genes, PEER.factors.df) %>% data.frame
