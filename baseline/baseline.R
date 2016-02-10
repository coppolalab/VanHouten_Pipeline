#Functional programming
library(magrittr)
library(purrr)
library(functional)

#String operations
library(stringr)

#Plotting
library(ggplot2)
library(extrafont)
library(Cairo)
library(heatmap.plus)
library(gplots) #for heatmap.2
library(RColorBrewer)

#Reading and writing tables
library(readr)
library(openxlsx)
library(readODS)

#For DE analysis
library(Biobase)
library(marray)
library(limma)
library(MASS)
library(matrixStats)
library(WGCNA)

#For batch correction and PEER
library(sva)
library(peer)

#Data arrangement
library(reshape2)
library(plyr)
library(dplyr)

#Boxplot function
gen.boxplot <- function(filename, dataset, targetset, maintext, ylabtext)
{
    dataset %<>% t %>% data.frame
    dataset.addvars <- mutate(dataset, Sample.Status = rownames(dataset), Status = targetset$Status)
    dataset.m <- melt(dataset.addvars, id = c("Sample.Status", "Status"))
    p <- ggplot(dataset.m, aes(x = Sample.Status, y = value, fill = factor(Status))) + geom_boxplot() + theme_bw()
    #p <- p + scale_fill_manual(values = colorscheme)
    p <- p + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1.0, size = 5))     
    p <- p + ggtitle(maintext) + ylab(ylabtext) + xlab("Sample") + theme(axis.text.x = element_text(size = 3))
    p <- p + theme(legend.position = "none", panel.grid.major = element_blank(), panel.grid.minor = element_blank())
    ggsave(filename = filename, plot = p, family = "Oxygen", width = 20 , height = 8)
}

#Heatmap bar function
gen.heatmapbars <- function(diagnosis.colorscheme, targetset)
{
    diagnosis.heatmap <- data.frame("Status" = levels(targetset$Status), "Diagnosis.Color" = diagnosis.colors)
    colorscheme <- data.frame("Status" = targetset$Status) %>% join(diagnosis.heatmap)
    colorscheme <- as.matrix(subset(colorscheme, select = c("Diagnosis.Color")))
    return(colorscheme)
}

#Heatmap function
gen.heatmap <- function(filename, dataset, targetset, heatmap.bars, maintitle)
{
    intensities1.cor <- cor(dataset)
    CairoPDF(filename, family = "Oxygen-Regular", width = 10, height = 10)
    heatmap.plus(intensities1.cor, col = heat.colors(40), scale = "none", cexCol = 0.4, cexRow = 0.4, main = maintitle)
    dev.off()
}

#IAC detection of outliers vix 
gen.IACcluster <- function(filename, dataset, maintitle)
{
    IAC = cor(dataset, use = "p")
    cluster1 = hclust(as.dist(1 - IAC))
    CairoPDF(filename, family = "Oxygen-Sans", width = 13, height = 10)
    plot(cluster1, main = paste(maintitle, " (no = ", dim(IAC)[2], ")"))
    dev.off()
    return(IAC)
}

#Create plot of standard deviations of all interarray correlations.  
gen.sdplot <- function(filename, dataset, maintitle)
{
    meanIAC <- apply(dataset, 2, mean)
    sdCorr <- sd(meanIAC)
    numbersd <- (meanIAC - mean(meanIAC)) / sdCorr
    numbersd.plot <- data.frame("Sample.Status" = names(numbersd), "Sample.Num" = seq(1:length(numbersd)), "Z.score" = numbersd)

    p <- ggplot(numbersd.plot, aes(x = Sample.Num, y = Z.score, label = Sample.Status) )
    p <- p + geom_text(size = 4, colour = "red")
    p <- p + geom_hline(aes(yintercept = -2)) + geom_hline(yintercept = -3) 
    p <- p + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
    p <- p + xlab("Sample Number") + ylab("Z score") + ggtitle(maintitle)
    CairoPDF(filename, family = "Oxygen-Regular", width = 10, height = 10)
    print(p)
    dev.off()
    return(numbersd)
}

#Median absolute deviation standardization function
standardize <- function(dataset)
{
    rowmed <- apply(dataset, 1, median)
    rowmad <- apply(dataset, 1, mad)
    rv <- sweep(dataset, 1, rowmed)
    rv <- sweep(rv, 1, rowmad, "/")
    return(rv)
}

#Create heatmap of top genes
gen.topgenes <- function(filename, dataset, heatmap.bars, maintitle, rowmads, num.genes)
{
    top.genes.names <- rowmads[1:num.genes]
    top.genes.intensities <- dataset[top.genes.names,]
    top.genes.dist <- dist(t(standardize(top.genes.intensities)))
    top.genes.clust <- hclust(top.genes.dist)
    top.genes.matrix <- as.matrix(top.genes.dist)

    CairoPDF(filename, family = "Oxygen-Regular", width = 10, height = 10)
    heatmap.plus(top.genes.matrix, col = rev(heat.colors(75)), distfun = function (x) as.dist(x), main = maintitle, scale = "none", cexRow = 0.08, cexCol = 0.08)
    dev.off()
    return(top.genes.dist)
}

#MDS function - may need work
gen.pca <- function(filename, dataset, targetset, colorscheme, variablename)
{
    dataset.plot <- data.frame(rownames(dataset$points), dataset$points)
    target.data <- data.frame(targetset$SampleID, factor(targetset[,variablename]))
    colnames(target.data) <- c("Sample.Status", variablename)
    colnames(dataset.plot) <- c("Sample.Status", "Component.1", "Component.2")
    dataset.plot <- merge(dataset.plot, target.data)
    p <- ggplot(dataset.plot, aes_string(x = "Component.1", y = "Component.2", col = variablename)) + geom_point() 
    p <- p + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
    p <- p + scale_fill_manual(values = colorscheme) + xlim(-80,80)
    p <- p + xlab("Component 1") + ylab("Component 2") + ggtitle(variablename)
    CairoPDF(file = filename, height = 7, width = 7)
    print(p)
    dev.off()
}

#Run statistical cutoff tests
gen.decide <- function(test, fit.object, write.results)
{
    results <- decideTests(fit.object, adjust.method = test[1], p = as.numeric(test[2])) #Run test at specified cutoff
    if(write.results == TRUE)
    {
        write.fit(file = paste("./fit_", test[1], ".tsv", sep = ""), fit.object, adjust = test[1], results = results)
    }
    num.genes <- length(which(apply(results, 1, function (x) any(x, na.rm = T))))  #Make this better
    mysum <- summary(results)[-2,] #Eliminate the row for no change in expression
    #mysum[2,] <- -(mysum[2,])
    #names(mysum) <- "Patient_vs._Control"
    mysum <- data.frame("Test" = paste(test[1], " p<", test[2], sep = ""), "Num" = paste(num.genes, "Genes", sep = " "), "Direction" = c("positive", "negative"), Patient_vs._Control = mysum)
    return(mysum)
}

#Plot statistical cutoff tests
gen.decideplot <- function(filename, decide.plot)
{
    decide.plot$variable <- str_replace_all(decide.plot$variable, "_", " ")
    p <- ggplot()
    p <- p + geom_bar(data = subset(decide.plot, Direction == "positive"),  aes(x = variable, y = value), stat = "identity", colour = "black", fill = "red", position = "dodge")   
    p <- p + geom_text(data = subset(decide.plot, Direction == "positive"), stat = "identity", size = 4, aes(x = variable, y = value, ymax = max(value) + 110, label = value), hjust = -0.3, position = position_dodge(width = 1))
    p <- p + geom_bar(data = subset(decide.plot, Direction == "negative"),  aes(x = variable, y = -value), stat = "identity", colour = "black", fill = "green", position = "dodge") 
    p <- p + geom_text(data = subset(decide.plot, Direction == "negative"), stat = "identity", size = 4, aes(x = variable, y = -value, ymax = min(-value) - 110, label = abs(value)), hjust = 1.3, position = position_dodge(width = 1))
    if (length(unique(decide.plot$Test)) > 1)
    {
        p <- p + facet_grid(Num + Test ~ .) 
        #p <- p + ggtitle("Threshold Selection")
    }
    else
    {
        p <- p + ggtitle(paste(decide.plot$Test, "\n", decide.plot$Num))
    }
    p <- p + theme_bw() + coord_flip() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) 
    p <- p + theme(axis.title.y = element_blank(), axis.ticks.y = element_blank(), axis.text.y = element_text(hjust = 0)) + ylab("Differentially Expressed Genes")
    CairoPDF(filename, family = "Oxygen", width = 8.2, height = 7)
    print(p)
    dev.off()
}

gen.pval.hist <- function(filename, fit.pvals)
{
    colnames(fit.pvals) <- c("Patient_vs._Control")
    fit.pvals.plot <- melt(fit.pvals)
    fit.pvals.plot$Contrasts <- str_replace_all(fit.pvals.plot$Contrasts, "_", " ")
    p <- ggplot(fit.pvals.plot, aes(x = value)) + geom_histogram(binwidth = 1/80) + facet_grid(. ~ Contrasts)
    p <- p + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
    p <- p + ggtitle("P-value distribution across contrasts") + theme(axis.title.y = element_blank()) + xlab("p-value")
    CairoPDF(filename, height = 7, width = 21)
    print(p)
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

#Calculate ratios.  Really needs work!
gen.ratios <- function(dataset, targetset)
{
    all.samples <- data.frame(dataset)
    all.controls <- all.samples[,targetset$Status == "Control"]
    all.controls.means <- rowMeans(all.controls)
    all.patients <- all.samples[,targetset$Status == "Patient"]

    all.coefficients.pat.cont <- all.patients - all.controls.means

    all.coefficients <- data.frame("Probe" = rownames(all.coefficients.pat.cont), all.coefficients.pat.cont)
    all.samples <- data.frame("Probe" = rownames(all.samples), all.samples)
    colnames(all.samples)[2:length(all.samples)] <- paste(colnames(all.samples[2:length(all.samples)]), "expr", sep = ".") 
    ratio.exp <- merge(all.coefficients, all.samples)
    return(ratio.exp)
}

#Generate fit object
gen.fit <- function(dataset, model.design)
{
    fit <- lmFit(dataset, model.design)
    contrasts.anova <- makeContrasts(time1.patient_vs_time1.control = Patient - Control, levels = model.design)
    fit2.anova <- contrasts.fit(fit, contrasts.anova)
    fitb <- eBayes(fit2.anova)
}

gen.workbook <- function(dataset, filename)
{
    pval.cols <- colnames(dataset) %>% str_detect("p.value") %>% which
    coef.cols <- colnames(dataset) %>% str_detect("Coef.") %>% which
    colnames(dataset)[coef.cols] %<>% str_replace("Coef.", "")

    wb <- createWorkbook()
    addWorksheet(wb = wb, sheetName = "Sheet 1", gridLines = TRUE)
    writeDataTable(wb = wb, sheet = 1, x = dataset, withFilter = FALSE)
    sig.pvalues <- createStyle(fontColour = "red")
    conditionalFormatting(wb, 1, cols = pval.cols, rows = 1:nrow(dataset), rule = "<0.05", style = sig.pvalues)
    conditionalFormatting(wb, 1, cols = coef.cols, rows = 1:nrow(dataset), style = c("#63BE7B", "white", "red"), type = "colourScale")
    setColWidths(wb, 1, cols = 1:3, widths = "auto")
    #setColWidths(wb, 1, cols = 2, widths = 15)
    setColWidths(wb, 1, cols = 4, widths = 45)
    setColWidths(wb, 1, cols = 5:7, widths = 15)
    pageSetup(wb, 1, orientation = "landscape", fitToWidth = TRUE)
    freezePane(wb, 1, firstRow = TRUE)
    showGridLines(wb, 1, showGridLines = TRUE)
    modifyBaseFont(wb, fontSize = 10.5, fontName = "Oxygen")

    saveWorkbook(wb, filename, overwrite = TRUE) 
}

#Create genelists
gen.tables <- function(dataset, results, annot, ratio.exp, suffix)
{
    treat.de <- data.frame("ProbeName" = rownames(results), dataset)
    anovalist <- apply(results, 1, any, na.rm = T) %>% which
    treat.de.anova <- treat.de[anovalist,]
    annot %<>% select(-Row, -Col)
    fitsel.ratio <- merge(treat.de.anova, annot) %>% merge(ratio.exp, by.x = "ProbeName", by.y = "Probe")
    fitsel.ratio <- filter(fitsel.ratio, !duplicated(fitsel.ratio))
    fitsel.ratio2 <- select(fitsel.ratio, ProbeName, SystematicName, GeneName, Description, Coef, contains("p.value"), contains("time1"), t, A, matches("GSM")) 

    fitsel.ratio.all <- merge(treat.de, annot) %>% merge(ratio.exp, by.x = "ProbeName", by.y = "Probe")
    fitsel.ratio.all <- filter(fitsel.ratio.all, !duplicated(fitsel.ratio.all))
    fitsel.ratio.all <- select(fitsel.ratio.all, ProbeName, SystematicName, GeneName, Coef, contains("p.value"), contains("time1"), t, A, Description, matches("GSM")) %>% arrange(desc(t))

    coef.cols <- colnames(fitsel.ratio2) %>% str_detect("Coef") %>% which
    colnames(fitsel.ratio2)[coef.cols] <- "Coef.Patient vs. Control"

    pco.abs <- as.numeric(fitsel.ratio2$"Coef.Patient vs. Control") %>% abs
    fitsel.ratio2 %<>% mutate(pco.abs) %>% arrange(desc(pco.abs)) %>% select(-pco.abs)
    gen.workbook(fitsel.ratio2, paste("./significant_geneList_", suffix, "_time1.xlsx", sep = ""))
    fitsel.return <- fitsel.ratio2

    ##Patient - Control
    #fitsel.pco <- select(fitsel.ratio2, Accession, Symbol, Definition, matches("Coef.Patient vs. Control")) %>% mutate(pco.abs) %>% arrange(desc(pco.abs)) %>% select(-pco.abs)
    #gen.small.workbook(fitsel.pco, paste("./significant_geneList_", suffix, "_time1_pco.xlsx", sep = ""))

    #write_csv(fitsel.ratio.all, path = paste("./complete_genelist_time1_", suffix, ".csv", sep = ""))
    return(fitsel.return)
}

gen.small.workbook <- function(dataset, filename)
{
    coef.cols <- colnames(dataset) %>% str_detect("Coef") %>% which
    pval.cols <- colnames(dataset) %>% str_detect("p.value") %>% which
    #colnames(dataset)[coef.cols] %<>% str_replace("Coef.", "")

    wb <- createWorkbook()
    addWorksheet(wb = wb, sheetName = "Sheet 1", gridLines = TRUE)
    writeDataTable(wb = wb, sheet = 1, x = dataset, withFilter = FALSE)
    sig.pvalues <- createStyle(fontColour = "red")
    conditionalFormatting(wb, 1, cols = coef.cols, rows = 1:nrow(dataset), style = c("#63BE7B", "white", "red"), type = "colourScale")
    conditionalFormatting(wb, 1, cols = pval.cols, rows = 1:nrow(dataset), rule = "<0.05", style = sig.pvalues)
    #setColWidths(wb, 1, cols = 2, widths = "auto")
    #setColWidths(wb, 1, cols = 1, widths = "auto")
    #setColWidths(wb, 1, cols = 3, widths = 45)
    #setColWidths(wb, 1, cols = 4:6, widths = 15)
    pageSetup(wb, 1, orientation = "landscape", fitToWidth = TRUE)
    freezePane(wb, 1, firstRow = TRUE)
    showGridLines(wb, 1, showGridLines = TRUE)
    modifyBaseFont(wb, fontSize = 10.5, fontName = "Oxygen")

    saveWorkbook(wb, filename, overwrite = TRUE) 
}
    
#Make anova objects
gen.anova <- function(dataset, suffix)
{
    plot.pc <- filter(dataset, time1.patient_vs_time1.control != 0) %>% select(-matches("expr")) %>% select(matches("GSM"))
    gen.anova.heatmap(paste("./9_anova_heatmap_patient_vs_control", suffix, sep = "_"), plot.pc, "Patients vs. Controls")
}

#Generate anova heatmaps
gen.anova.heatmap <- function(filename, dataset, maintitle)
{ CairoPDF(filename, width = 10, height = 10)
    heatmap.2(as.matrix(dataset), col = rev(redgreen(48)), breaks=(c(-3, -2.5, -2, -1.5, seq(-1, 1, 0.05), 1.5, 2, 2.5, 3)), trace = "none", cexCol = 0.7, labRow = "", keysize = 0.9)
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
    width.dynamic <- 3 + (1 * ncol(text.matrix))
    height.dynamic <- 3 + (1 * nrow(text.matrix))
    cex.dynamic <- 0.5 + ((30 - nrow(text.matrix)))/30
    CairoPDF(filename, width = width.dynamic, height = height.dynamic)
    par(mar = c(8, 8, 3, 3))
    labeledHeatmap(Matrix = cor.dataset, xLabels = x.names, yLabels = y.names, ySymbols = y.names, yColorLabels = TRUE, colors = greenWhiteRed(50), textMatrix = text.matrix, setStdMargins = F, cex.text = cex.dynamic, cex.lab = cex.dynamic, zlim = c(-1,1), main = maintitle)
    dev.off()
}


targets <- read.ods("../Van Houten pheno.ods") %>% data.frame
targets[targets == "NA"] <- NA
colnames(targets) <- targets[1,]
targets <- targets[-1,]
targets$Status %<>% factor
targets$Sex %<>% factor
targets$Line %<>% factor

targets.final <- filter(targets, Line == 0) %>% filter(!is.na(Sex))
targets.final$SampleID <- paste(targets.final$Prefix, targets.final$ID, sep = "")
save(targets.final, file = "./save/targets.final.rda")

intensities <- read.maimages(source = "agilent", path = "../GSE11204_RAW/RAW", ext="txt")
intensities.bg <- backgroundCorrect(intensities, method = "normexp", offset = 50)
intensities.wa <- normalizeWithinArrays(intensities.bg, method = "loess")
intensities.ba <- normalizeBetweenArrays(intensities.wa, method = "quantile")
intensities.avg <- avereps(intensities.wa, intensities.wa$genes$ProbeName)
intensities.avg.norm <- avereps(intensities.ba, intensities.ba$genes$ProbeName)
save(intensities.avg.norm, file = "./save/intensities.avg.norm.rda")

intensities.mat <- intensities.avg$A[-(1:2),]
intensities.mat.norm <- intensities.avg.norm$A[-(1:2),]
intensities1 <- select(data.frame(intensities.mat), one_of(targets.final$SampleID)) %>% as.matrix #%>% normalizeBetweenArrays(method = "quantile")
intensities1.norm <- select(data.frame(intensities.mat.norm), one_of(targets.final$SampleID)) %>% as.matrix %>% normalizeBetweenArrays(method = "quantile")
save(intensities1.norm, file = "./save/intensities1.norm.rda")

gen.boxplot("1_base_signal.jpg", intensities1, targets.final, "Signal intensity not normalized", "Intensity")
gen.boxplot("1_base_signal_norm.jpg", intensities1.norm, targets.final, "Signal intensity normalized", "Intensity")

#Make heatmap of all gene intensities
diagnosis.colors <- c("magenta", "green")
heatmap.bars1 <- gen.heatmapbars(diagnosis.colors, targets.final) 
gen.heatmap("2_base_heatmap", intensities1.norm, targets.final, heatmap.bars1, "Clustering Based on Inter-Array Pearson Coefficient")

#IAC detection of outliers
IAC <- gen.IACcluster("2_base_IAC.pdf", intensities1.norm, "All samples")
gen.sdplot("2_base_IAC_sd", IAC, "label-1 samples")

#Top 500 and 1000 genes
intensities1.mads <- apply(intensities1.norm, 1, mad)
intensities1.ordered <- order(intensities1.mads, decreasing = TRUE)

top500.dist <- gen.topgenes("4_base_heatmap_mostVariable_500", intensities1.norm, heatmap.bars1, "Clustering Based on the Top 500 Most Variable Genes", intensities1.ordered, 500)
top1000.dist <- gen.topgenes("4_base_heatmap_mostVariable_1000", intensities1.norm, heatmap.bars1, "Clustering Based on the Top 1000 Most Variable Genes", intensities1.ordered, 1000)

#Principal components analysis
cm1 <- cmdscale(top1000.dist, eig = TRUE)
gen.pca("4_base_MDS_Status", cm1, targets.final, diagnosis.colors, "Status")

targets.final$Age %<>% as.numeric
pc.points <- data.frame(cm1$points, targets.final$SampleID, targets.final$Status, targets.final$Sex, targets.final$Age)
colnames(pc.points) <- c("Component.1", "Component.2", "SampleID", "Status", "Sex", "Age")
pc.patients <- filter(pc.points, Status == "Patient")
min.component1 <- min(pc.patients$Component.1)
max.component2 <- max(pc.patients$Component.2)

new.points <- filter(pc.points, Component.1 > min.component1 & Component.2 < max.component2)
new.points.key <- paste(new.points$SampleID, collapse = "|")
old.points <- filter(pc.points, Component.1 < min.component1 | Component.2 > max.component2)
control.points <- filter(new.points, Status == "Control")

targets.final %<>% filter(grepl(new.points.key, SampleID))
intensities1.norm %<>% data.frame %>% select(one_of(targets.final$SampleID)) %>% as.matrix %>% normalizeBetweenArrays(method = "quantile")

model.status <- model.matrix( ~ 0 + factor(targets.final$Status) )
colnames(model.status) <- c("Control", "Patient")
model.status.reduce <- model.status[,-1]

model.sex <- model.matrix( ~ 0 + factor(targets.final$Sex) )
colnames(model.sex) <- c("Female", "Male")
model.sex.reduce <- model.sex[,-1]

model.PEER <- cbind(model.status.reduce, Male = model.sex.reduce, Age = log2(as.numeric(targets.final$Age)))
colnames(model.PEER)[1] <- "Patient"
save(model.PEER, file = "./save/model.PEER.rda")

gen.peer(4, intensities1.norm, TRUE, model.PEER)
model.PEER_covariate <- read_csv("./factor_4.csv") %>% select(-(X1:X4))
rownames(model.PEER_covariate) <- colnames(intensities1.norm)
colnames(model.PEER_covariate) <- paste("X", 1:ncol(model.PEER_covariate), sep = "")

PEER.weights <- read_csv("./weight_4.csv") %>% select(-(X1:X4))
PEER.weights.sums <- colSums(abs(PEER.weights)) %>% data.frame
PEER.weights.sums$Factor <- 1:nrow(PEER.weights.sums)
colnames(PEER.weights.sums)[1] <- "Weight"

p <- ggplot(PEER.weights.sums, aes(x = factor(Factor), y = as.numeric(Weight), group = 1)) + geom_line(color = "blue") 
p <- p + theme_bw() + xlab("Factor") + ylab("Weigtht")
CairoPDF("./PEER_weights", height = 4, width = 6)
print(p)
dev.off()

targets.final$GAA1 %<>% as.numeric
targets.final$GAA2 %<>% as.numeric 
targets1.gaa <- select(targets.final, SampleID, GAA1) %>% filter(!is.na(GAA1))
cor.gaa <- gen.cor(model.PEER_covariate, targets1.gaa)

targets.final$Onset %<>% as.numeric 
targets1.onset <- select(targets.final, SampleID, Onset) %>% filter(!is.na(Onset)) 
cor.onset <- gen.cor(model.PEER_covariate, targets1.onset)

PEER.traits.all <- cbind(cor.gaa, cor.onset) %>% data.frame
PEER.traits.pval <- select(PEER.traits.all, contains("p.value")) %>% as.matrix
PEER.traits.cor <- select(PEER.traits.all, -contains("p.value")) %>% as.matrix

text.matrix.PEER <- paste(signif(PEER.traits.cor, 2), '\n(', signif(PEER.traits.pval, 1), ')', sep = '')
dim(text.matrix.PEER) <- dim(PEER.traits.cor)
gen.text.heatmap(PEER.traits.cor, text.matrix.PEER, colnames(PEER.traits.cor), rownames(PEER.traits.cor), "", "PEER factor-trait relationships")

PEER.trait.out <- data.frame(Factor = rownames(PEER.traits.cor), PEER.traits.cor, PEER.traits.pval)
write_csv(PEER.trait.out, "PEER_trait_cor.csv")

#Calculate ratios for use in tables
ratio.exp <- gen.ratios(intensities1.norm, targets.final)
save(ratio.exp, file = "./save/ratio.exp.rda")

#Linear model fitting
model.cov <- cbind(model.status, Male = model.sex.reduce, Age = log2(as.numeric(targets.final$Age)))
model.full <- cbind(model.cov, model.PEER_covariate)

fit.object <- gen.fit(intensities1.norm, model.cov)
save(fit.object, file = "./save/fit.object.rda")

#Generate statisical cutoff
decide <- list(c("fdr", 0.01), c("fdr", 0.05), c("fdr", 0.1), c("none", 0.001), c("none", 0.005), c("none", 0.01))
decide.plot <- ldply(decide, gen.decide, fit.object, FALSE) %>% melt(id.vars = c("Test", "Num", "Direction"))
gen.decideplot("./7_threshold_selection", decide.plot)

results <- decideTests(fit.object, adjust.method = "none", p = 0.005)
results.fdr <- decideTests(fit.object, adjust.method = "fdr", p = 0.01)

#Make tables
de.object <- read_tsv("./fit_none.tsv")
de.object.fdr <- read_tsv("./fit_fdr.tsv")
fit.selection <- gen.tables(de.object, results, intensities.ba$genes, ratio.exp, "pLess005")
fit.selection.fdr <- gen.tables(de.object.fdr, results.fdr, intensities.ba$genes, ratio.exp, "fdrLess001")

#P-value histogram
gen.pval.hist("./8_hist_pvalue", fit.object$p.value)

#Anova heatmaps
gen.anova(fit.selection, "none")
gen.anova(fit.selection.fdr, "fdr")

ourgenes <- read_csv("../../FRDA project/baseline/complete_genelist_time1_fdrLess01.csv")
colnames(ourgenes) %<>% str_replace_all("\\.vs\\.\\.", "_vs_")
#ourgenes %<>% select(-Res.time1.patient_vs_time1.control)
ourgenes.reduce <- select(ourgenes, Accession:Definition, Coef.time1.patient_vs_time1.control, p.value.time1.patient_vs_time1.control, p.value.adj.time1.patient_vs_time1.control, Res.time1.patient_vs_time1.control) #%>% filter(Res.time1.patient_vs_time1.control != 0)
ourgenes.join <- select(ourgenes.reduce, Symbol, Res.time1.patient_vs_time1.control)
save(ourgenes.reduce, file = "./save/ourgenes.reduce.rda")

complete.genes <- read_csv("./complete_genelist_time1_fdrLess001.csv")
theirgenes.join <- select(complete.genes, GeneName, time1.patient_vs_time1.control)
test.join <- merge(theirgenes.join, ourgenes.join, by.x = "GeneName", by.y = "Symbol")

#test.join %<>% select(-GeneName)
colnames(test.join)[-1] <- c("VanHouten", "Coppola")
#sumV <- summary(test.join)
#v <- c("VanHouten", "Coppola") 
CairoPDF("overlap", width = 6, height = 6)
vennDiagram(results, names = v, main = "", include = c("up", "down"), counts.col=c(2,3), cex = 0.8)
dev.off()

fitsel.reduce <- select(complete.genes, GeneName, matches("Coef"), p.value, p.value.adj)
combined.genes <- merge(ourgenes, fitsel.reduce, by.x = "Symbol", by.y = "GeneName", all = FALSE)
colnames(combined.genes)[4:9] <- c("Coppola.Coef", "Coppola.p.value", "Coppola.p.value.adj", "VanHouten.Coef", "VanHouten.p.value", "VanHouten.p.value.adj")
combined.genes %<>% select(Accession, Symbol, Definition, Coppola.Coef, VanHouten.Coef, Coppola.p.value, Coppola.p.value.adj, VanHouten.p.value, VanHouten.p.value.adj)
gen.small.workbook(combined.genes, "combined.genes.xlsx")
