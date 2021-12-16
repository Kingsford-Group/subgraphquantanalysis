library(ggplot2)
library(cowplot)
library(assertthat)
library(DESeq2)
library(tximport)

options(stringsAsFactors = FALSE)
args <- commandArgs(trailingOnly = FALSE)

# Meta file for SRA ID
codedir_split <- unlist( strsplit( sub("--file=", "", args[grep("--file=", args)]), "/") )
if (length(codedir_split) > 2) {
        codedir <- paste0(codedir_split[1:(length(codedir_split)-2)], sep="/")
} else if (length(codedir_split) == 2) {
        codedir <- "./"
} else if (length(codedir_split) == 1) {
        codedir <- "../"
}
metafile <- paste0(codedir,"/data/Metadata_mcf10.txt",sep="")

# output directory
idx <- grep("--args", args)
folder <- paste0(args[idx+1], "/MCF10/", sep="")
print(paste0("Meta data file = ", metafile, sep=""))
print(paste0("Output folder = ", folder, sep=""))

# read Meta file
meta <- read.table(metafile, header=F, sep="\t")
meta <- meta[, c(5, 8, 10)]
colnames(meta) <- c("ID", "treatment", "time")

# read IValue
iv <- read.table(paste0(folder, "IValue_mean_ext.txt", sep=""), header=F, sep="\t")
colnames(iv) <- c("Name", "Group1Mean", "Group2Mean", "IV", "IV25")

# DE detection on transcript level
files <- paste0(folder, meta$ID, "/quant.sf", sep="")
names(files) <- meta$ID
txi <- tximport(files, type = "salmon", txOut = TRUE)
dds <- DESeqDataSetFromTximport(txi, meta, ~treatment)
dds <- DESeq(dds)
res <- results(dds, alpha=0.01)
print(sum(res$padj < 0.01, na.rm=TRUE))
res <- res[!is.na(res$padj),]
res$IV <- iv[match(rownames(res), iv$Name), "IV"]
res$IV25 <- iv[match(rownames(res), iv$Name), "IV25"]

# DE detection on gene level
gene_trans_map <- read.table(paste0(folder, "../gencode.v26.Gene_Trans_Map.txt", sep=""), header=F, sep="\t")
colnames(gene_trans_map) <- c("GENEID", "TXNAME")
gene_trans_map <- gene_trans_map[, c("TXNAME", "GENEID")]
txi <- tximport(files, type = "salmon", tx2gene = gene_trans_map)
dds_gene <- DESeqDataSetFromTximport(txi, meta, ~treatment)
dds_gene <- DESeq(dds_gene)
res_gene <- results(dds_gene, alpha=0.01)
res_gene <- res_gene[!is.na(res_gene$padj),]
print(sum(res_gene$padj < 0.01, na.rm=TRUE))

related_trans <- res[rownames(res) %in% gene_trans_map[gene_trans_map$GENEID %in% rownames(res_gene[res_gene$padj < 0.01, ]), "TXNAME"], ]
print(nrow(related_trans[related_trans$padj < 0.01, ]))
print(related_trans[related_trans$padj < 0.01, "IV25"])

# draw example curve
large_font <- 18

t <- read.table(paste0(folder, "IValue_curve_example.txt", sep=""), header=F, sep="\t")
colnames(t) <- c("Name", "Group", "Sample", "reference_proportion", "lb", "ub")
df <- data.frame(Name=rep(t$Name, 2), Group=rep(t$Group, 2), Sample=rep(t$Sample, 2), reference_proportion=rep(t$reference_proportion, 2), Expression=c(t$lb, t$ub), Type=c(rep("lower bound", nrow(t)), rep("upper bound", nrow(t))))

tname <- "ENST00000509980.5"
g <- "EGF Stimulation + DMSO"
tmpdf <- data.frame(Group=g,  reference_proportion=df[df$Name == tname & df$Sample == "Mean" & df$Type=="lower bound" & df$Group==g, "reference_proportion"], lb=df[df$Name == tname & df$Sample == "Mean" & df$Type=="lower bound" & df$Group==g, "Expression"], ub=df[df$Name == tname & df$Sample == "Mean" & df$Type=="upper bound" & df$Group==g, "Expression"])
g <- "no EGF + DMSO"
tmpdf <- rbind(tmpdf, data.frame(Group=g, reference_proportion=df[df$Name == tname & df$Sample == "Mean" & df$Type=="lower bound" & df$Group==g, "reference_proportion"], lb=df[df$Name == tname & df$Sample == "Mean" & df$Type=="lower bound" & df$Group==g, "Expression"], ub=df[df$Name == tname & df$Sample == "Mean" & df$Type=="upper bound" & df$Group==g, "Expression"]) )
p3 <- ggplot(df[df$Name == tname & df$Sample == "Mean", ]) +
	geom_ribbon_pattern(data=tmpdf, aes(x=1-reference_proportion, ymin=lb, ymax=ub, pattern=Group), alpha=0.5, fill='white', color="white") +
	geom_line(aes(x=1-reference_proportion, y=Expression, group=interaction(Group, Sample, Type), linetype=Type), size=1.2) + 
	geom_vline(xintercept=1-0.0423675453975355, linetype="dotted") + geom_vline(xintercept=1-0.0954330296051542, linetype="dotdash") + theme_bw() + labs(title = paste0("Ranges of optima of ", tname, sep=""), x = "reference completeness", y = "normalized abundance") + 
	theme(legend.title=element_blank()) + theme(axis.text.x = element_text(size = large_font), axis.title.x = element_text(size = large_font), axis.text.y = element_text(size = large_font), axis.title.y = element_text(size = large_font), legend.text=element_text(size=large_font))

tname <- "ENST00000308394.8"
g <- "EGF Stimulation + DMSO"
tmpdf <- data.frame(Group=g,  reference_proportion=df[df$Name == tname & df$Sample == "Mean" & df$Type=="lower bound" & df$Group==g, "reference_proportion"], lb=df[df$Name == tname & df$Sample == "Mean" & df$Type=="lower bound" & df$Group==g, "Expression"], ub=df[df$Name == tname & df$Sample == "Mean" & df$Type=="upper bound" & df$Group==g, "Expression"])
g <- "no EGF + DMSO"
tmpdf <- rbind(tmpdf, data.frame(Group=g, reference_proportion=df[df$Name == tname & df$Sample == "Mean" & df$Type=="lower bound" & df$Group==g, "reference_proportion"], lb=df[df$Name == tname & df$Sample == "Mean" & df$Type=="lower bound" & df$Group==g, "Expression"], ub=df[df$Name == tname & df$Sample == "Mean" & df$Type=="upper bound" & df$Group==g, "Expression"]) )
p4 <- ggplot(df[df$Name == tname & df$Sample == "Mean", ]) + 
	geom_ribbon_pattern(data=tmpdf, aes(x=1-reference_proportion, ymin=lb, ymax=ub, pattern=Group), alpha=0.5, fill='white', color="white") +
	geom_line(aes(x=1-reference_proportion, y=Expression, group=interaction(Group, Sample, Type), linetype=Type), size=1.2) + 
	geom_vline(xintercept=1-0.0356286092384825, linetype="dotted") + geom_vline(xintercept=1-0.0695868455633252, linetype="dotdash") + theme_bw() + labs(title = paste0("Ranges of optima of ", tname, sep=""), x = "reference completeness", y = "normalized abundance") + 
	theme(legend.title=element_blank()) + theme(axis.text.x = element_text(size = large_font), axis.title.x = element_text(size = large_font), axis.text.y = element_text(size = large_font), axis.title.y = element_text(size = large_font), legend.text=element_text(size=large_font))

tname <- "ENST00000443868.6"
g <- "EGF Stimulation + DMSO"
tmpdf <- data.frame(Group=g,  reference_proportion=df[df$Name == tname & df$Sample == "Mean" & df$Type=="lower bound" & df$Group==g, "reference_proportion"], lb=df[df$Name == tname & df$Sample == "Mean" & df$Type=="lower bound" & df$Group==g, "Expression"], ub=df[df$Name == tname & df$Sample == "Mean" & df$Type=="upper bound" & df$Group==g, "Expression"])
g <- "no EGF + DMSO"
tmpdf <- rbind(tmpdf, data.frame(Group=g, reference_proportion=df[df$Name == tname & df$Sample == "Mean" & df$Type=="lower bound" & df$Group==g, "reference_proportion"], lb=df[df$Name == tname & df$Sample == "Mean" & df$Type=="lower bound" & df$Group==g, "Expression"], ub=df[df$Name == tname & df$Sample == "Mean" & df$Type=="upper bound" & df$Group==g, "Expression"]) )
p5 <- ggplot(df[df$Name == tname & df$Sample == "Mean", ]) + 
	geom_ribbon_pattern(data=tmpdf, aes(x=1-reference_proportion, ymin=lb, ymax=ub, pattern=Group), alpha=0.5, fill='white', color="white") + 
	geom_line(aes(x=1-reference_proportion, y=Expression, group=interaction(Group, Sample, Type), linetype=Type), size=1.2) + 
	geom_vline(xintercept=1-0.0107085707561145, linetype="dotted") + geom_vline(xintercept=1-0.0180927714007716, linetype="dotdash") + theme_bw() + labs(title = paste0("Ranges of optima of ", tname, sep=""), x = "reference completeness", y = "normalized abundance") + 
	theme(legend.title=element_blank()) + theme(axis.text.x = element_text(size = large_font), axis.title.x = element_text(size = large_font), axis.text.y = element_text(size = large_font), axis.title.y = element_text(size = large_font), legend.text=element_text(size=large_font))

tname <- "ENST00000292433.3"
g <- "EGF Stimulation + DMSO"
tmpdf <- data.frame(Group=g,  reference_proportion=df[df$Name == tname & df$Sample == "Mean" & df$Type=="lower bound" & df$Group==g, "reference_proportion"], lb=df[df$Name == tname & df$Sample == "Mean" & df$Type=="lower bound" & df$Group==g, "Expression"], ub=df[df$Name == tname & df$Sample == "Mean" & df$Type=="upper bound" & df$Group==g, "Expression"])
g <- "no EGF + DMSO"
tmpdf <- rbind(tmpdf, data.frame(Group=g, reference_proportion=df[df$Name == tname & df$Sample == "Mean" & df$Type=="lower bound" & df$Group==g, "reference_proportion"], lb=df[df$Name == tname & df$Sample == "Mean" & df$Type=="lower bound" & df$Group==g, "Expression"], ub=df[df$Name == tname & df$Sample == "Mean" & df$Type=="upper bound" & df$Group==g, "Expression"]) )
p6 <- ggplot(df[df$Name == tname & df$Sample == "Mean", ]) +
	geom_ribbon_pattern(data=tmpdf, aes(x=1-reference_proportion, ymin=lb, ymax=ub, pattern=Group), alpha=0.5, fill='white', color="white") +
	geom_line(aes(x=1-reference_proportion, y=Expression, group=interaction(Group, Sample, Type), linetype=Type), size=1.2) +
	geom_vline(xintercept=1-0.618471, linetype="dotted") + geom_vline(xintercept=1-0.80263, linetype="dotdash") + theme_bw() + labs(title = paste0("Ranges of optima of ", tname, sep=""), x = "reference completeness", y = "normalized abundance") +
	theme(legend.title=element_blank()) + theme(axis.text.x = element_text(size = large_font), axis.title.x = element_text(size = large_font), axis.text.y = element_text(size = large_font), axis.title.y = element_text(size = large_font), legend.text=element_text(size=large_font))

p2 <- ggplot(as.data.frame(res[res$padj < 0.01, ])) + geom_histogram(aes(x=1-IV25)) + theme_cowplot() + labs(x = "reference completeness", title = "Number of DE transcripts that are unreliable \nunder the reference completeness parameter") + 
	theme(axis.text.x = element_text(size = large_font), axis.title.x = element_text(size = large_font), axis.text.y = element_text(size = large_font), axis.title.y = element_text(size = large_font), legend.text=element_text(size=large_font))

ptmp <- plot_grid(p3 + theme(legend.position="none"), p4 + theme(legend.position="none"), p5 + theme(legend.position="none"), p6 + theme(legend.position="none"), p2, nrow=2, labels="AUTO", label_x=0.05, label_size=large_font)
p <- plot_grid(ptmp, get_legend(p3 + theme(legend.position="bottom")), nrow=2, rel_heights = c(2, 0.1))
save_plot(paste0(folder, "IValue_mcf10.pdf", sep=""), p, base_aspect_ratio = 1.85, base_height=9)

