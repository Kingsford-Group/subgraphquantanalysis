library(ggplot2)
library(cowplot)
library(assertthat)
library(DESeq2)
library(tximport)

options(stringsAsFactors = FALSE)
args <- commandArgs(trailingOnly = FALSE)

codedir_split <- unlist( strsplit( sub("--file=", "", args[grep("--file=", args)]), "/") )
if (length(codedir_split) > 2) {
	        codedir <- paste0(codedir_split[1:(length(codedir_split)-2)], sep="/")
} else if (length(codedir_split) == 2) {
	        codedir <- "./"
} else if (length(codedir_split) == 1) {
	        codedir <- "../"
}
metafile <- paste0(codedir,"/data/Metadata_immune29.txt",sep="")

# output directory
idx <- grep("--args", args)
folder <- paste0(args[idx+1], "/Immune29/", sep="")
print(paste0("Meta data file = ", metafile, sep=""))
print(paste0("Output folder = ", folder, sep=""))

# read Meta file
meta <- read.table(metafile, header=F, sep=",")
meta <- meta[, c(1, 30)]
colnames(meta) <- c("ID", "celltype")

# read ivalue
iv <- read.table(paste0(folder, "IValue_mean_ext.txt", sep=""), header=F, sep="\t")
colnames(iv) <- c("Name", "Group1Mean", "Group2Mean", "IV", "IV25")

# DE detection on transcript level
files <- paste0(folder, meta$ID, "/quant.sf", sep="")
names(files) <- meta$ID
txi <- tximport(files, type = "salmon", txOut = TRUE)
dds <- DESeqDataSetFromTximport(txi, meta, ~celltype)
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
dds_gene <- DESeqDataSetFromTximport(txi, meta, ~celltype)
dds_gene <- DESeq(dds_gene)
res_gene <- results(dds_gene, alpha=0.01)
res_gene <- res_gene[!is.na(res_gene$padj),]
print(sum(res_gene$padj < 0.01, na.rm=TRUE))

related_trans <- res[rownames(res) %in% gene_trans_map[gene_trans_map$GENEID %in% rownames(res_gene[res_gene$padj < 0.01, ]), "TXNAME"], ]
print(nrow(related_trans[related_trans$padj < 0.01, ]))
print(related_trans[related_trans$padj < 0.01, "IV25"])

# plot example curve
t <- read.table(paste0(folder, "IValue_curve_example.txt", sep=""), header=F, sep="\t")
colnames(t) <- c("Name", "Group", "Sample", "reference_proportion", "lb", "ub")
df <- data.frame(Name=rep(t$Name, 2), Group=rep(t$Group, 2), Sample=rep(t$Sample, 2), reference_proportion=rep(t$reference_proportion, 2), Expression=c(t$lb, t$ub), Type=c(rep("lower bound", nrow(t)), rep("upper bound", nrow(t))))

large_font <- 18

tname <- "ENST00000258412.7"
# data frame
g <- "Effector memory CD8 T cells"
tmpdf <- data.frame(Group=g, reference_proportion=df[df$Name == tname & df$Sample == "Mean" & df$Type=="lower bound" & df$Group==g, "reference_proportion"], lb=df[df$Name == tname & df$Sample == "Mean" & df$Type=="lower bound" & df$Group==g, "Expression"], ub=df[df$Name == tname & df$Sample == "Mean" & df$Type=="upper bound" & df$Group==g, "Expression"])
g <- "Naive CD8 T cells"
tmpdf <- rbind(tmpdf, data.frame(Group=g, reference_proportion=df[df$Name == tname & df$Sample == "Mean" & df$Type=="lower bound" & df$Group==g, "reference_proportion"], lb=df[df$Name == tname & df$Sample == "Mean" & df$Type=="lower bound" & df$Group==g, "Expression"], ub=df[df$Name == tname & df$Sample == "Mean" & df$Type=="upper bound" & df$Group==g, "Expression"]) )
# plot
p1 <- ggplot(df[df$Name == tname & df$Sample == "Mean", ]) + 
	geom_ribbon_pattern(data=tmpdf, aes(x=1-reference_proportion, ymin=lb, ymax=ub, pattern=Group), alpha=0.5, fill='white', color="white") + 
	geom_line(aes(x=1-reference_proportion, y=Expression, group=interaction(Group, Sample, Type), linetype=Type), size=1.2) + 
	geom_vline(xintercept=1-0.04499289308663804, linetype="dotted") + geom_vline(xintercept=1-0.04687254015675888, linetype="dotdash") + theme_bw() + labs(title = paste0("Ranges of optima of ", tname, sep=""), x = "reference completeness", y = "normalized abundance") + 
	theme(legend.title=element_blank()) + theme(axis.text.x = element_text(size = large_font), axis.title.x = element_text(size = large_font), axis.text.y = element_text(size = large_font), axis.title.y = element_text(size = large_font), legend.text=element_text(size=large_font))

tname <- "ENST00000252725.9"
# data frame
g <- "Effector memory CD8 T cells"
tmpdf <- data.frame(Group=g, reference_proportion=df[df$Name == tname & df$Sample == "Mean" & df$Type=="lower bound" & df$Group==g, "reference_proportion"], lb=df[df$Name == tname & df$Sample == "Mean" & df$Type=="lower bound" & df$Group==g, "Expression"], ub=df[df$Name == tname & df$Sample == "Mean" & df$Type=="upper bound" & df$Group==g, "Expression"])
g <- "Naive CD8 T cells"
tmpdf <- rbind(tmpdf, data.frame(Group=g, reference_proportion=df[df$Name == tname & df$Sample == "Mean" & df$Type=="lower bound" & df$Group==g, "reference_proportion"], lb=df[df$Name == tname & df$Sample == "Mean" & df$Type=="lower bound" & df$Group==g, "Expression"], ub=df[df$Name == tname & df$Sample == "Mean" & df$Type=="upper bound" & df$Group==g, "Expression"]) )
# plot
p2 <- ggplot(df[df$Name == tname & df$Sample == "Mean", ]) + 
	geom_ribbon_pattern(data=tmpdf, aes(x=1-reference_proportion, ymin=lb, ymax=ub, pattern=Group), alpha=0.5, fill='white', color="white") + 
	geom_line(aes(x=1-reference_proportion, y=Expression, group=interaction(Group, Sample, Type), linetype=Type), size=1.2) + 
	geom_vline(xintercept=1-0.03342894511714156, linetype="dotted") + geom_vline(xintercept=1-0.08858369915708197, linetype="dotdash") + theme_bw() + labs(title = paste0("Ranges of optima of ", tname, sep=""), x = "reference completeness", y = "normalized abundance") + 
	theme(legend.title=element_blank()) + theme(axis.text.x = element_text(size = large_font), axis.title.x = element_text(size = large_font), axis.text.y = element_text(size = large_font), axis.title.y = element_text(size = large_font), legend.text=element_text(size=large_font))

tname <- "ENST00000619423.4"
# data frame
g <- "Effector memory CD8 T cells"
tmpdf <- data.frame(Group=g, reference_proportion=df[df$Name == tname & df$Sample == "Mean" & df$Type=="lower bound" & df$Group==g, "reference_proportion"], lb=df[df$Name == tname & df$Sample == "Mean" & df$Type=="lower bound" & df$Group==g, "Expression"], ub=df[df$Name == tname & df$Sample == "Mean" & df$Type=="upper bound" & df$Group==g, "Expression"])
g <- "Naive CD8 T cells"
tmpdf <- rbind(tmpdf, data.frame(Group=g, reference_proportion=df[df$Name == tname & df$Sample == "Mean" & df$Type=="lower bound" & df$Group==g, "reference_proportion"], lb=df[df$Name == tname & df$Sample == "Mean" & df$Type=="lower bound" & df$Group==g, "Expression"], ub=df[df$Name == tname & df$Sample == "Mean" & df$Type=="upper bound" & df$Group==g, "Expression"]) )
# plot
p3 <- ggplot(df[df$Name == tname & df$Sample == "Mean", ]) + 
	geom_ribbon_pattern(data=tmpdf, aes(x=1-reference_proportion, ymin=lb, ymax=ub, pattern=Group), alpha=0.5, fill='white', color="white") + 
	geom_line(aes(x=1-reference_proportion, y=Expression, group=interaction(Group, Sample, Type), linetype=Type), size=1.2) + 
	geom_vline(xintercept=1-0.27387597150657944, linetype="dotted") + geom_vline(xintercept=1-0.5877837534929233, linetype="dotdash") + theme_bw() + labs(title = paste0("Ranges of optima of ", tname, sep=""), x = "reference completeness", y = "normalized abundance") + 
	theme(legend.title=element_blank()) + theme(axis.text.x = element_text(size = large_font), axis.title.x = element_text(size = large_font), axis.text.y = element_text(size = large_font), axis.title.y = element_text(size = large_font), legend.text=element_text(size=large_font))

tname <- "ENST00000435064.5"
# data Frame
g <- "Effector memory CD8 T cells"
tmpdf <- data.frame(Group=g, reference_proportion=df[df$Name == tname & df$Sample == "Mean" & df$Type=="lower bound" & df$Group==g, "reference_proportion"], lb=df[df$Name == tname & df$Sample == "Mean" & df$Type=="lower bound" & df$Group==g, "Expression"], ub=df[df$Name == tname & df$Sample == "Mean" & df$Type=="upper bound" & df$Group==g, "Expression"])
g <- "Naive CD8 T cells"
tmpdf <- rbind(tmpdf, data.frame(Group=g, reference_proportion=df[df$Name == tname & df$Sample == "Mean" & df$Type=="lower bound" & df$Group==g, "reference_proportion"], lb=df[df$Name == tname & df$Sample == "Mean" & df$Type=="lower bound" & df$Group==g, "Expression"], ub=df[df$Name == tname & df$Sample == "Mean" & df$Type=="upper bound" & df$Group==g, "Expression"]) )
# plot
p4 <- ggplot(df[df$Name == tname & df$Sample == "Mean", ]) + 
	geom_ribbon_pattern(data=tmpdf, aes(x=1-reference_proportion, ymin=lb, ymax=ub, pattern=Group), alpha=0.5, fill='white', color="white") + 
	geom_line(aes(x=1-reference_proportion, y=Expression, group=interaction(Group, Sample, Type), linetype=Type), size=1.2) + 
	geom_vline(xintercept=1-0.5312576093476583, linetype="dotted") + geom_vline(xintercept=1-0.6897569071633847, linetype="dotdash") + theme_bw() + labs(title = paste0("Ranges of optima of ", tname, sep=""), x = "reference completeness", y = "normalized abundance") + 
	theme(legend.title=element_blank()) + theme(axis.text.x = element_text(size = large_font), axis.title.x = element_text(size = large_font), axis.text.y = element_text(size = large_font), axis.title.y = element_text(size = large_font), legend.text=element_text(size=large_font))

tname <- "ENST00000371733.7"
# data frame
g <- "Effector memory CD8 T cells"
tmpdf <- data.frame(Group=g, reference_proportion=df[df$Name == tname & df$Sample == "Mean" & df$Type=="lower bound" & df$Group==g, "reference_proportion"], lb=df[df$Name == tname & df$Sample == "Mean" & df$Type=="lower bound" & df$Group==g, "Expression"], ub=df[df$Name == tname & df$Sample == "Mean" & df$Type=="upper bound" & df$Group==g, "Expression"])
g <- "Naive CD8 T cells"
tmpdf <- rbind(tmpdf, data.frame(Group=g, reference_proportion=df[df$Name == tname & df$Sample == "Mean" & df$Type=="lower bound" & df$Group==g, "reference_proportion"], lb=df[df$Name == tname & df$Sample == "Mean" & df$Type=="lower bound" & df$Group==g, "Expression"], ub=df[df$Name == tname & df$Sample == "Mean" & df$Type=="upper bound" & df$Group==g, "Expression"]) )
# plot
p5 <- ggplot(df[df$Name == tname & df$Sample == "Mean", ]) + 
	geom_ribbon_pattern(data=tmpdf, aes(x=1-reference_proportion, ymin=lb, ymax=ub, pattern=Group), alpha=0.5, fill='white', color="white") + 
	geom_line(aes(x=1-reference_proportion, y=Expression, group=interaction(Group, Sample, Type), linetype=Type), size=1.2) + 
	geom_vline(xintercept=1-0.5755250675766428, linetype="dotted") + geom_vline(xintercept=1-0.8502603740959019, linetype="dotdash") + theme_bw() + labs(title = paste0("Ranges of optima of ", tname, sep=""), x = "reference completeness", y = "normalized abundance") + 
	theme(legend.title=element_blank()) + theme(axis.text.x = element_text(size = large_font), axis.title.x = element_text(size = large_font), axis.text.y = element_text(size = large_font), axis.title.y = element_text(size = large_font), legend.text=element_text(size=large_font))

p6 <- ggplot(as.data.frame(res[res$padj < 0.01, ])) + geom_histogram(aes(x=1-IV25)) + theme_cowplot() + labs(x = "reference completeness", title = "Number of DE transcripts that are unreliable \nunder the reference completeness parameter") + 
	theme(axis.text.x = element_text(size = large_font), axis.title.x = element_text(size = large_font), axis.text.y = element_text(size = large_font), axis.title.y = element_text(size = large_font), legend.text=element_text(size=large_font))

ptmp <- plot_grid(p1 + theme(legend.position="none"), p2 + theme(legend.position="none"), p3 + theme(legend.position="none"), p4 + theme(legend.position="none"), p5 + theme(legend.position="none"), p6, nrow=3, labels="AUTO", label_x=0.03, label_size=large_font)

p <- plot_grid(ptmp, get_legend(p1 + theme(legend.position="bottom")), nrow=2, rel_heights = c(3, 0.1))
save_plot(paste0(folder, "IValue_immune29.pdf", sep=""), p, base_aspect_ratio = 0.85, base_height=14)

