library("phyloseq") #1.24.2
library(dplyr)
library("ggplot2") #3.0.0
library("tidyverse") #1.2.1
library(qiime2R) #0.99.1
library(microbiomeSeq) #0.1
library(DESeq2) #1.20.0
library("yingtools2") #install_github("ying14/yingtools2")) 0.0.0.89
library(data.table)

#R3.5.1


#required inputs:
#1. metadata in utf8
#2. feature table
#3. taxonomy
#4. rooted tree

outputPath = "H:\\OneDrive\\School_Work\\_Beluga\\beluga_all_paired\\correlations"
dir.create(file.path(outputPath), showWarnings = FALSE)
setwd(file.path(outputPath))

#metadata = read_tsv("D:\\OneDrive\\School_Work\\_Beluga\\beluga_all_merged\\csv\\metadata.tsv")
metadata = read_tsv("H:\\OneDrive\\School_Work\\_Beluga\\metadatas\\metadata.ansi.contaminantsOnly.short.quartiles.tsv")
#colnames(metadata)

featureTable = read_qza("H:\\OneDrive\\School_Work\\_Beluga\\beluga_all_paired\\scp\\beluga_all_paired.featuretable.qza", "temp")
#colnames(featureTable$data)

taxonomy = read_qza("H:\\OneDrive\\School_Work\\_Beluga\\beluga_all_paired\\scp\\beluga_all_paired.taxonomy.vsearch.qza", "temp")
#colnames(taxonomy$data)
taxtable = taxonomy$data %>% as.tibble() %>% separate(Taxon, sep = ";", c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")) #convert the table into a tabular split version
#View(taxtable)
tree = read_qza("H:\\OneDrive\\School_Work\\_Beluga\\beluga_all_paired\\scp\\beluga_all_paired.rootedTree.qza", "temp")
#View(metadata)
#shannon = read_qza("D:\\\\OneDrive\\\\Desktop\\\\_Beluga\\\\Beluga\\\workflow.p1p2\\phyloseq\\shannon_vector.qza", "temp")
#pcoa = read_qza("D:\\OneDrive\\Desktop\\_Beluga\\Beluga\\workflow.p1p2\\phyloseq\\unweighted_unifrac_pcoa_results.qza", "temp")

phyloseqObj = phyloseq(
    otu_table(featureTable$data, taxa_are_rows = T),
    phy_tree(tree$data),
    tax_table(as.data.frame(taxtable) %>% select(-Confidence) %>% column_to_rownames("Feature.ID") %>% as.matrix()),
    sample_data(metadata %>% as.data.frame() %>% column_to_rownames("sampleid"))
    )


phyloseqObj_percent = transform_sample_counts(phyloseqObj, function(x) 100 * x / sum(x))
#taxa_are_rows(phyloseqObj)
tPhyloseqObj = t(phyloseqObj)
#taxa_are_rows(tPhyloseqObj)
#tax_table(tPhyloseqObj)

#lefse
lefseWrapper <- function(phy, class, subclass = NA, subject = NA, anova.alpha = 0.05,
                wilcoxon.alpha = 0.05, lda.cutoff = 2, wilcoxon.within.subclass = FALSE,
                one.against.one = FALSE, mult.test.correction = 0, make.lefse.plots = FALSE,
                by_otus = FALSE, levels = rank_names(phy)) {
    keepvars <- c(class, subclass, subject, "sample")
    keepvars <- unique(keepvars[!is.na(keepvars)])

    samp <- get.samp(phy)[, keepvars]

    if (by_otus) {
        otu <- get.otu.melt(phy, sample_data = FALSE)

        otu.levels <- otu %>%
            mutate(taxon = otu) %>%
            group_by(sample, taxon) %>%
            summarize(pctseqs = sum(pctseqs)) %>%
            mutate(taxon = gsub(" ", "_", taxon))
    } else {
        otu <- get.otu.melt(phy, sample_data = FALSE)

        otu.list <- lapply(1:length(levels), function(i) {
            lvls <- levels[1:i]
            lvl <- levels[i]
            otu.level <- otu
            otu.level$taxon <- do.call(paste, c(lapply(lvls, function(l) otu[[l]]), sep = "|"))
            otu.level$rank <- lvl
            otu.level2 <- otu.level %>%
              group_by(sample, taxon, rank) %>%
              summarize(pctseqs = sum(pctseqs)) %>%
              ungroup()
            return(otu.level2)
        })
        otu.levels <- bind_rows(otu.list) %>%
            mutate(taxon = gsub(" ", "_", taxon))
    }
    otu.tbl <- otu.levels %>%
          dcast(sample ~ taxon, value.var = "pctseqs", fill = 0) %>%
          left_join(samp, by = "sample") %>%
          select_(.dots = c(keepvars, lazyeval::interp(~everything())))
    if (is.na(subject) | subject != "sample") {
        otu.tbl <- otu.tbl %>% select(-sample)
    }
    tbl <- otu.tbl %>% t()
    write.table(tbl, "lefse.txt", quote = FALSE, sep = "\t",
                    col.names = FALSE)
    opt.class <- paste("-c", which(keepvars %in% class))
    opt.subclass <- ifelse(is.na(subclass), "", paste("-s", which(keepvars %in% subclass)))
    opt.subject <- ifelse(is.na(subject), "", paste("-u", which(keepvars %in% subject)))
    format.command <- paste("bash -c \"./format_input.py lefse.txt lefse.in", opt.class, opt.subclass, opt.subject, "-o 1000000\"")
    print(format.command)
    system(format.command)
    lefse.command <- paste("bash -c \"./run_lefse.py lefse.in lefse.res",
                               "-a", anova.alpha, "-w", wilcoxon.alpha, "-l", lda.cutoff,
                               "-e", as.numeric(wilcoxon.within.subclass), "-y", as.numeric(one.against.one),
                               "-s", mult.test.correction, "\"")
    print(lefse.command)
    system(lefse.command)
    print("Wrote lefse.res")
    lefse.out <- read.table("lefse.res", header = FALSE, sep = "\t") %>%
          rename(taxon = V1, log.max.pct = V2, direction = V3, lda = V4, p.value = V5)
    if (make.lefse.plots) {
        system("bash -c \"./plot_res.py lefse.res lefse_lda.png\"")
        print("Wrote lefse_lda.png")
        system("bash -c \"./plot_cladogram.py lefse.res lefse_clado.pdf --format pdf\"")
        print("Wrote lefse_clado.pdf")
    }
    return(lefse.out)
}

setwd("H:\\OneDrive\\ProjectMicrobiome\\lefse")
lefse.tbl <- lefseWrapper(phyloseqObj, class = "Penta_Quartile", subclass = "Sex", make.lefse.plots = TRUE)


#plot_bar(phyloseqObj_percent, fill = "Order", x = "SampleSource", facet_grid = ~Sex) + geom_bar(aes(fill = Order, width = 0.75), stat = "identity", position = "stack")

#plot_tree(phyloseqObj, color="Phylum")
#plot_bar(phyloseqObj, fill = "Order")

#begin useful functions
#####env correlation

plot_taxa_env <- function(df) {
    p <- ggplot2::ggplot(aes(x = Env, y = Taxa, fill = Correlation), data = df)
    p <- p + ggplot2::geom_tile() + scale_fill_gradient2(low = "#2C7BB6", mid = "white", high = "#D7191C")
    p <- p + ggplot2::theme(axis.text.x = element_text(size = 20, angle = 90, hjust = 1, vjust = 0.5))
    p <- p + ggplot2::theme(axis.text.y = element_text(size = 20, angle = 0, hjust = 1, vjust = 0.5))
    p <- p + ggplot2::geom_text(aes(label = Significance), color = "black", size = 8) + labs(y = NULL, x = NULL)
    p <- p + ggplot2::facet_grid(. ~ Env, drop = TRUE, scale = "free", space = "free_x")
    p <- p + ggplot2::xlab("Groups")
    p <- p + ggplot2::theme(strip.background = element_rect(fill = "white"))
    return(p)
}

#phyloseqObj_phylum = taxa_level(tPhyloseqObj, "Genus")
#colnames(tax_table(tPhyloseqObj))

taxa_level <- function(physeq, which_level) {
    if (taxa_are_rows(physeq)) {
        physeq <- t(physeq)
    }
    #physeq = tPhyloseqObj
    #which_level = "Genus"
    OTU <- otu_table(physeq)
    SAM <- sample_data(physeq)
    OTU_taxonomy <- tax_table(physeq)
    new_abund_table <- NULL
    if (which_level == "Otus") {
        OTU_tree <- phy_tree(physeq)
        new_abund_table <- OTU
    } else {
        colNum = grep(which_level, colnames(OTU_taxonomy))
        subsetDF = data.frame(na.omit(unique(OTU_taxonomy[, (colNum - 1):colNum])))
        subsetDF = mutate(subsetDF, name = rownames(subsetDF))
        subsetDF = mutate(subsetDF, merged = paste(trimws(subsetDF[, 1]), " ", trimws(subsetDF[, 2])))
        subsetDF$merged = gsub(paste0("D_", as.character(colNum - 2), "__"), "", subsetDF$merged)
        subsetDF$merged = gsub(paste0("D_", as.character(colNum - 1), "__"), "", subsetDF$merged)
        subsetDF$merged[subsetDF$merged == "   "] = "Unclassified"
        subsetDF = subset(subsetDF, subsetDF$merged != "Unclassified")
        #        subsetDF = subset(subsetDF, unique(subsetDF$merged))
        #subsetDF$merged[subsetDF$merged == "NA   NA"] = "Unclassified"
        #subsetDF[, 2][subsetDF[, 2] == "NA"] = "Unclassified"
        rownames(subsetDF) = subsetDF$name
        list = subsetDF[, 2]
        list2 = subsetDF$merged
        #print(list2)

        #list <- na.omit(unique(OTU_taxonomy[, which_level]))

        #View(list2)
        new_abund_table <- NULL
        index = 0
        for (i in list) {
            index = index + 1
            #print(i)
            rt <- na.omit(rownames(OTU_taxonomy)[OTU_taxonomy[, which_level] == i])
            #print(rt)
            tmp <- data.frame(rowSums(OTU[, rt]))
            #print(tmp)
            colnames(tmp) = list2[index]
            #if (i == "") { colnames(tmp) <- c("__Unknowns__") } else { colnames(tmp) <- list2[index] }
            if (is.null(new_abund_table)) { new_abund_table <- tmp } else { new_abund_table <- cbind(tmp, new_abund_table) }
            }
        #View(new_abund_table)
        #colnames(new_abund_table)
    }
    OTU <- as.data.frame(as(new_abund_table, "matrix"))
    #print(OTU)
    #View(OTU)
    #Convert the data to phyloseq format
    OTU = otu_table(as.matrix(OTU), taxa_are_rows = FALSE)
    TAX = tax_table(as.matrix(OTU_taxonomy))
    SAM = sample_data(SAM)

    #reconstruct the phyloseq object
    physeq <- NULL
    if (which_level == "Otus") {
        physeq <- merge_phyloseq(phyloseq(OTU, TAX), SAM, midpoint(OTU_tree))
    } else {
        physeq <- merge_phyloseq(phyloseq(OTU), SAM)
    }
    return(physeq)
}

tables.correlate <- function(table1, table2, groups = NULL, method) {
    df <- NULL
    for (i in colnames(table1)) {
        for (j in colnames(table2)) {
            if (!is.null(groups)) {
                for (k in unique(groups)) {
                    a <- table1[groups == k, i, drop = F]
                    b <- table2[groups == k, j, drop = F]
                    tmp <- c(i, j, cor(a[complete.cases(b),], b[complete.cases(b),], use = "everything", method = method), cor.test(a[complete.cases(b),], b[complete.cases(b),], method = method)$p.value, k)
                    if (is.null(df)) { df <- tmp } else { df <- rbind(df, tmp) }
                    }
            }
            else {
                a <- table1[, i, drop = F]
                b <- table2[, j, drop = F]
                tmp <- c(i, j, cor(a[complete.cases(b),], b[complete.cases(b),], use = "everything", method = method), cor.test(a[complete.cases(b),], b[complete.cases(b),], method = method)$p.value)
                if (is.null(df)) { df <- tmp } else { df <- rbind(df, tmp) }
                }
        }
    }
    df <- data.frame(row.names = NULL, df)
    return(df)
}

# df is a data frame
p.adjust.cor <- function(df, adjustment = 1, padjust.method = "BH") {
    if (adjustment == 1) {
        df$AdjPvalue <- df$Pvalue
    } else if (adjustment == 2) {
        for (i in unique(df$Env)) {
            for (j in unique(df$Type)) {
                sel <- df$Env == i & df$Type == j
                df$AdjPvalue[sel] <- p.adjust(df$Pvalue[sel], method = padjust.method)
            }
        }
    } else if (adjustment == 3) {
        for (i in unique(df$Taxa)) {
            for (j in unique(df$Type)) {
                sel <- df$Taxa == i & df$Type == j
                df$AdjPvalue[sel] <- p.adjust(df$Pvalue[sel], method = padjust.method)
            }
        }
    } else if (adjustment == 4) {
        for (i in unique(df$Taxa)) {
            sel <- df$Taxa == i
            df$AdjPvalue[sel] <- p.adjust(df$Pvalue[sel], method = padjust.method)
        }
    } else if (adjustment == 5) {
        for (i in unique(df$Env)) {
            sel <- df$Env == i
            df$AdjPvalue[sel] <- p.adjust(df$Pvalue[sel], method = padjust.method)
        }
    }
    return(df)
}

make_env_cor_graph <- function(obj, level) {
    #level = "Family"
    #obj = tPhyloseqObj
    print(level)
    phyloseqObj_phylum = taxa_level(obj, level)
    #corr = taxa.env.correlation(tPhyloseqObj, NULL, method = "pearson", padjust.method = "BH")

    #grouping_column = "BDE7"
    abund_table <- otu_table(phyloseqObj_phylum)
    meta_table <- data.frame(sample_data(phyloseqObj_phylum))
    #groups <- meta_table[, grouping_column]
    #groups

    mt_env <- meta_table[, sapply(meta_table, is.numeric)]
    abund_table_filt <- abund_table[rownames(mt_env),]
    abund_table_filt <- abund_table_filt[, order(colSums(abund_table_filt), decreasing = TRUE)]
    taxa_list <- colnames(abund_table_filt)
    taxa_list <- taxa_list[!grepl("Unknown", taxa_list)]
    abund_table_filt <- data.frame(abund_table_filt[, colnames(abund_table_filt) %in% taxa_list])
    #df <- tables.correlate(abund_table_filt, mt_env, groups, "pearson")

    table1 = abund_table_filt
    table2 = mt_env

    method = "pearson"
    groups = NULL
    df <- NULL
    for (i in colnames(table1)) {
        for (j in colnames(table2)) {
            if (!is.null(groups)) {
                for (k in unique(groups)) {
                    a <- table1[groups == k, i, drop = F]
                    #print(a)
                    b <- table2[groups == k, j, drop = F]
                    #print(b)
                
                    cora = a[complete.cases(b),]
                    corb = b[complete.cases(b),]
                    #print(cora)
                    #print(corb)
                    #tmp <- c(i, j, cor(a[complete.cases(b),], b[complete.cases(b),], use = "everything", method = method), cor.test(a[complete.cases(b),], b[complete.cases(b),], method = method)$p.value, k)
                    #if (is.null(df)) { df <- tmp } else { df <- rbind(df, tmp) }
                }
            }
            else {
                #print(i)
                #print(j)
                a <- table1[, i, drop = F]
                b <- table2[, j, drop = F]
                #print(a)
                #print(b)
                tmp <- c(i, j, cor(a[complete.cases(b),], b[complete.cases(b),], use = "everything", method = method), cor.test(a[complete.cases(b),], b[complete.cases(b),], method = method)$p.value)
                if (is.null(df)) { df <- tmp } else { df <- rbind(df, tmp) }
                }
        }
    }

    View(df)

    df <- data.frame(row.names = NULL, df)
    #df
    #View(df)
    colnames(df) <- c("Taxa", "Env", "Correlation", "Pvalue")
    df$Pvalue <- as.numeric(as.character(df$Pvalue))
    df$Correlation <- as.numeric(as.character(df$Correlation))
    df$AdjPvalue <- rep(0, dim(df)[1])
    df <- p.adjust.cor(df, 1, "BH")
    df$Significance <- cut(df$AdjPvalue, breaks = c(-Inf, 0.001, 0.01, 0.05, Inf), label = c("***", "**", "*", ""))
    df <- df[complete.cases(df),]

    write.table(df, file = paste0("envCorrelation.", level, ".All.tsv"), quote = FALSE, sep = '\t')
    plot_taxa_env(df)
    ggsave(paste0("envCorrelation.", level, ".All.png"), last_plot(), units = "in", width = 15 + 0.5 * ncol(table2), height = 2 + 0.5 * ncol(table1), dpi = 160, limitsize = FALSE)

    sig = c(0.001, 0.01, 0.05)

    for (i in 1:length(sig)) {
        print(sig[i])
        #i = 3
        df.clean = df %>% filter(AdjPvalue <= sig[i])
        df.sig = df %>% filter(Taxa %in% df.clean$Taxa)
        if (nrow(df.sig) > 0) {
            write.table(df.sig, file = paste0("envCorrelation.", level, ".SignificantAt", as.character(sig[i]), ".tsv"), quote = FALSE, sep = '\t')
            plot_taxa_env(df.sig)
            ggsave(paste0("envCorrelation.", level, ".SignificantAt", as.character(sig[i]), ".png"), last_plot(), units = "in", width = 15 + 0.5 * ncol(table2), height = 2 + 0.5 * nrow(unique(df.sig["Taxa"])), dpi = 320, limitsize = FALSE, device = "png")
        }
    }
}

level = c("Phylum", "Class", "Order", "Family", "Genus", "Species")
level = c("Family")
for (i in 1:length(level)) {
    make_env_cor_graph(tPhyloseqObj, level[i])
}

#end env cor


##begin test functions
#subset = prune_taxa(rownames(sigtab), phyloseqObj_percent) #   subset_samples(phyloseqObj_percent, BDE49 == " 11 ")

#plot_heatmap(subset, method = NULL, sample.label = "BDE7", sample.order = "BDE7", taxa.label = "Family")
#ggsave(paste0(cell, ".png"), last_plot())


#deseq1obj = phyloseq_to_deseq2(phyloseqObj, ~ BDE7)
#gm_mean = function(x, na.rm = TRUE) {
    #exp(sum(log(x[x > 0]), na.rm = na.rm) / length(x))
#}
#geoMeans = apply(counts(deseq1obj), 1, gm_mean)
#deseq1obj = estimateSizeFactors(deseq1obj, geoMeans = geoMeans)
#deseq1obj = DESeq(deseq1obj, fitType = "local")

#res = results(deseq1obj, cooksCutoff = FALSE)
#res
#alpha = 0.05
#sigtab = res[which(res$padj < alpha),]
#sigtab
#sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(phyloseqObj)[rownames(sigtab),], "matrix"))


#theme_set(theme_bw())
#sigtabgen = subset(sigtab, !is.na(Genus))
## Phylum order
#x = tapply(sigtabgen$log2FoldChange, sigtabgen$Phylum, function(x) max(x))
#x = sort(x, TRUE)
#sigtabgen$Phylum = factor(as.character(sigtabgen$Phylum), levels = names(x))
## Genus order
#x = tapply(sigtabgen$log2FoldChange, sigtabgen$Genus, function(x) max(x))
#x = sort(x, TRUE)
#sigtabgen$Genus = factor(as.character(sigtabgen$Genus), levels = names(x))
#ggplot(sigtabgen, aes(y = Genus, x = log2FoldChange, color = Phylum)) +
  #geom_vline(xintercept = 0.0, color = "gray", size = 0.5) +
  #geom_point(size = 6) +
  #theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust = 0.5))



#taxa_are_rows(phyloseqObj)
#phyloseqObj = t(phyloseqObj)
#phyloseqObj

#p = plot_anova_diversity(phyloseqObj, method = c("richness", "simpson", "shannon"), grouping_column = "Sex", pValueCutoff = 0.05)
#p = plot_richness(phyloseqObj, x="Weightanalyzed(mg)", color="Weightanalyzed(mg)", measures = "Chaos1", title = "alpha div")

#print(p)

#physeq <- normalise_data(phyloseqObj, norm.method = "relative")
#taxa_level
#p <- plot_taxa(physeq, grouping_column = "Weightanalyzed(mg)", method = "hellinger", number.taxa = 21, filename = NULL)
#print(p)

#p <- plot_anova_env(physeq, grouping_column = "Sex", select.variables = c("BDE7", "Weightanalyzed(mg)"))

#columns = read.table(text = 'TotalBDE')
##columns = read.table(text='Lipides	BDE7	BDE10	BDE17	BDE28orPBT	BDE47	BDE49	BDE66	BDE77	BDE85	BDE99	BDE100	BDE126	BDE138	BDE139	BDE140	BDE153	BDE154orBB153	BDE183orDec604	BDE184	BDE196	BDE197or204	BDE201	BDE203	BDE207	BDE208	BDE209	Penta	Octa	Deca	percent_Penta	percent_Octa	percent_Deca	percent_Other	PBDE_tot	Dec604CB	BEHTBPorsynDP	antiDP	DP_tot	fanti	HBB	PBEB	HFRs	Dechloranes', header=FALSE)
#columns
#dir.create("sig", showWarnings = FALSE)
#dir.create("sig\\phylum", showWarnings = FALSE)
#dir.create("sig\\family", showWarnings = FALSE)
#dir.create("sig\\genus", showWarnings = FALSE)
#dir.create("all", showWarnings = FALSE)
#dir.create("all\\phylum", showWarnings = FALSE)
#dir.create("all\\family", showWarnings = FALSE)
#dir.create("all\\genus", showWarnings = FALSE)
#print(ncol(columns))
#for (col in 1:ncol(columns)) {
    #cell = paste0("~", as.character(columns[1, col]))
    #cell = as.formula(cell)
    #View(cell)

    #cellName = as.character(columns[1, col])
    #print(cellName)

    #deseq1obj = phyloseq_to_deseq2(phyloseqObj, cell)
    #gm_mean = function(x, na.rm = TRUE) {
        #exp(sum(log(x[x > 0]), na.rm = na.rm) / length(x))
    #}
    #geoMeans = apply(counts(deseq1obj), 1, gm_mean)
    #deseq1obj = estimateSizeFactors(deseq1obj, geoMeans = geoMeans)
    #deseq1obj = DESeq(deseq1obj, fitType = "local")

    #res = results(deseq1obj, cooksCutoff = FALSE)
    #alpha = 0.05
    #sigtab = res[which(res$padj < alpha),]
    #nrow(sigtab)
    #if (nrow(sigtab) > 0) {
        #sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(phyloseqObj)[rownames(sigtab),], "matrix"))
        #subset = prune_taxa(rownames(sigtab), phyloseqObj_percent)
        #plot_heatmap(subset, method = NULL, sample.label = cellName, sample.order = cellName, taxa.label = "Phylum")
        #ggsave(paste0("sig\\phylum\\_", cellName, ".phylum.png"), last_plot(), units = "in", width = 15, height = 15, dpi = 320)
        #plot_heatmap(subset, method = NULL, sample.label = cellName, sample.order = cellName, taxa.label = "Family")
        #ggsave(paste0("sig\\family\\_", cellName, ".family.png"), last_plot(), units = "in", width = 15, height = 15, dpi = 320)
        #plot_heatmap(subset, method = NULL, sample.label = cellName, sample.order = cellName, taxa.label = "Genus")
        #ggsave(paste0("sig\\genus\\_", cellName, ".genus.png"), last_plot(), units = "in", width = 15, height = 15, dpi = 320)
        #theme_set(theme_bw())
        #sigtabgen = subset(sigtab, !is.na(Genus))
        ## Phylum order
        #x = tapply(sigtabgen$log2FoldChange, sigtabgen$Phylum, function(x) max(x))
        #x = sort(x, TRUE)
        #sigtabgen$Phylum = factor(as.character(sigtabgen$Phylum), levels = names(x))
        ## Genus order
        #x = tapply(sigtabgen$log2FoldChange, sigtabgen$Genus, function(x) max(x))
        #x = sort(x, TRUE)
        #sigtabgen$Genus = factor(as.character(sigtabgen$Genus), levels = names(x))
        #ggplot(sigtabgen, aes(y = Genus, x = log2FoldChange, color = Phylum)) +
        #geom_vline(xintercept = 0.0, color = "gray", size = 0.5) +
        #geom_point(size = 6) +
        #theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust = 0.5))
        #ggsave(paste0("sig\\_", cellName, ".diff.png"), last_plot(), units = "in", width = 15, height = 15, dpi = 320)

    #}
    #tax_table(phyloseqObj_percent)
    #plot_heatmap(phyloseqObj_percent, method = NULL, sample.label = cellName, sample.order = cellName, taxa.label = "Phylum")
    #ggsave(paste0("all\\phylum\\", cellName, ".phylum.png"), last_plot(), units = "in", width = 15, height = 15, dpi = 320)
    #plot_heatmap(phyloseqObj_percent, method = NULL, sample.label = cellName, sample.order = cellName, taxa.label = "Family")
    #ggsave(paste0("all\\family\\", cellName, ".family.png"), last_plot(), units = "in", width = 15, height = 15, dpi = 320)
    #plot_heatmap(phyloseqObj_percent, method = NULL, sample.label = cellName, sample.order = cellName, taxa.label = "Genus", taxa.order = "")
    #ggsave(paste0("all\\genus\\", cellName, ".genus.png"), last_plot(), units = "in", width = 15, height = 15, dpi = 320)
#}