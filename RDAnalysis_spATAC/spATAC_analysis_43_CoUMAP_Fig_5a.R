rm(list = ls())

work_wd <- "/share/analysisdata/zhangx/CAESAR/SpATACME/"
setwd(work_wd)
source("functions.R")

figs_wd <- paste0(work_wd, "Figs/")

load("cols_cts_q.rda")

ana_wd <- paste0(work_wd, "spATACME/")
setwd(ana_wd)



################################################################################
# co-embedding plot
################################################################################
# find genes
library(Matrix)

sigList <- c(
    sc = lapply("sig_list_scRNA.rda", function(file) {
        load(paste0(ana_wd, file))
        sig_list[[1]]
    }),
    sp = lapply("df_sig_list_CoFASTu.rda", function(file) {
        load(paste0(ana_wd, file))
        df_sig_list
    })
)

ct_counts <- c(
    sc = lapply("seu_sc.rda", function(file) {
        load(paste0(ana_wd, file))
        table(seu_sc$merged_clusters)
    }),
    sp = lapply("seu.rda", function(file) {
        load(paste0(ana_wd, file))
        table(seu$pred_CoFASTu)
    })
)


ct_ratio <- lapply(ct_counts, function(counts) {
    counts / sum(counts)
})

cts <- setdiff(names(sigList[["sp"]]), "unassigned")

sg_overlap <- setNames(
    lapply(cts, function(ct) {
        sig <- lapply(sigList["sc"], function(sig_i) {
            sig_i_ct <- sig_i[[ct]]
            sig_i_ct <- sig_i_ct[sig_i_ct$expr.prop > 0.1, ]
            sig_i_ct$gene
        })
        gene_here <- Reduce(intersect, sig)

        mitochondrial_pattern <- "^(mt-|Mt-)"
        ribosomal_pattern <- "^(Rps|Rpl)"
        mitochondrial_genes <- grep(mitochondrial_pattern, gene_here, value = TRUE)
        ribosomal_genes <- grep(ribosomal_pattern, gene_here, value = TRUE)
        gene_here <- setdiff(gene_here, union(mitochondrial_genes, ribosomal_genes))
        if (length(gene_here) == 0) {
            return(NULL)
        }
        sig <- lapply(sig, function(sig_i_ct) {
            sig_i_ct[sig_i_ct %in% gene_here]
        })
        # id <- Reduce(`+`, lapply(sig, function(sig_i_ct) {
        #     match(gene_here, sig_i_ct)
        # }))
        # sig_overlap <- gene_here[order(id)[1:min(2, length(id))]]

        id2 <- Reduce(`+`, lapply(1:length(sig), function(i) {
            match(gene_here, sig[[i]]) * ct_ratio[[i]][ct]
        }))
        sig_overlap <- gene_here[order(id2)[1:min(2, length(id2))]]
        sig_overlap

        # lapply(sigList, function(sig_i) {
        #     sig_i[[ct]][sig_overlap, ]
        # })
    }),
    cts
)
save(sg_overlap, file = paste0(work_wd, "sg_overlap.rda"))

cols_coumap <- c(cols[cts], "unassigned" = "#737373")




pred <- "CAESARunasg"
celltypes <- names(cols_coumap)
Idents(seu) <- factor(seu[[pred]][, 1], levels = c(celltypes, "unassigned"))

cts_sig <- names(ct_counts[["sp"]][ct_counts[["sp"]] > 50])
cts_sig <- setdiff(cts_sig, "unassigned")

coumap_gene <- Reduce(rbind, lapply(
    cts_sig, function(ct) {
        gene <- sg_overlap[[ct]][1]
        sig_sp <- sigList[["sp"]][[ct]]
        sig_sp[match(gene, sig_sp$gene), ]
    }
))

seu <- CoUMAP(
    seu,
    reduction = "caesar", reduction.name = "CoUMAP",
    gene.set = coumap_gene
)


plot3 <- CoUMAP.plot(
    seu,
    reduction = "CoUMAP", gene_txtdata = df_coumap_gene,
    cols = c("gene" = "#000000", cols_coumap),
    pt_size = 0.8, pt_text_size = 4, alpha = 0.9
) + guides(
    shape = guide_legend(override.aes = list(size = 3)),
    color = guide_legend(ncol = 1, override.aes = list(size = 3))
)
ggsave(
    file = paste0(figs_wd, "CoUMAP_", coumap.name, "_", pred, "_E11_overlap.png"),
    plot = plot3, width = 9, height = 5, units = "in", dpi = 200
)




## seu_sc
Idents(seu_sc) <- seu_sc$merged_clusters
coumap.name.sc <- "UMAP"
cols_coumap_sc <- cols[names(sigList[["sc"]])]
seu_sc <- coembedding_umap_here_sc(
    seu_sc,
    reduction = "ncfm", reduction.name = coumap.name.sc,
    gene.set = unique(topnsig_overlap$gene)
)


plot1 <- CoUMAP.plot(
    seu_sc[, sample(seq_len(ncol(seu_sc)), ncol(seu_sc))],
    reduction = coumap.name.sc, gene_txtdata = topnsig_overlap,
    cols = c("gene" = "#000000", cols_coumap_sc),
    pt_size = 0.5, pt_text_size = 4, alpha = 0.8
) + guides(
    shape = guide_legend(override.aes = list(size = 3)),
    color = guide_legend(ncol = 1, override.aes = list(size = 3))
)
ggsave(
    file = paste0(figs_wd, "Step2_CoUMAP_", coumap.name.sc, "_scRNA_overlap.png"),
    plot = plot1, width = 12, height = 10, units = "in", dpi = 200
)














