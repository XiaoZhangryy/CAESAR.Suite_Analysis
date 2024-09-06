rm(list = ls())
library(Seurat)
message("Seurat version is ", packageVersion("Seurat"))
library(ProFAST)
library(ggplot2)
library(dplyr)
library(CelliD)
library(data.table)
library(CAESAR.Suite)

work_wd <- "/share/analysisdata/zhangx/CAESAR/SpATACME/"
dir.create(work_wd, showWarnings = FALSE)
figs_wd <- paste0(work_wd, "Figs/")
dir.create(figs_wd, showWarnings = FALSE)

setwd(work_wd)

load("cols_cts_q.rda")

ana_wd <- paste0(work_wd, "spATACME11/")
dir.create(ana_wd, showWarnings = FALSE)
setwd(ana_wd)

load("seu.rda")

################################################################################
# co-embedding plot
################################################################################
# find genes
library(Matrix)

sigList <- c(
    sc = lapply("sg_sc_List.rda", function(file) {
        load(paste0(work_wd, file))
        sg_sc_List
    }),
    sp = lapply("sg_list.rda", function(file) {
        load(paste0(ana_wd, file))
        sg_list
    })
)

ct_counts <- c(
    sc = lapply("seu_sc.rda", function(file) {
        load(paste0(work_wd, file))
        table(seu_sc$merged_clusters)
    }),
    sp = lapply("seu.rda", function(file) {
        load(paste0(ana_wd, file))
        table(seu$CAESARunasg)
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

        id2 <- Reduce(`+`, lapply(1:length(sig), function(i) {
            match(gene_here, sig[[i]]) * ct_ratio[[i]][ct]
        }))
        sig_overlap <- gene_here[order(id2)[1:min(2, length(id2))]]
        sig_overlap
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

topnsig_overlap <- Reduce(rbind, lapply(
    cts_sig, function(ct) {
        gene <- sg_overlap[[ct]][1]
        sig_sp <- sigList[["sp"]][[ct]]
        sig_sp[match(gene, sig_sp$gene), ]
    }
))

coumap.name <- "UMAP"
seu <- CoUMAP(
    seu,
    reduction = "caesar", reduction.name = coumap.name,
    gene.set = unique(topnsig_overlap$gene)
)


plot3 <- CoUMAP.plot(
    seu,
    reduction = coumap.name, gene_txtdata = topnsig_overlap,
    cols = c("gene" = "#000000", cols_coumap),
    pt_size = 0.8, pt_text_size = 4, alpha = 0.9
) + guides(
    shape = guide_legend(override.aes = list(size = 3)),
    color = guide_legend(ncol = 1, override.aes = list(size = 3))
)
ggsave(
    file = paste0(figs_wd, "Step2_CoUMAP_", coumap.name, "_", pred, "_E11_overlap.png"),
    plot = plot3, width = 9, height = 5, units = "in", dpi = 200
)


Idents(seu_sc) <- seu_sc$merged_clusters
coumap.name.sc <- "UMAP"
cols_coumap_sc <- cols[names(sigList[["sc"]])]
seu_sc <- CoUMAP(
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


