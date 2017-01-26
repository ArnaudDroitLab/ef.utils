library(eric.utils)
library(testthat)
library(ENCODExplorer)
load("C:/Users/UL/Desktop/Stage/GitProstate/sb_lab/eric.utils/R/sysdata.rda")

# Unit tests
guillimin.dir = system.file("extdata", "guillimin", package="eric.utils")
guillimin.meta = import_into_grl(input.dir=guillimin.dir, file.format="narrow", dir.type="mugqic")
stopifnot(length(guillimin.meta)==3)
stopifnot(all(unlist(lapply(guillimin.meta, length))==100))

guillimin.nometa = import_into_grl(input.dir=guillimin.dir, file.format="narrow", dir.type="mugqic", discard.metadata=TRUE)
stopifnot(ncol(GenomicRanges::mcols(guillimin.nometa[[1]]))==0)

plain.dir = system.file("extdata", "plain", package="eric.utils")
plain.narrow = import_into_grl(input.dir=file.path(plain.dir, "narrow"), file.format="narrow")
plain.bed = import_into_grl(input.dir=file.path(plain.dir, "bed"), file.format="bed", file.ext = "narrowPeak")

######################################
# Test build_intersect
######################################
intersect.test  = build_intersect(plain.narrow)
# General tests
    test_that("build_intersect_type", expect_type(intersect.test, "list"))
    # Test the form of returned Regions
    test_that("build_intersect_regions type", expect_type(intersect.test$Regions, "S4"))
    test_that("build_intersect_regions class", expect_s4_class(intersect.test$Regions, "GRanges"))
    test_that("build_intersect_regions size", expect_true(length(plain.narrow@unlistData@elementMetadata@listData$name)
                                                          >= length(intersect.test$Regions)))
    # Test the content of returned Regions
    test_that("build_intersect_regions", expect_equal(intersect.test$Regions, GenomicRanges::reduce(unlist(plain.narrow))))
    # Test the form of returned Matrix
    test_that("build_intersect_matrix", expect_type(intersect.test$Matrix, "double"))
    test_that("build_intersect_ncols matrix", expect_equal(length(unique(plain.narrow@partitioning@NAMES)),
                                                           ncol(intersect.test$Matrix)))
    test_that("build_intersect_colnames matrix", expect_equal(unique(plain.narrow@partitioning@NAMES),
                                                           colnames(intersect.test$Matrix)))
    test_that("build_intersect_nrow matrix", expect_equal(length(intersect.test$Regions), nrow(intersect.test$Matrix)))
    # Test the content of returned Matrix
    test_that("build_intersect_no empty row", expect_equal(0, sum(apply(intersect.test$Matrix, 1, sum) == 0)))
    test_that("build_intersect_correct overlap", expect_equal(length(plain.narrow@unlistData@elementMetadata@listData$name),
                                                              sum(apply(intersect.test$Matrix, 1, sum))))
    # Test the form and content of returned Names
    test_that("build_intersect_names", expect_type(intersect.test$Names, "character"))
    test_that("build_intersect_names length", expect_equal(ncol(intersect.test$Matrix), length(intersect.test$Names)))
    test_that("build_intersect_names", expect_equal(colnames(intersect.test$Matrix), intersect.test$Names))
    # Test the form and content of returned Length
    test_that("build_intersect_length", expect_type(intersect.test$Length, "integer"))
    test_that("build_intersect_length", expect_equal(ncol(intersect.test$Matrix), intersect.test$Length))

# Specific tests

######################################
# Test pairwise_overlap
######################################
overlap.test = pairwise_overlap(intersect.test)
# General tests
      # Test the form of returned object
      test_that("pairwise_overlap", expect_type(overlap.test, "double"))
      test_that("pairwise_overlap_ncol", expect_equal(length(intersect.test$Names), ncol(overlap.test)))
      test_that("pairwise_overlap_nrow", expect_equal(ncol(overlap.test), nrow(overlap.test)))
      # Test the form of returned object
      test_that("pairwise_overlap_colnames", expect_equal(intersect.test$Names, colnames(overlap.test)))
      test_that("pairwise_overlap_rownames", expect_equal(intersect.test$Names, row.names(overlap.test)))
      test_that("pairwise_overlap_values", expect_true(sum(apply(overlap.test, 1, sum)) < ncol(overlap.test)^2))

# Specific tests
      stopifnot(overlap.test[1,1] == 1)
      stopifnot(ncol(overlap.test)==3 && nrow(overlap.test)==3)
      pairwise_overlap(intersect.test, filename = "pairwise_overlap.txt")
      test_that("pairwise_overlap_output file", expect_true(file.exists("pairwise_overlap.txt")))
      file.remove("pairwise_overlap.txt")

#####################################
# Test intersect_overlap
#####################################
intersect_overlap.test <- intersect_overlap(intersect.test)
# General tests
      # Test the form of returned object
      test_that("intersect_overlap", expect_type(intersect_overlap.test, "S4"))
      test_that("intersect_overlap_class", expect_s4_class(intersect_overlap.test, "GRanges"))

# Specific tests
      stopifnot(length(intersect_overlap(intersect.test))==38)
      stopifnot(length(intersect_overlap(intersect.test, names=c("GM12878_CDK9_Rep2_2WCE_Narrow_Peaks", "GM12878_MED1_Rep1_2WCE_Narrow_Peaks")))==44)
      stopifnot(length(intersect_overlap(intersect.test, names=c("GM12878_CDK9_Rep2_2WCE_Narrow_Peaks", "GM12878_MED1_Rep1_2WCE_Narrow_Peaks"), exclusive=TRUE))==6)

      # When considering all factors, exclusive should be == inclusive.
      stopifnot(intersect_overlap(intersect.test, exclusive=TRUE)==intersect_overlap(intersect.test))

#####################################
# Test intersect_venn_plot
#####################################
# Specific tests
      intersect.test.2 <- intersect.test
      intersect.test.2$Length <- 7
      test_that("intersect_venn_plot_length > 5", expect_error(intersect_venn_plot(intersect.test.2)))

      intersect_venn_plot(intersect.test, file="Test.tiff")
      stopifnot(file.exists("Test.tiff"))
      file.remove("Test.tiff")

#####################################
# Test select_annotations
#####################################
# General tests
      test_that("select_annotations_incorrect genome", expect_error(select_annotations("genome")))
      test_that("select_annotations_type", expect_type(select_annotations("hg19"), "list"))
      test_that("select_annotations_length", expect_equal(6, length(select_annotations("hg19"))))

# Specific tests
      select_annotations.test <- select_annotations("hg19")
      data(PWMLogn.hg19.MotifDb.Hsap, package="PWMEnrich.Hsapiens.background")
      test_that("select_annotations_TxDb", expect_equal(TxDb.Hsapiens.UCSC.hg19.knownGene::TxDb.Hsapiens.UCSC.hg19.knownGene,
                                                        select_annotations.test[[1]]))
      test_that("select_annotations_OrgDbStr", expect_equal("org.Hs.eg.db", select_annotations.test[[2]]))
      test_that("select_annotations_OrgDb", expect_equal(org.Hs.eg.db::org.Hs.eg.db, select_annotations.test[[3]]))
      test_that("select_annotations_PWMBG", expect_equal(PWMLogn.hg19.MotifDb.Hsap, select_annotations.test[[4]]))
      test_that("select_annotations_BSGenome", expect_equal(BSgenome.Hsapiens.UCSC.hg19::BSgenome.Hsapiens.UCSC.hg19, select_annotations.test[[5]]))
      test_that("select_annotations_KEGG", expect_equal(hs.keggs, select_annotations.test[[6]]))

#####################################
# Test annotate_venn_center
#####################################
# General tests
      annotate_venn_center.test <- annotate_venn_center(intersect.test, select_annotations.test)
      test_that("annotate_venn_center_type", expect_type(intersect_overlap.test, "S4"))
      #test_that("annotate_venn_center_class", expect_s4_class(intersect_overlap.test, "csAnno"))

# Specific tests
      annotate_venn_center.test <- annotate_venn_center(intersect.test, select_annotations.test, "annotate_venn_center.txt")
      test_that("annotate_venn_center_output file", expect_true(file.exists("annotate_venn_center.txt")))
      file.remove("annotate_venn_center.txt")

  #####################################
# Test annotate_venn_exclusive
#####################################
# General tests
      annotate_venn_exclusive.test <- annotate_venn_exclusive(intersect.test, select_annotations.test)
      test_that("annotate_venn_exclusive_type", expect_type(annotate_venn_exclusive.test, "list"))

# Specific tests

#####################################
# Test build_intersect_all
#####################################
# General tests
      build_intersect_all.test <- build_intersect_all(plain.narrow, select_annotations.test, "test")
      test_that("build_intersect_all_type", expect_type(build_intersect_all.test$Regions, "S4"))
      test_that("build_intersect_all_class", expect_s4_class(build_intersect_all.test$Regions, "GRanges"))
      test_that("build_intersect_all_values", expect_equal(build_intersect_all.test, intersect.test))

# Specific tests
      test_that("build_intersect_all_output file", expect_true(dir.exists("output")))
      unlink("output")

#####################################
# Test annotate_region
#####################################
# General tests
      annotate_region.test <- annotate_region(intersect.test$Regions, select_annotations.test)
      test_that("annotate_region_type", expect_type(annotate_region.test, "S4"))
      test_that("annotate_region_class", expect_s4_class(annotate_region.test, "csAnno"))
      test_that("annotate_region_length", expect_equal(length(intersect.test$Regions@seqnames), length(annotate_region.test@anno@seqnames)))

# Specific tests
      annotate_region.df <- as.data.frame(annotate_region.test)
      test_that("annotate_region_lign1 symbol", expect_equal("LOC100288069", annotate_region.df$SYMBOL[1]))
      test_that("annotate_region_lign10 transcriptID", expect_equal("uc001acj.4", annotate_region.df$transcriptId[10]))
      test_that("annotate_region_lign106, symbol", expect_equal("KLHL21", annotate_region.df$SYMBOL[106]))

      annotate_region.test <- annotate_region(intersect.test$Regions, select_annotations.test, filename = "annotate.txt")
      test_that("annotate_region_output file", expect_true(file.exists("annotate.txt")))
      file.remove("annotate.txt")

#####################################
# Test motif.enrichement
#####################################
# General tests
      motif_enrichment.test <- motif_enrichment(intersect.test$Regions, select_annotations.test)
      test_that("motif_enrichment_type", expect_type(motif_enrichment.test, "list"))
      test_that("motif_enrichment_length", expect_equal(3, length(motif_enrichment.test)))
      test_that("motif_enrichment_Report", expect_equal(length(select_annotations.test$PWMBG$pwms),
                                                        length(motif_enrichment.test$Report@pwms)))

# Specific tests
      test_that("motif_enrichment_Regions", expect_equal(176, length(motif_enrichment.test$Region)))

      dir.create("motif_enrichment", recursive = TRUE)
      motif_enrichment.test <- motif_enrichment(intersect.test$Regions, select_annotations.test, file.label = "motif_enrichment/")
      test_that("motif_enrichment_output file", expect_true(file.exists("motif_enrichment/ MotifEnrichment.txt")))
      unlink("motif_enrichment", recursive = TRUE)

#####################################
# Test get_promoters
#####################################
# General tests
      genes <- unique(annotate_region.df$geneId)
      get_promoters.test <- get_promoters(genes, select_annotations.test)
      test_that("get_promoters_type", expect_type(get_promoters.test, "S4"))
      test_that("get_promoters_class", expect_s4_class(get_promoters.test, "GRanges"))
      test_that("get_promoters_length", expect_equal(length(genes), length(get_promoters.test@seqnames)))

# Specific tests


#####################################
# Test motif_enrichment_genes
#####################################
# General tests
      motif_enrichment_genes.test <- motif_enrichment_genes(genes, select_annotations.test)
      test_that("motif_enrichment_genes_type", expect_type(motif_enrichment_genes.test, "list"))
      test_that("motif_enrichment_genes_length", expect_equal(3, length(motif_enrichment_genes.test)))
      test_that("motif_enrichment_genes_Report", expect_equal(length(select_annotations.test$PWMBG$pwms),
                                                        length(motif_enrichment_genes.test$Report@pwms)))
      test_that("motif_enrichment_genes_Regions", expect_equal(length(genes), length(motif_enrichment_genes.test$Region)))

# Specific tests
      dir.create("motif_enrichment_genes", recursive = TRUE)
      motif_enrichment_genes.test <- motif_enrichment_genes(genes, select_annotations.test, file.label = "motif_enrichment_genes/")
      test_that("motif_enrichment_output file", expect_true(file.exists("motif_enrichment_genes/ MotifEnrichment.txt")))
      unlink("motif_enrichment_genes", recursive = TRUE)

#####################################
# Test kegg_enrichment
#####################################
# General tests
      kegg_enrichment.test <- kegg_enrichment(genes, select_annotations.test)
      test_that("kegg_enrichment_type", expect_type(kegg_enrichment.test, "list"))
      test_that("kegg_enrichment_ncol", expect_equal(8, ncol(kegg_enrichment.test)))

# Specific tests
      test_that("kegg_enrichment_nrow", expect_equal(length(names(select_annotations.test$KEGG$kg.sets[select_annotations.test$KEGG$sigmet.idx])),
                                                     nrow(kegg_enrichment.test)))
      kegg_enrichment.test.diseases <- kegg_enrichment(genes, select_annotations.test, diseases = TRUE)
      test_that("kegg_enrichment_nrow diseases", expect_equal(length(names(select_annotations.test$KEGG$kg.sets[select_annotations.test$KEGG$dise.idx])),
                                                              nrow(kegg_enrichment.test.diseases)))

      kegg_enrichment.test <- kegg_enrichment(genes, select_annotations.test, filename = "kegg_enrichment.txt")
      test_that("kegg_enrichment_output file", expect_true(file.exists("kegg_enrichment.txt")))

#####################################
# Test gene_from_regions
#####################################
# General tests
      gene_from_regions.test <- gene_from_regions(intersect.test$Region, select_annotations.test)
      test_that("gene_from_regions_type", expect_type(gene_from_regions.test, "character"))
      test_that("gene_from_regions_length", expect_true(length(gene_from_regions.test) <= length(motif_enrichment.test$Region@seqnames)))
      test_that("gene_from_regions_values by default", expect_equal(length(gene_from_regions.test), sum(gene_from_regions.test %in% genes)))

      test_that("gene_from_regions_invalid argument", expect_error(gene_from_regions(intersect.test$Regions, select_annotations.test, region.types = "")))

      gene_from_regions.test <- gene_from_regions(intersect.test$Regions, select_annotations.test, region.types = c("Promoter", "Exon"))

      gene_from_regions.test <- gene_from_regions(intersect.test$Regions, select_annotations.test, region.types = "All")
      test_that("gene_from_regions_values all", expect_equal(gene_from_regions.test, genes))

# Specific tests

####################################
# Test kegg_enrichment_regions
####################################
# General tests
      kegg_enrichment_regions.test <- kegg_enrichment_regions(intersect.test$Regions, select_annotations.test)
      test_that("kegg_enrichment_regions_type", expect_type(kegg_enrichment_regions.test, "list"))
      test_that("kegg_enrichment_regions_ncol", expect_equal(8, ncol(kegg_enrichment_regions.test)))

# Specific tests
      test_that("kegg_enrichment_regions_nrow", expect_equal(length(names(select_annotations.test$KEGG$kg.sets[select_annotations.test$KEGG$sigmet.idx])),
                                                     nrow(kegg_enrichment_regions.test)))
      test_that("kegg_enrichment_regions_correspondance with kegg_enrichment_regions",
                expect_equal(kegg_enrichment.test, kegg_enrichment_regions.test))
      kegg_enrichment_regions.test.diseases <- kegg_enrichment_regions(intersect.test$Regions, select_annotations.test, diseases = TRUE)

####################################
# Test characterize_region
####################################
# General tests
      characterize_region.test <- characterize_region(intersect.test$Regions, select_annotations.test, output.dir = NULL)
      test_that("characterize_region_type", expect_type(characterize_region.test, "list"))
      test_that("characterize_region_length", expect_equal(4, length(characterize_region.test)))
      test_that("characterize_region_annotation type", expect_type(characterize_region.test[["Annotation"]], "S4"))
      test_that("characterize_region_annotation class", expect_s4_class(characterize_region.test[["Annotation"]], "csAnno"))
      test_that("characterize_region_motif", expect_type(characterize_region.test[["Motif"]], "list"))
      test_that("characterize_region_motif length", expect_equal(3, length(characterize_region.test[["Motif"]])))
      test_that("characterize_region_KEGG.sig", expect_type(characterize_region.test[["KEGG.sig"]], "list"))
      test_that("characterize_region_KEGG.dis", expect_type(characterize_region.test[["KEGG.dis"]], "list"))

# Specific tests
      test_that("characterize_region_annotation values", expect_equal(characterize_region.test[["Annotation"]],
                                                                      annotate_region.test))
      test_that("characterize_region_motif values", expect_equal(characterize_region.test[["Motif"]],
                                                                   motif_enrichment.test))
      test_that("characterize_region_kegg values", expect_equal(characterize_region.test[["KEGG.sig"]],
                                                                  kegg_enrichment_regions.test))
      test_that("characterize_region_kegg diseases values", expect_equal(characterize_region.test[["KEGG.dis"]],
                                                                           kegg_enrichment_regions.test.diseases))

      characterize_region.test <- characterize_region(intersect.test$Regions, select_annotations.test)
      test_that("characterize_region_output file annotation", expect_true(file.exists("output/Annotations.txt")))
      test_that("characterize_region_output file motif", expect_true(file.exists("output/Motifs MotifEnrichment.txt")))
      test_that("characterize_region_output file KEGG.sig", expect_true(file.exists("output/KEGG signalisation and metabolism.txt")))
      test_that("characterize_region_output file KEGG.dis", expect_true(file.exists("output/KEGG diseases.txt")))
      unlink("output", recursive = TRUE)

####################################
# Test characterize_gene_set
####################################
# General tests
      characterize_gene_set.test <- characterize_gene_set(genes, select_annotations.test, output.dir = NULL)
      test_that("characterize_gene_set_type", expect_type(characterize_gene_set.test, "list"))
      test_that("characterize_gene_set_length", expect_equal(3, length(characterize_gene_set.test)))
      test_that("characterize_gene_set_motif", expect_type(characterize_gene_set.test[["Motif"]], "list"))
      test_that("characterize_gene_set_motif length", expect_equal(3, length(characterize_gene_set.test[["Motif"]])))
      test_that("characterize_gene_set_KEGG.sig", expect_type(characterize_gene_set.test[["KEGG.sig"]], "list"))
      test_that("characterize_gene_set_KEGG.dis", expect_type(characterize_gene_set.test[["KEGG.dis"]], "list"))

# Specific tests
      test_that("characterize_gene_set_motif values", expect_equal(characterize_gene_set.test[["Motif"]],
                                                                   motif_enrichment_genes.test))
      test_that("characterize_gene_set_kegg values", expect_equal(characterize_gene_set.test[["KEGG.sig"]],
                                                                  kegg_enrichment.test))
      test_that("characterize_gene_set_kegg diseases values", expect_equal(characterize_gene_set.test[["KEGG.dis"]],
                                                                  kegg_enrichment.test.diseases))

      characterize_gene_set.test <- characterize_gene_set(genes, select_annotations.test)
      test_that("characterize_gene_set_output file motif", expect_true(file.exists("output/Motifs MotifEnrichment.txt")))
      test_that("characterize_gene_set_output file KEGG.sig", expect_true(file.exists("output/KEGG signalisation and metabolism.txt")))
      test_that("characterize_gene_set_output file KEGG.dis", expect_true(file.exists("output/KEGG diseases.txt")))
      unlink("output", recursive = TRUE)

####################################
# Test import.and.discard.metadata
####################################
# General tests


# Specific tests

####################################
# Test import_into_grl
####################################
# General tests


# Specific tests
      test_that("import_into_grl_nrow", expect_equal(length(plain.bed@unlistData@seqnames), 100*3))

####################################
# Test default.download.filter
####################################
# General tests


# Specific tests

####################################
# Test default_download_filter_rna
####################################
# General tests


# Specific tests

####################################
# Test download_encode_chip
####################################
# General tests
      download_encode_chip.test <- download_encode_chip("MCF-7", "hg19")

# Specific tests

####################################
# Test download_encode_rna
####################################
# General tests
      download_encode_rna.test <- download_encode_rna("MCF-7", "hg19")
      test_that("download_encode_rna_values", expect_false(0 == sum(apply(download_encode_rna.test, 1, sum))))

# Specific tests

