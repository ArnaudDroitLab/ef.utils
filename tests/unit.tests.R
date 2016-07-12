library(eric.utils)
library(testthat)
library(ENCODExplorer)
load("C:/Users/UL/Desktop/Stage/GitProstate/sb_lab/eric.utils/R/sysdata.rda")

# Unit tests
guillimin.dir = system.file("extdata", "guillimin", package="eric.utils")
guillimin.meta = import.into.grl(input.dir=guillimin.dir, file.format="narrow", dir.type="mugqic")
stopifnot(length(guillimin.meta)==3)
stopifnot(all(unlist(lapply(guillimin.meta, length))==100))

guillimin.nometa = import.into.grl(input.dir=guillimin.dir, file.format="narrow", dir.type="mugqic", discard.metadata=TRUE)
stopifnot(ncol(GenomicRanges::mcols(guillimin.nometa[[1]]))==0)

plain.dir = system.file("extdata", "plain", package="eric.utils")
plain.narrow = import.into.grl(input.dir=file.path(plain.dir, "narrow"), file.format="narrow")
plain.bed = import.into.grl(input.dir=file.path(plain.dir, "bed"), file.format="bed", file.ext = "narrowPeak")

######################################
# Test build.intersect
######################################
intersect.test  = build.intersect(plain.narrow)
# General tests
    test_that("build.intersect_type", expect_type(intersect.test, "list"))
    # Test the form of returned Regions
    test_that("build.intersect_regions type", expect_type(intersect.test$Regions, "S4"))
    test_that("build.intersect_regions class", expect_s4_class(intersect.test$Regions, "GRanges"))
    test_that("build.intersect_regions size", expect_true(length(plain.narrow@unlistData@elementMetadata@listData$name)
                                                          >= length(intersect.test$Regions)))
    # Test the content of returned Regions
    test_that("build.intersect_regions", expect_equal(intersect.test$Regions, GenomicRanges::reduce(unlist(plain.narrow))))
    # Test the form of returned Matrix
    test_that("build.intersect_matrix", expect_type(intersect.test$Matrix, "double"))
    test_that("build.intersect_ncols matrix", expect_equal(length(unique(plain.narrow@partitioning@NAMES)),
                                                           ncol(intersect.test$Matrix)))
    test_that("build.intersect_colnames matrix", expect_equal(unique(plain.narrow@partitioning@NAMES),
                                                           colnames(intersect.test$Matrix)))
    test_that("build.intersect_nrow matrix", expect_equal(length(intersect.test$Regions), nrow(intersect.test$Matrix)))
    # Test the content of returned Matrix
    test_that("build.intersect_no empty row", expect_equal(0, sum(apply(intersect.test$Matrix, 1, sum) == 0)))
    test_that("build.intersect_correct overlap", expect_equal(length(plain.narrow@unlistData@elementMetadata@listData$name),
                                                              sum(apply(intersect.test$Matrix, 1, sum))))
    # Test the form and content of returned Names
    test_that("build.intersect_names", expect_type(intersect.test$Names, "character"))
    test_that("build.intersect_names length", expect_equal(ncol(intersect.test$Matrix), length(intersect.test$Names)))
    test_that("build.intersect_names", expect_equal(colnames(intersect.test$Matrix), intersect.test$Names))
    # Test the form and content of returned Length
    test_that("build.intersect_length", expect_type(intersect.test$Length, "integer"))
    test_that("build.intersect_length", expect_equal(ncol(intersect.test$Matrix), intersect.test$Length))

# Specific tests

######################################
# Test pairwise.overlap
######################################
overlap.test = pairwise.overlap(intersect.test)
# General tests
      # Test the form of returned object
      test_that("pairwise.overlap", expect_type(overlap.test, "double"))
      test_that("pairwise.overlap_ncol", expect_equal(length(intersect.test$Names), ncol(overlap.test)))
      test_that("pairwise.overlap_nrow", expect_equal(ncol(overlap.test), nrow(overlap.test)))
      # Test the form of returned object
      test_that("pairwise.overlap_colnames", expect_equal(intersect.test$Names, colnames(overlap.test)))
      test_that("pairwise.overlap_rownames", expect_equal(intersect.test$Names, row.names(overlap.test)))
      test_that("pairwise.overlap_values", expect_true(sum(apply(overlap.test, 1, sum)) < ncol(overlap.test)^2))

# Specific tests
      stopifnot(overlap.test[1,1] == 1)
      stopifnot(ncol(overlap.test)==3 && nrow(overlap.test)==3)
      pairwise.overlap(intersect.test, filename = "pairwise.overlap.txt")
      test_that("pairwise.overlap_output file", expect_true(file.exists("pairwise.overlap.txt")))
      file.remove("pairwise.overlap.txt")

#####################################
# Test intersect.overlap
#####################################
intersect.overlap.test <- intersect.overlap(intersect.test)
# General tests
      # Test the form of returned object
      test_that("intersect.overlap", expect_type(intersect.overlap.test, "S4"))
      test_that("intersect.overlap_class", expect_s4_class(intersect.overlap.test, "GRanges"))

# Specific tests
      stopifnot(length(intersect.overlap(intersect.test))==38)
      stopifnot(length(intersect.overlap(intersect.test, names=c("GM12878_CDK9_Rep2_2WCE_Narrow_Peaks", "GM12878_MED1_Rep1_2WCE_Narrow_Peaks")))==44)
      stopifnot(length(intersect.overlap(intersect.test, names=c("GM12878_CDK9_Rep2_2WCE_Narrow_Peaks", "GM12878_MED1_Rep1_2WCE_Narrow_Peaks"), exclusive=TRUE))==6)

      # When considering all factors, exclusive should be == inclusive.
      stopifnot(intersect.overlap(intersect.test, exclusive=TRUE)==intersect.overlap(intersect.test))

#####################################
# Test plot.intersect.venn
#####################################
# Specific tests
      intersect.test.2 <- intersect.test
      intersect.test.2$Length <- 7
      test_that("plot.intersect.venn_length > 5", expect_error(plot.intersect.venn(intersect.test.2)))

      plot.intersect.venn(intersect.test, file="Test.tiff")
      stopifnot(file.exists("Test.tiff"))
      file.remove("Test.tiff")

#####################################
# Test select.annotations
#####################################
# General tests
      test_that("select.annotations_incorrect genome", expect_error(select.annotations("genome")))
      test_that("select.annotations_type", expect_type(select.annotations("hg19"), "list"))
      test_that("select.annotations_length", expect_equal(6, length(select.annotations("hg19"))))

# Specific tests
      select.annotations.test <- select.annotations("hg19")
      data(PWMLogn.hg19.MotifDb.Hsap, package="PWMEnrich.Hsapiens.background")
      test_that("select.annotations_TxDb", expect_equal(TxDb.Hsapiens.UCSC.hg19.knownGene::TxDb.Hsapiens.UCSC.hg19.knownGene,
                                                        select.annotations.test[[1]]))
      test_that("select.annotations_OrgDbStr", expect_equal("org.Hs.eg.db", select.annotations.test[[2]]))
      test_that("select.annotations_OrgDb", expect_equal(org.Hs.eg.db::org.Hs.eg.db, select.annotations.test[[3]]))
      test_that("select.annotations_PWMBG", expect_equal(PWMLogn.hg19.MotifDb.Hsap, select.annotations.test[[4]]))
      test_that("select.annotations_BSGenome", expect_equal(BSgenome.Hsapiens.UCSC.hg19::BSgenome.Hsapiens.UCSC.hg19, select.annotations.test[[5]]))
      test_that("select.annotations_KEGG", expect_equal(hs.keggs, select.annotations.test[[6]]))

#####################################
# Test annotate.venn.center
#####################################
# General tests
      annotate.venn.center.test <- annotate.venn.center(intersect.test, select.annotations.test)
      test_that("annotate.venn.center_type", expect_type(intersect.overlap.test, "S4"))
      #test_that("annotate.venn.center_class", expect_s4_class(intersect.overlap.test, "csAnno"))

# Specific tests
      annotate.venn.center.test <- annotate.venn.center(intersect.test, select.annotations.test, "annotate.venn.center.txt")
      test_that("annotate.venn.center_output file", expect_true(file.exists("annotate.venn.center.txt")))
      file.remove("annotate.venn.center.txt")

  #####################################
# Test annotate.venn.exclusive
#####################################
# General tests
      annotate.venn.exclusive.test <- annotate.venn.exclusive(intersect.test, select.annotations.test)
      test_that("annotate.venn.exclusive_type", expect_type(annotate.venn.exclusive.test, "list"))

# Specific tests

#####################################
# Test build.intersect.all
#####################################
# General tests
      build.intersect.all.test <- build.intersect.all(plain.narrow, select.annotations.test, "test")
      test_that("build.intersect.all_type", expect_type(build.intersect.all.test$Regions, "S4"))
      test_that("build.intersect.all_class", expect_s4_class(build.intersect.all.test$Regions, "GRanges"))
      test_that("build.intersect.all_values", expect_equal(build.intersect.all.test, intersect.test))

# Specific tests
      test_that("build.intersect.all_output file", expect_true(dir.exists("output")))
      unlink("output")

#####################################
# Test annotate.region
#####################################
# General tests
      annotate.region.test <- annotate.region(intersect.test$Regions, select.annotations.test)
      test_that("annotate.region_type", expect_type(annotate.region.test, "S4"))
      test_that("annotate.region_class", expect_s4_class(annotate.region.test, "csAnno"))
      test_that("annotate.region_length", expect_equal(length(intersect.test$Regions@seqnames), length(annotate.region.test@anno@seqnames)))

# Specific tests
      annotate.region.df <- as.data.frame(annotate.region.test)
      test_that("annotate.region_lign1 symbol", expect_equal("LOC100288069", annotate.region.df$SYMBOL[1]))
      test_that("annotate.region_lign10 transcriptID", expect_equal("uc001acj.4", annotate.region.df$transcriptId[10]))
      test_that("annotate.region_lign106, symbol", expect_equal("KLHL21", annotate.region.df$SYMBOL[106]))

      annotate.region.test <- annotate.region(intersect.test$Regions, select.annotations.test, filename = "annotate.txt")
      test_that("annotate.region_output file", expect_true(file.exists("annotate.txt")))
      file.remove("annotate.txt")

#####################################
# Test motif.enrichement
#####################################
# General tests
      motif.enrichment.test <- motif.enrichment(intersect.test$Regions, select.annotations.test)
      test_that("motif.enrichment_type", expect_type(motif.enrichment.test, "list"))
      test_that("motif.enrichment_length", expect_equal(3, length(motif.enrichment.test)))
      test_that("motif.enrichment_Report", expect_equal(length(select.annotations.test$PWMBG$pwms),
                                                        length(motif.enrichment.test$Report@pwms)))

# Specific tests
      test_that("motif.enrichment_Regions", expect_equal(176, length(motif.enrichment.test$Region)))

      dir.create("motif.enrichment", recursive = TRUE)
      motif.enrichment.test <- motif.enrichment(intersect.test$Regions, select.annotations.test, file.label = "motif.enrichment/")
      test_that("motif.enrichment_output file", expect_true(file.exists("motif.enrichment/ MotifEnrichment.txt")))
      unlink("motif.enrichment", recursive = TRUE)

#####################################
# Test get.promoters
#####################################
# General tests
      genes <- unique(annotate.region.df$geneId)
      get.promoters.test <- get.promoters(genes, select.annotations.test)
      test_that("get.promoters_type", expect_type(get.promoters.test, "S4"))
      test_that("get.promoters_class", expect_s4_class(get.promoters.test, "GRanges"))
      test_that("get.promoters_length", expect_equal(length(genes), length(get.promoters.test@seqnames)))

# Specific tests


#####################################
# Test motif.enrichment.genes
#####################################
# General tests
      motif.enrichment.genes.test <- motif.enrichment.genes(genes, select.annotations.test)
      test_that("motif.enrichment.genes_type", expect_type(motif.enrichment.genes.test, "list"))
      test_that("motif.enrichment.genes_length", expect_equal(3, length(motif.enrichment.genes.test)))
      test_that("motif.enrichment.genes_Report", expect_equal(length(select.annotations.test$PWMBG$pwms),
                                                        length(motif.enrichment.genes.test$Report@pwms)))
      test_that("motif.enrichment.genes_Regions", expect_equal(length(genes), length(motif.enrichment.genes.test$Region)))

# Specific tests
      dir.create("motif.enrichment.genes", recursive = TRUE)
      motif.enrichment.genes.test <- motif.enrichment.genes(genes, select.annotations.test, file.label = "motif.enrichment.genes/")
      test_that("motif.enrichment_output file", expect_true(file.exists("motif.enrichment.genes/ MotifEnrichment.txt")))
      unlink("motif.enrichment.genes", recursive = TRUE)

#####################################
# Test kegg.enrichment
#####################################
# General tests
      kegg.enrichment.test <- kegg.enrichment(genes, select.annotations.test)
      test_that("kegg.enrichment_type", expect_type(kegg.enrichment.test, "list"))
      test_that("kegg.enrichment_ncol", expect_equal(8, ncol(kegg.enrichment.test)))

# Specific tests
      test_that("kegg.enrichment_nrow", expect_equal(length(names(select.annotations.test$KEGG$kg.sets[select.annotations.test$KEGG$sigmet.idx])),
                                                     nrow(kegg.enrichment.test)))
      kegg.enrichment.test.diseases <- kegg.enrichment(genes, select.annotations.test, diseases = TRUE)
      test_that("kegg.enrichment_nrow diseases", expect_equal(length(names(select.annotations.test$KEGG$kg.sets[select.annotations.test$KEGG$dise.idx])),
                                                              nrow(kegg.enrichment.test.diseases)))

      kegg.enrichment.test <- kegg.enrichment(genes, select.annotations.test, filename = "kegg.enrichment.txt")
      test_that("kegg.enrichment_output file", expect_true(file.exists("kegg.enrichment.txt")))

#####################################
# Test gene.from.regions
#####################################
# General tests
      gene.from.regions.test <- gene.from.regions(intersect.test$Region, select.annotations.test)
      test_that("gene.from.regions_type", expect_type(gene.from.regions.test, "character"))
      test_that("gene.from.regions_length", expect_true(length(gene.from.regions.test) <= length(motif.enrichment.test$Region@seqnames)))
      test_that("gene.from.regions_values by default", expect_equal(length(gene.from.regions.test), sum(gene.from.regions.test %in% genes)))

      test_that("gene.from.regions_invalid argument", expect_error(gene.from.regions(intersect.test$Regions, select.annotations.test, region.types = "")))

      gene.from.regions.test <- gene.from.regions(intersect.test$Regions, select.annotations.test, region.types = c("Promoter", "Exon"))

      gene.from.regions.test <- gene.from.regions(intersect.test$Regions, select.annotations.test, region.types = "All")
      test_that("gene.from.regions_values all", expect_equal(gene.from.regions.test, genes))

# Specific tests

####################################
# Test kegg.enrichment.regions
####################################
# General tests
      kegg.enrichment.regions.test <- kegg.enrichment.regions(intersect.test$Regions, select.annotations.test)
      test_that("kegg.enrichment.regions_type", expect_type(kegg.enrichment.regions.test, "list"))
      test_that("kegg.enrichment.regions_ncol", expect_equal(8, ncol(kegg.enrichment.regions.test)))

# Specific tests
      test_that("kegg.enrichment.regions_nrow", expect_equal(length(names(select.annotations.test$KEGG$kg.sets[select.annotations.test$KEGG$sigmet.idx])),
                                                     nrow(kegg.enrichment.regions.test)))
      test_that("kegg.enrichment.regions_correspondance with kegg.enrichment.regions",
                expect_equal(kegg.enrichment.test, kegg.enrichment.regions.test))
      kegg.enrichment.regions.test.diseases <- kegg.enrichment.regions(intersect.test$Regions, select.annotations.test, diseases = TRUE)

####################################
# Test characterize.region
####################################
# General tests
      characterize.region.test <- characterize.region(intersect.test$Regions, select.annotations.test, output.dir = NULL)
      test_that("characterize.region_type", expect_type(characterize.region.test, "list"))
      test_that("characterize.region_length", expect_equal(4, length(characterize.region.test)))
      test_that("characterize.region_annotation type", expect_type(characterize.region.test[["Annotation"]], "S4"))
      test_that("characterize.region_annotation class", expect_s4_class(characterize.region.test[["Annotation"]], "csAnno"))
      test_that("characterize.region_motif", expect_type(characterize.region.test[["Motif"]], "list"))
      test_that("characterize.region_motif length", expect_equal(3, length(characterize.region.test[["Motif"]])))
      test_that("characterize.region_KEGG.sig", expect_type(characterize.region.test[["KEGG.sig"]], "list"))
      test_that("characterize.region_KEGG.dis", expect_type(characterize.region.test[["KEGG.dis"]], "list"))

# Specific tests
      test_that("characterize.region_annotation values", expect_equal(characterize.region.test[["Annotation"]],
                                                                      annotate.region.test))
      test_that("characterize.region_motif values", expect_equal(characterize.region.test[["Motif"]],
                                                                   motif.enrichment.test))
      test_that("characterize.region_kegg values", expect_equal(characterize.region.test[["KEGG.sig"]],
                                                                  kegg.enrichment.regions.test))
      test_that("characterize.region_kegg diseases values", expect_equal(characterize.region.test[["KEGG.dis"]],
                                                                           kegg.enrichment.regions.test.diseases))

      characterize.region.test <- characterize.region(intersect.test$Regions, select.annotations.test)
      test_that("characterize.region_output file annotation", expect_true(file.exists("output/Annotations.txt")))
      test_that("characterize.region_output file motif", expect_true(file.exists("output/Motifs MotifEnrichment.txt")))
      test_that("characterize.region_output file KEGG.sig", expect_true(file.exists("output/KEGG signalisation and metabolism.txt")))
      test_that("characterize.region_output file KEGG.dis", expect_true(file.exists("output/KEGG diseases.txt")))
      unlink("output", recursive = TRUE)

####################################
# Test characterize.gene.set
####################################
# General tests
      characterize.gene.set.test <- characterize.gene.set(genes, select.annotations.test, output.dir = NULL)
      test_that("characterize.gene.set_type", expect_type(characterize.gene.set.test, "list"))
      test_that("characterize.gene.set_length", expect_equal(3, length(characterize.gene.set.test)))
      test_that("characterize.gene.set_motif", expect_type(characterize.gene.set.test[["Motif"]], "list"))
      test_that("characterize.gene.set_motif length", expect_equal(3, length(characterize.gene.set.test[["Motif"]])))
      test_that("characterize.gene.set_KEGG.sig", expect_type(characterize.gene.set.test[["KEGG.sig"]], "list"))
      test_that("characterize.gene.set_KEGG.dis", expect_type(characterize.gene.set.test[["KEGG.dis"]], "list"))

# Specific tests
      test_that("characterize.gene.set_motif values", expect_equal(characterize.gene.set.test[["Motif"]],
                                                                   motif.enrichment.genes.test))
      test_that("characterize.gene.set_kegg values", expect_equal(characterize.gene.set.test[["KEGG.sig"]],
                                                                  kegg.enrichment.test))
      test_that("characterize.gene.set_kegg diseases values", expect_equal(characterize.gene.set.test[["KEGG.dis"]],
                                                                  kegg.enrichment.test.diseases))

      characterize.gene.set.test <- characterize.gene.set(genes, select.annotations.test)
      test_that("characterize.gene.set_output file motif", expect_true(file.exists("output/Motifs MotifEnrichment.txt")))
      test_that("characterize.gene.set_output file KEGG.sig", expect_true(file.exists("output/KEGG signalisation and metabolism.txt")))
      test_that("characterize.gene.set_output file KEGG.dis", expect_true(file.exists("output/KEGG diseases.txt")))
      unlink("output", recursive = TRUE)

####################################
# Test import.and.discard.metadata
####################################
# General tests


# Specific tests

####################################
# Test import.into.grl
####################################
# General tests


# Specific tests
      test_that("import.into.grl_nrow", expect_equal(length(plain.bed@unlistData@seqnames), 100*3))

####################################
# Test default.download.filter
####################################
# General tests


# Specific tests

####################################
# Test default.download.filter.rna
####################################
# General tests


# Specific tests

####################################
# Test download.encode.chip
####################################
# General tests
      download.encode.chip.test <- download.encode.chip("MCF-7", "hg19")

# Specific tests

####################################
# Test download.encode.rna
####################################
# General tests
      download.encode.rna.test <- download.encode.rna("MCF-7", "hg19")
      test_that("download.encode.rna_values", expect_false(0 == sum(apply(download.encode.rna.test, 1, sum))))

# Specific tests

