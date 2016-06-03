library(eric.utils)

# Unit tests
guillimin.dir = system.file("extdata", "guillimin", package="eric.utils")
guillimin.meta = import.into.grl(input.dir=guillimin.dir, file.format="narrow", dir.type="mugqic")
stopifnot(length(guillimin.meta)==3)
stopifnot(all(unlist(lapply(guillimin.meta, length))==100))

guillimin.nometa = import.into.grl(input.dir=guillimin.dir, file.format="narrow", dir.type="mugqic", discard.metadata=TRUE)
stopifnot(ncol(GenomicRanges::mcols(guillimin.nometa[[1]]))==0)

plain.dir = system.file("extdata", "plain", package="eric.utils")
plain.narrow = import.into.grl(input.dir=file.path(plain.dir, "narrow"), file.format="narrow")
plain.bed = import.into.grl(input.dir=file.path(plain.dir, "bed"), file.format="bed")

intersect.test  = build.intersect(plain.narrow)
overlap.test = pairwise.overlap(intersect.test)
stopifnot(overlap.test[1,1] == 1)
stopifnot(ncol(overlap.test)==3 && nrow(overlap.test)==3)

stopifnot(length(intersect.overlap(intersect.test)==38))
stopifnot(length(intersect.overlap(intersect.test, names=c("GM12878_CDK9_Rep2_2WCE_Narrow_Peaks", "GM12878_MED1_Rep1_2WCE_Narrow_Peaks")))==44)
stopifnot(length(intersect.overlap(intersect.test, names=c("GM12878_CDK9_Rep2_2WCE_Narrow_Peaks", "GM12878_MED1_Rep1_2WCE_Narrow_Peaks"), exclusive=TRUE))==6)

# When considering all factors, exclusive should be == inclusive.
stopifnot(intersect.overlap(intersect.test, exclusive=TRUE)==intersect.overlap(intersect.test))

plot.intersect.venn(intersect.test, file="Test.tiff")
stopifnot(file.exists("Test.tiff"))
file.remove("Test.tiff")