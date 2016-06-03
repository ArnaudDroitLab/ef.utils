library(eric.utils)



# Unit tests
guillimin.meta = import.into.grl(input.dir="/rap/fhq-091-aa/Working_Directory/Imene/GM12878/Guillimin", file.format="narrow", dir.type="mugqic")
guillimin.nometa = import.into.grl(input.dir="/rap/fhq-091-aa/Working_Directory/Imene/GM12878/Guillimin", file.format="narrow", dir.type="mugqic", discard.metadata=TRUE)

plain.narrow = import.into.grl(input.dir="/lustre2/rap/fhq-091-aa/Working_Directory/Imene/GM12878/Analyses/input", file.format="narrow")
plain.bed = import.into.grl(input.dir="/lustre2/rap/fhq-091-aa/Working_Directory/Imene/GM12878/Analyses/input", file.format="bed")

intersect.test  = build.intersect(plain.meta)