pkgname <- "MTTK"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
library('MTTK')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
base::assign(".old_wd", base::getwd(), pos = 'CheckExEnv')
cleanEx()
nameEx("MTTKExample")
### * MTTKExample

flush(stderr()); flush(stdout())

### Name: MTTKExample
### Title: Example 'MTTKExperiment'
### Aliases: MTTKExample
### Keywords: datasets

### ** Examples

data("MTTKExample")
MTTKExample
genomeData(MTTKExample)




cleanEx()
nameEx("aggregateByLink")
### * aggregateByLink

flush(stderr()); flush(stdout())

### Name: aggregateByLink
### Title: Aggregate Assays Along a Link Path
### Aliases: aggregateByLink

### ** Examples

x <- makeExampleMTTKExperiment()

aggregateByLink(x, path = "gene_to_genome")
aggregateByLink(x, path = c("gene_to_ko", "ko_to_module"))




cleanEx()
nameEx("aggregateToGenome")
### * aggregateToGenome

flush(stderr()); flush(stdout())

### Name: aggregateToGenome
### Title: Aggregate Assays to the Genome Level
### Aliases: aggregateToGenome

### ** Examples

x <- makeExampleMTTKExperiment()
aggregateToGenome(x)




cleanEx()
nameEx("makeExampleMTTKExperiment")
### * makeExampleMTTKExperiment

flush(stderr()); flush(stdout())

### Name: makeExampleMTTKExperiment
### Title: Create the Packaged Example 'MTTKExperiment'
### Aliases: makeExampleMTTKExperiment

### ** Examples

x <- makeExampleMTTKExperiment()
x
names(links(x))




cleanEx()
nameEx("summarizeActivity")
### * summarizeActivity

flush(stderr()); flush(stdout())

### Name: summarizeActivity
### Title: Summarize RNA-over-DNA Activity
### Aliases: summarizeActivity

### ** Examples

x <- makeExampleMTTKExperiment()

summarizeActivity(x, by = "feature")
summarizeActivity(x, by = "genome")




### * <FOOTER>
###
cleanEx()
options(digits = 7L)
base::cat("Time elapsed: ", proc.time() - base::get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')
