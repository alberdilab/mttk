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
