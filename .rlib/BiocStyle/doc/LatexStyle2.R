### R code from vignette source 'LatexStyle2.Rnw'

###################################################
### code chunk number 1: style
###################################################
BiocStyle::latex()


###################################################
### code chunk number 2: fig
###################################################
par(mar=c(4,4,0.5,0.5))
plot(cars)


###################################################
### code chunk number 3: widefig
###################################################
par(mar=c(4,4,0.5,0.5))
plot(cars)


###################################################
### code chunk number 4: smallfig
###################################################
par(mar=c(4,4,0.5,0.5))
plot(cars)


###################################################
### code chunk number 5: squarefig
###################################################
par(mar=c(4,4,0.5,0.5))
plot(cars)


###################################################
### code chunk number 6: figureexample
###################################################
par(mar=c(4,4,0.5,0.5))
v = seq(0, 60i, length=1000)
plot(abs(v)*exp(v), type="l", col="Royalblue")


###################################################
### code chunk number 7: sessionInfo
###################################################
toLatex(sessionInfo())


###################################################
### code chunk number 8: latexPkgs
###################################################
bioconductor.sty = readLines(BiocStyle:::bioconductor.sty)
pattern = "^\\\\RequirePackage(\\[.*\\])?\\{([[:alnum:]]+)\\}.*"
pkgs = sub(pattern, "\\2", grep(pattern, bioconductor.sty, value = TRUE))
# add hyperref which is not captured by the regex
pkgs = c(pkgs, "hyperref")
# list sorted and unique names
pkgs = sort(unique(pkgs))
latexpkgs = paste0("\\\\latex{", pkgs, "}", collapse=", ")


