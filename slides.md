---
title: "Some fussing about SummarizedExperiment"
author: Tim Triche, Jr.
date: July 20th, 2015
output: ioslides_presentation

---

# boilerplate

This is all very nice...

```R 
library(SummarizedExperiment)
nrows <- 200; ncols <- 6
counts <- matrix(runif(nrows * ncols, 1, 1e4), nrows)
rowRanges <- GRanges(rep(c("chr1", "chr2"), c(50, 150)),
                     IRanges(floor(runif(200, 1e5, 1e6)), width=100),
                     strand=sample(c("+", "-"), 200, TRUE),
                     feature_id=sprintf("ID%03d", 1:200))
colData <- DataFrame(Treatment=rep(c("ChIP", "Input"), 3),
                     row.names=LETTERS[1:6])
rse <- SummarizedExperiment(assays=SimpleList(counts=counts),
                            rowRanges=rowRanges, colData=colData)
```

---

# wat (part I) 

This is perhaps less nice:

```R
is(SummarizedExperiment(), "SummarizedExperiment")
```

---

# wat (part II)

OK so maybe the previous example was silly, what about this?

```R
is(SummarizedExperiment(assays=SimpleList(counts=counts),
                        rowRanges=rowRanges, colData=colData),
   "SummarizedExperiment")
```

---

# wat (part III)

Fine, what about the previous-release behavior?

```R
is(GenomicRanges:::SummarizedExperiment(counts, rowRanges, colData),
   "SummarizedExperiment")
```

---

# This is a little silly

[It's not consistent with BioC deprecation guidelines, either](http://www.bioconductor.org/developers/how-to/deprecation/http://www.bioconductor.org/developers/how-to/deprecation/)

```R 
se <- GenomicRanges:::SummarizedExperiment(counts, rowRanges, colData)
rowData(se)
```

---

# No worries, this is just for internal use

Unless anyone reads the [Nature Methods paper](http://www.nature.com/nmeth/journal/v12/n2/fig_tab/nmeth.3252_F2.html)

![Doh](images/natmeth.jpg)


