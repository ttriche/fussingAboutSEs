---
title: "Some fussing about SummarizedExperiment"
author: Tim Triche, Jr.
date: July 20th, 2015
output: ioslides_presentation

---

# A bedrock data structure

So hey, we have a complicated experiment and we'd like to summarize it. 

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

Said experiment involves a lot of flat priors.

---

# wat (part I) 

However, apparently that's not what the function does?

```R
is(SummarizedExperiment(), "SummarizedExperiment")
```
```
FALSE 
```

---

# wat (part II)

OK so maybe the previous example was silly, what about this?

```R
is(SummarizedExperiment(assays=SimpleList(counts=counts),
                        rowRanges=rowRanges, colData=colData),
   "SummarizedExperiment")
```
```
FALSE 
```

---

# wat (part III)

Fine, what about the previous-release behavior?

```R
is(GenomicRanges:::SummarizedExperiment(counts, rowRanges, colData),
   "SummarizedExperiment")
```
```
FALSE 
```
---

# No worries, this is just for internal use

Unless anyone reads the [Nature Methods paper](http://www.nature.com/nmeth/journal/v12/n2/fig_tab/nmeth.3252_F2.html)

![D'oh](images/natmeth.jpg)

---

# I'm so excited about this thing I saw in Nature Methods!!1

```R 
se <- GenomicRanges:::SummarizedExperiment(counts, rowRanges, colData)
rowData(se)
```
```
Error: 'rowData' is defunct.
Use 'rowRanges' instead.
See help("Defunct")
```

[That's not consistent with BioC deprecation guidelines](http://www.bioconductor.org/developers/how-to/deprecation/http://www.bioconductor.org/developers/how-to/deprecation/)
Speaking as a chump who actually followed the rules, this is a drag.
(You could always argue that WORKING protocols go in Nature Protocols, however.)

---

# Bonus!  For Nature Methods readers

Proving once again that anything published is already out of date:  

```R
exptData(se)
```
```
List of length 0
Warning message:
'exptData' is deprecated.
Use 'metadata' instead.
See help("Deprecated") 
```

But at least it's not defunct() yet like that other one. 
Maybe the thing to do here is issue a corrigendum?
That always makes everything better

---

# Fin 
