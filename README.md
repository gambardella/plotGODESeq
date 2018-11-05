# plotGODESeq
R function to integrate differential gene expression and Gene Ontology over-representation in a 2D plot. Every aspect of the plot can be configured. The initial idea came from the function GOBubble in the package GOplot by Wencke Walter <http://wencke.github.io/>.

At the moment, there is only one function in the script, plotGODESeq(). It takes two mandatory inputs, a file containing Gene Ontology enrichment data, and a file containing differential gene expression data. Note that the function works better if the dataset is limited, in particular the number of GO terms. It is useful to analyse the effect of a perturbation, chemical or genetic, or to compare two cell types that are not too dissimilar. Comparing samples that exhibit several thousands of differentially expressed genes, resulting in thousands of enriched GO terms, will not only slow the function to a halt, it is also useless (GO enrichment should not be used in these conditions anyway).

A large variety of other arguments can be used to customise the plot, but none are mandatory. 

```R
plotGODESeq(goenrich_data, 
            deseq_data, 
            maxFDR = 0.05, 
            collapse = 1, 
            color="l2fc", 
            lowCol = "blue",
            midCol = "yellow",
            highCol = "#FF0000",
            extrawidth=1,
            centered=FALSE,
            leghoffset=-0.5,
            legvoffset=1,
            label = "description",
            scale = 4,
            maxFDRLab = 0.0001,
            minZscoreLab = 1,
            wrap = 10)
```

## Input

The Gene Ontology enrichment data must be a dataframe containing at least the columns: `ID` - the identifiers of the GO term, `description`- the description of the term, `Enrich` - the ratio of observed genes over expected genes annotated by the GO term, `FDR` - the False Discovery Rate (a.k.a. adjusted p-value), computed e.g. with the Benjamini-Hochberg correction, and `genes` - the list of observed genes annotated by the GO term, as in:

| ID        | description | Enrich | FDR | genes |
| --------- | ------------- | ----- | --- | --- |
| GO:0000123| interesting process | 1.58 | 1e-05 | DRD2;ADORA2B;PENK | 

Any other column can be present. It will not be taken into account. The order of columns does not matter.

The differential expression data must be a dataframe which rownames are the gene symbols, from the same namespace as the `genes` column of the GO enrichment data above. In addition, one column must be named `log2FoldChange`, containing the quantitative difference of expression between two conditions. Any other column can be present. It will not be taken into account. The order of columns does not matter.

## Configuration arguments

### *maxFDR* 
Only Gene Ontology terms with False Discovery Rate (a.k. adjusted p-values) smaller than this value will be considered for the plot. **Default is "0.05"**.

### *collapse*
If two many terms are plotted, making the plot unreadable, one can use the *reduce_overlap* function of the package GOplot to merge GO terms together when they share a certain fraction of annotated genes. Values can be either "FALSE", disabling the collpasing, or a value between 0 and 1, representing the proportion of genes in common for two GO terms to be merged. Beware, collapse = FALSE and collapse = "1" provide different results. The former does not call *reduce_overlap*, and all GO terms are plotted, while the latter will merge all terms with identical annotated gene sets. **Default is FALSE**.  

### *color*
Color scheme used to color the bubbles. Values can be "zscore" or "l2fc". The zscore is defined as ( nb(genes up) - nb(genes down) )/sqrt( nb(genes up) - nb(genes down) ). "l2fc" uses the mean of the log2(fold change) of all the genes annotated by the GO term. **Default is "zscore"**.

### *lowCol*
Color of the bubbles for the low values of the color scheme (left positions). **Default is "darkcyan"**. 

### *midCol*
Color of the bubbles for the meddium values of the color scheme (middle positions). **Default is "yellow"**. 

### *highCol*
Color of the bubbles for the high values of the color scheme (right positions). **Default is "darkred"**. 

### *scale*
Scale the bubbles radius. E.g. a value of "2" will increase the radius twofold while a value of "0.5" with decrease the radius twofold. **Default is "0.5"**. 

### *label*
Which information is used to label the bubbles. Values can be "id" for the bubbles to be labeled with the Identifier of the GO term (e.g. GO:0000001) or "description" for the bubbles to be labeled with the full description of the GO term (e.g. "mitochondrion inheritance"). **Default is "id"**. 

### *maxFDRLab*
Maximum Padj value at which bubbles are labelled. Bubbles with FDR or Padj above this value, but under maxFDR, will be plotted but not labeled. **Default is the value given to *maxFDR***.

### *minZscoreLab*
Minimum absolute zscore value at which bubbles are labeled. I.e. only bubbles with a zscore under -minZscoreLab or over minZscoreLab will be labeled. **Default is "0"**.

### *extrawidth*
Width added to the horizontal axis, on the left of min(zscore) and right of max(zscore). If the "centered" is set to TRUE, it is added on both side to max(abs(min(zscore)),abs(max(zscore))). **Default is "1"**.

### *centered*
If set to TRUE, the X axis will be centered on a zscore of 0. If FALSE, the X axis is chosen to balance the distribution of bubbles. **Default is FALSE**.  

### *leghoffset*
Value added to the horizontal position of the legend. **Default is "0"**. 

### *legvoffset*
Value added to the vertical position of the legend. **Default is "0"**. 

### *wrap*
Width of label's lines in characters. **Default is "15"**.