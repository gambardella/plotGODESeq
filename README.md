# plotGODESeq
R function to integrate differential gene expression and Gene Ontology over-representation in a 2D plot. Every aspect of the plot can be configured. The initial idea came from the function GOBubble in the package GOplot by Wencke Walter <http://wencke.github.io/>.

At the moment, there is only one function in the script, plotGODESeq(). It takes two mandatory inputs, a file containing Gene Ontology enrichment data, and a file containing differential gene expression data. A large variety of other arguments can be used to customise the plot, but none are mandatory. 

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

The Gene Ontology enrichment data must be dataframe containing at least the columns: `ID` - the identifiers of the GO term, `description`- the description of the term, `Enrich` - the ratio of observed genes over expected genes annotated by the GO term, `FDR` - the False Discovery Rate (a.k.a. adjusted p-value), computed e.g. with the Benjamini-Hochberg correction, and `genes` - the list of observed genes annotated by the GO term, as in:

| ID        | description | Enrich | FDR | genes |
| --------- | ------------- | ----- | --- | --- |
| GO:0000123| interesting process | 1.58 | 1e-05 | DRD2;ADORA2B;PENK | 

Any other column can be present. It will not be taken into account. The order of columns does not matter.

The differential expression data must be a dataframe which rownames are the gene symbols, from the same namespace as the `genes` column of the GO enrichment data above. In addition, one column must be named `log2FoldChange`, containing the quantitative difference of expression between two conditions. Any other column can be present. It will not be taken into account. The order of columns does not matter.

**WARNING:** The current code does not work if the color scheme chosen for the bubbles is based on variable, `l2fc` or `zscore`, that do not contain negative **and** positive values. 

## Configuration arguments
