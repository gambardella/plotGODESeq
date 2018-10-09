# plotGODESeq
R function to integrate differential gene expression and Gene Ontology over-representation in a 2D plot. Every aspect of the plot can be configured. The initial idea came from the function GOBubble in the package GOplot by Wencke Walter <http://wencke.github.io/>.

At the moment, there is only one function in the script, plotGODESeq(). It takes two mandatory inputs, a file containing Gene Ontology enrichment data, and a file containing differential gene expression data. A large variety of arguments can be used to customise the plot, but none are mandatory. E.g

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
            wrap = 10
```
