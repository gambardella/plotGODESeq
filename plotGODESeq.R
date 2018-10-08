# plotGODESeq
# Copyright (C) 2018 Nicolas Le Novère
# 
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with this program. If not, see <http://www.gnu.org/licenses/>.

# This function aims at integrating Gene Ontology overepresentation analysis and differential gene expression.

plotGODESeq <- function(goenrich,
                        deseq, 
                        maxFDR,
                        collapse,
                        color,
                        scale,
                        label,
                        maxFDRLab,
                        minZscoreLab,
                        extrawidth,
                        centered,
                        leghoffset,
                        legvoffset,
                        wrap){
  
###### Input:

###    
# label: label of the bubbles. Can take values of "id" (default) or "description"
# char
# NB: Needs to be checked before goenrich, because value is used in checking goenrich
  if( missing(label) ){ 
    label = "id"
    cat("No label scheme was specified. GO term IDs will be used.\n")
  }
  if( !( label %in% c("id","description") ) ){ 
    label = "id"
    cat("The specified label scheme was neither \"id\" nor \"description\". GO IDs will be used.\n")
  }   
  
###  
# goenrich: result of over-representation analysis. The format is based on 
# dataframe(
# no rownames
# ID          : GO term ID (Factor)
# description : GO term description (character. Hopefully without ",")
# Enrich      : ratio between observed and expected annotated genes (numeric)
# PValue      : significance of the over-representation, from hypergeometric test (numeric) NOT STRICTLY NEEDED AT THE MOMENT
# FDR         : False Discovery Rate, a.k.a adjusted p-value, normally obtained using Benjamini-Hochberg correction (numeric)
# genes       : list of genes annotated by the GO term, separated with semicolumns (Factor)
# )
# What to do if no ORA input?
  if( missing(goenrich) ){
    stop("No Gene Ontology enrichment dataset provided.")
  }
  if( label == "id" & !( "ID" %in% colnames(goenrich) ) ){
    stop("Gene Ontology enrichment dataset does not provide term IDs, necessary to label the bubbles.")
  }
  if( label == "description" & !( "description" %in% colnames(goenrich) ) ){
      stop("Gene Ontology enrichment dataset does not provide term descriptions, necessary to label the bubbles.")
  }
  if( !( "Enrich" %in% colnames(goenrich) ) ){
    stop("Gene Ontology enrichment dataset does not provide a column \"Enrich\" with term enrichments, necessary for bubble size.")
  }
# THE FOLLOWING CODE IS NOT NEEDED AT THE MOMENT. KEPT FOR FUTURE USE  
#  if( !( "PValue" %in% colnames(goenrich) ) ){
#    stop("Gene Ontology enrichment dataset does not provide non corrected significance of term enrichment, necessary for ???.")
#  }
# FDR
  if( !( "FDR" %in% colnames(goenrich) ) ){
    stop("Gene Ontology enrichment dataset does not provide a column \"FDR\" with False Discovery Rates (a.k.a. adjusted P-values). They are necessary for plot Y axis.")
  }
# genes  
  if( !( "genes" %in% colnames(goenrich) ) ){
    stop("Gene Ontology enrichment dataset does not provide a column \"genes\" with the list of genes annotated by each GO term. This is necessary for computing the zscores used as X axis, and the mean differential expression, used to color bubbles.")
  }

###    
# color: color of the bubbles. Can take values of "zscore" (default) or "l2fc"
# char
# NB: Needs to be checked before deseq, because its value is used when checking deseq
  if( missing(color) ){ 
    color = "zscore"
    cat("No color scheme was specified. Zscores will be used.\n")
  }
  if( !( color %in% c("zscore","l2fc") ) ){ 
    color = "zscore"
    cat("The specified color scheme was neither \"zscore\" nor \"l2fc\". Zscores will be used.\n")
  }  
  
###
# deseq: result of DESeq2 analysis. We only use two type of information
# dataframe(
# rownames: Gene symbols
# log2FoldChange: Extent of the differential expression for each gene
# )
# What to do if no DESeq input?  
  if( missing(deseq) ){
    stop("No differential expression dataset provided.")
  }
  # L2FC 
  if( color == "l2fc" & !( "log2FoldChange" %in% colnames(deseq) ) ){
    stop("Differential Expression dataset does not provide a column \"log2FoldChange\" necessary to compute the plot color.")
  }

###  
# maxFDR: maximum value of the ORA False Discovery Rate that we consider for the plot 
# numeric
  if( missing(maxFDR) ){ 
    maxFDR = 0.05 
    cat("No maximum FDR provided for plotting. Default of 0.05 will be used.\n")
  }  

###
# collapse: value between 0 and 1, representing the proportion of genes in common for two GO terms to be merged. 
  if( missing(collapse) ){ collapse = 1  
    cat("No maximum GO term collapsing level is provided. Default of 1 will be used, i.e. only terms with 100% gene overlap will be collapsed in one bubble.\n")
  } 
  
###
# scale: scale the radius the bubbles. The bigger, the smaller the radius
# numeric 
  if( missing(scale) ){ scale = 2 }
  
###
# maxFDRLab: minimum Padj value at which bubbles are labelled  
# numeric
  if( missing(maxFDRLab) ){ 
    maxFDRLab = 0.05
    cat(paste("No maximum FDR provided for labelling. Default of",maxFDRLab,"will be used.\n"))
    }

###
# minZscoreLab: minimum absolute zscore value at which bubbles are labelled  
# numeric
  if( missing(minZscoreLab) ){ 
    minZscoreLab = 0 
    cat(paste("No minimum zscore provided for labelling. Default of",minZscoreLab,"will be used.\n"))
    }
      
###
# extrawidth: width added to the left of min(zscore) and right of max(zscore)
# if the "centered" is set to TRUE, it is added on both side to max(abs(min(zscore)),abs(max(zscore)))
# numeric
if( missing(extrawidth) ){ extrawidth = 1 }  

###
# centered: decides if the plot is symmetrical  
# Boolean
if ( missing(centered) ){ centered = FALSE }

### 
# value added to the horizontal position of the legend 
# numeric  
  if( missing(leghoffset)){leghoffset = 0 }
###
# value added to the vertical position of the legend
# numeric   
if( missing(legvoffset)){ legvoffset = 0 }
  
###
# wrap: width of label's lines
# numeric
if ( missing(wrap) ){ wrap = 15 }  
  
  ##############################
  # loaging required packages
  ##############################
  
  library(GOplot)   # necessary for the function collapsing GO terms that overlap
  library(tidyr)    # necessary for the function separate_rows
  library(SDMTools) # necessary for the function legen.gradient
  
  ##############################
  ## Preparation of the data
  ##############################
  
  # Replace the p-values and FDR that are zero so that they can be plotted on log-scale.

  # PVALUE IS NOT USED AT THE MOMENT. IF USED IN THE FUTURE, MOVE THE TEST TO INPUT PARSING ABOVE
  if ("PValue" %in% colnames(goenrich)){
    # Get the minimal non-0 p-value
    minpval <- min(goenrich[goenrich$PValue>0,]$PValue)
    # replace the 0 p-values by a tenth of the minimal ones
    goenrich[goenrich$PValue==0,]$PValue <- minpval/10
  }

  # Get the minimal non-0 FDR
  minfdr <- min(goenrich[goenrich$FDR>0,]$FDR)
  # replace the 0 FDR by a tenth of the minimal ones
  goenrich[goenrich$FDR==0,]$FDR <- minfdr/10
  
  # Prepare a subset of the GO enrichment. NB always use this subset, even if it is 100% of the entire results.
  goenrich_subset <- subset(goenrich, goenrich$FDR < maxFDR)
  
  ### Compute zscores, i.e. relative over or underexpression of Genes annotated by each term.
  
  getZscore <- function(goterm) {
    # get the list off all genes annotated by a term
    genes<-strsplit(as.character(unlist(goterm)),";")
    
    # using the DESeq results, get the number of genes that are upregulated or downregulated
    up <- nrow(deseq_data[rownames(deseq_data) %in% unlist(genes)&deseq$log2FoldChange>=0,])
    down <- nrow(deseq_data[rownames(deseq_data) %in% unlist(genes)&deseq$log2FoldChange<0,])
    
    # compute the Zscore
    zscore <- (up - down)/sqrt(up+down)
  }
  
  # Get the Zscore for each GO term
  resZscore<-lapply(goenrich_subset$genes,getZscore)
  
  # add zscore to the GO enrichment table
  goenrich_subset$zscore <- unlist(resZscore)
  
### Compute average gene expression enrichment

  getMeanL2FC <- function(goterm) {
  # get the list off all genes annotated by a term
    genes<-unlist(strsplit(as.character(unlist(goterm)),";"))
    
    checkGeneDESeq <- function(gene){
      if ( !( gene %in% rownames(deseq_data) ) ){
        cat(paste("gene",gene,"is in the GO enrichment table but not in the list of differential expression\n"))
      }
    }
    lapply(genes,checkGeneDESeq)
    
    # compute the mean of L2FC for all the genes annotated by the term
    meanL2FC <- mean(deseq_data[rownames(deseq_data) %in% genes,]$log2FoldChange)
  }
  
  # Get the mean L2FC for each GO term 
  resL2FC<-lapply(goenrich_subset$genes,getMeanL2FC)
  
  # add average gene expression enrichment to the GO enrichment table
  goenrich_subset$meanL2FC <- unlist(resL2FC)

  # rename the FDR into adj_pval and description into term to be compatible with the GOPlot package
  colnames(goenrich_subset)[colnames(goenrich_subset) %in% c("description","FDR")] <- c("term","adj_pval")
  
  # explodes each GO term row into as many rows as genes involved
  # This is necessary for collapsing terms using the package GOPlot
  # separate_rows is a function of the package tidyr
  goenrich_expand <- separate_rows(goenrich_subset, genes)  
  print(paste("Nb of GO terms in input file): ",nrow(goenrich)))
  print(paste("Nb of GO terms after FDR threshold): ",nrow(goenrich_subset)))

  ##############################
  ## Plot the data
  ##############################
  
  ## Use GOPlot package to merge together terms annotating similar genesets
  enrich_red <- reduce_overlap(goenrich_expand, overlap = collapse)
  print(paste("Nb of GO terms after collapsing): ",nrow(enrich_red)))
  
  ####### Prepare colours
  
  # build the low palette
  rcPallow <- colorRampPalette(c("darkcyan","yellow"))
  # build the high palette
  rcPalhigh <- colorRampPalette(c("yellow","darkred"))
  
  ####### Colors based on zscores
  if (color == "zscore"){
    print(paste("Chosen color for the bubble is: ",color))
    
    # Extract zscores under 0
    downreg <- enrich_red[enrich_red$zscore < 0,]$zscore
    
    # Extract zscores above 0
    upreg <- enrich_red[enrich_red$zscore >= 0,]$zscore
    
    # Compute the number of breaks
    maxzscore <- max(abs(downreg),abs(upreg))
    
    # I want 10 colors max on each side
    ndownbreak = ceiling(10* max(abs(downreg))/maxzscore)
    nupbreak = ceiling(10* max(abs(upreg))/maxzscore  )
    
    # Build color arrays using the zscore. Use ceiling and floor to ensure the right number in total
    datcollow <- rcPallow(ndownbreak)[as.numeric(cut(downreg,breaks = ndownbreak))]
    datcolhigh <- rcPalhigh(nupbreak)[as.numeric(cut(upreg,breaks = nupbreak))]
    
    # build df to contain all zscores and colors
    datcol <- data.frame(zscore = enrich_red$zscore,colour=NA)
    
    # Populate the colours
    datcol[datcol$zscore < 0,]$colour <- datcollow
    datcol[datcol$zscore >= 0,]$colour <- datcolhigh
  }
  
  ####### Colors based on L2FC
  if (color == "l2fc"){
    print(paste("Chosen color for the bubble is: ",color))

    # Extract L2FC under 0
    downreg <- enrich_red[enrich_red$meanL2FC < 0,]$meanL2FC
    
    # Extract L2FC above 0
    upreg <- enrich_red[enrich_red$meanL2FC >= 0,]$meanL2FC
    
    # Compute the number of breaks
    maxl2fc <- max(abs(downreg),abs(upreg))
    
    # I want 10 colors max on each side
    ndownbreak = ceiling(20* max(abs(downreg))/maxl2fc)
    nupbreak = ceiling(20* max(abs(upreg))/maxl2fc  )
    
    # Build color arrays using the L2FC. Use ceiling and floor to ensure the right number in total
    datcollow <- rcPallow(ndownbreak)[as.numeric(cut(downreg,breaks = ndownbreak))]
    datcolhigh <- rcPalhigh(nupbreak)[as.numeric(cut(upreg,breaks = nupbreak))]
    
    # build df to contain all L2FC and colors
    datcol <- data.frame(L2FC = enrich_red$meanL2FC,colour=NA)
    
    # Populate the colours
    datcol[datcol$L2FC < 0,]$colour <- datcollow
    datcol[datcol$L2FC >= 0,]$colour <- datcolhigh
  }
  
  ####### create width and height of plot

  if(centered == TRUE){
    xmin = -max(abs(min(enrich_red$zscore)),abs(max(enrich_red$zscore))) - extrawidth
    xmax = max(abs(min(enrich_red$zscore)),abs(max(enrich_red$zscore))) + extrawidth
      } else {
  xmin = min(enrich_red$zscore) - extrawidth
  xmax = max(enrich_red$zscore) + extrawidth
}
  ymin = min(-log10(enrich_red$adj_pval)) - 0.5
  ymax = max(-log10(enrich_red$adj_pval)) + 0.5
  
  ###### Plot!
  
  symbols(enrich_red$zscore, 
          -log10(enrich_red$adj_pval),
          circle = sqrt(enrich_red$Enrich/pi)/scale, # gives the ratius in fonction of enrichment. 
          inches=FALSE,
          fg="white",
          bg = datcol$colour,
          xlim=c(xmin,xmax),ylim=c(ymin,ymax),
          xlab="zscore",ylab="-log(adjusted p-value GO enrichment)"
  )
  ###### Create the Legend 
  
  # Create a color ramp
  rbPal <- colorRampPalette(c("darkred","yellow","darkcyan"))
  legcol <- rbPal(50)
  
  #points for the gradient legend
  
  pnts = cbind(x =c(xmin+0.5+leghoffset,xmin+0.75+leghoffset,xmin+0.75+leghoffset,xmin+0.5+leghoffset), 
               y =c(ymin+2.5+legvoffset,ymin+2.5+legvoffset,ymin+1+legvoffset,ymin+1+legvoffset))
  
  #create the gradient legend
  
  legend.gradient(pnts,
                  legcol,
                  c(sprintf("%.2f",min(enrich_red$meanL2FC)),sprintf("%.2f",max(enrich_red$meanL2FC))), 
                  title = "mean DESeq2 L2FC")
  
  # Legend for enrichment: black circle of surface "1". 
  # NB: I have to feed a vector of 1 element for the size to work. Don't ask. Symbols is crazy
  symbols(c(xmin+0.5+leghoffset+0.125), # move by 0.125 because gradient is 0.25 large
          c(ymin+3.5+legvoffset),
          circle = c(sqrt(1/pi))/scale, # gives the ratius in fonction of enrichment. 
          inches=FALSE,
          fg="white",
          bg = "black",
          add = TRUE
  )
  text(xmin+0.5+leghoffset+0.125+sqrt(1/pi)/2,ymin+3.5+legvoffset,label="Enrich = 1",pos=4)
    
  # Select the bubbles that will be labeled. 
  # We only select adj p-value less than maxFDRLab and zscore inferieur to -minZscoreLab or superieur to minZscoreLab
  tolabel <- subset(enrich_red, enrich_red$adj_pval < maxFDRLab | abs(enrich_red$zscore)>minZscoreLab)
  
  
  if (label == "id"){
    ####### GO ID as labels

    text(enrich_red$zscore, 
         -log10(enrich_red$adj_pval),
         enrich_red$ID,
         cex=0.75
    )
  } 
  
  if (label == "description"){
  ####### GO description as labels
  
  ## Build wrapped label
  
  # Core wrapping function
    wrap.it <- function(x, len)
    { 
      sapply(x, function(y) paste(strwrap(y, len), 
                                  collapse = "\n"), 
             USE.NAMES = FALSE)
    }
  
    # Call this function with a list or vector
    wrap.labels <- function(x, len)
    {
      if (is.list(x))
      {
        lapply(x, wrap.it, len)
      } else {
        wrap.it(x, len)
      }
    }
  
    labels <- wrap.labels(tolabel$term, wrap)
  
    text(tolabel$zscore, 
         -log10(tolabel$adj_pval),
         labels,
         cex=0.75
    )
  }
  
} ### END OF MAIN FUNCTION