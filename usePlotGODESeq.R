source('plotGODESeq.R')

# Load results from DESeq2
deseq_data <- read.table("DESeq-example.csv", sep=",",fill=T,header=T,row.names=1)

# load results from WebGestalt
goenrich_data <- read.table("GO-example.csv", sep="\t",fill=T,quote="\"",header=T)

# rename the columns to make them less weird and compatible with the GOPlot package
colnames(goenrich_data)[colnames(goenrich_data) %in% c("geneset","R","OverlapGene_UserID")] <- c("ID","Enrich","genes")

# remove commas from descriptions, because they suck
goenrich_data$description <- gsub(',',"",goenrich_data$description)

plotGODESeq(goenrich_data,
            deseq_data,
            maxFDR = 1e-8,
            collapse = 0.9,
            color="l2fc",
            lowCol = "darkcyan",
            midCol = "yellow",
            highCol = "darkred",
            extrawidth=1,
            centered=0.5,
            leghoffset=-0.5,
            legvoffset=1,
            label = "description",
            scale = 0.8,
            maxFDRLab = 1e-12,
            minZscoreLab = 2,
            wrap = 15
            )
