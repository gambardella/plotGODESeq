source('plotGODESeq.R')

# load results from WebGestalt
goenrich_data <- read.table("GO-example2.csv",sep="\t",fill=T,quote="\"",header=T)

# rename the columns to make them less weird and compatible with the GOPlot package
colnames(goenrich_data)[colnames(goenrich_data) %in% c("geneset","R","OverlapGene_UserID")] <- c("ID","Enrich","genes")

# remove commas from descriptions, because they suck
goenrich_data$description <- gsub(',',"",goenrich_data$description)

# Load results from DESeq2
deseq_data <- read.table("DESeq-example2.csv", sep=",",fill=T,header=T,row.names=1)


plotGODESeq(goenrich_data,
            deseq_data,
            maxFDR = 0.05,
            collapse = 1,
            color="l2fc",
            lowCol = "deepskyblue4",
            midCol = "#DDDDDD",
            highCol = "firebrick",
            extrawidth=3,
            centered=FALSE,
            fixed_ymax = 10,
            leghoffset=-0.5,
            legvoffset=1.5,
            label = "description",
            scale = 0.7,
            maxFDRLab = 1e-12,
            minZscoreLab = 2.5,
            wrap = 15
            )
