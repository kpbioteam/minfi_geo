require(GEOquery, quietly = TRUE)
require(minfi, quietly = TRUE)
require(IlluminaHumanMethylation450kanno.ilmn12.hg19, quietly = TRUE)

options(warn = -1)

args <- commandArgs(trailingOnly = TRUE)
input1 = args[1] 
input2 = args[2]
input3 = args[3]
output = args[4] 

gset <- getGEO(input1)

 if(length(gset)==0) stop("Empty list retrieved from GEO.")
    if(length(gset)>1){
        warning("More than one ExpressionSet found:\n",names(gset),"\nUsing entry ",input2)
        gset <- gset[[input2]]
    } else gset <- gset[[1]]

platform <- annotation(gset)

gr <-  getLocations(IlluminaHumanMethylation450kanno.ilmn12.hg19,orderByLocation = TRUE)

locusNames <- names(gr)

sampleNames(gset) <- gset$title

common <- intersect(locusNames, featureNames(gset))
if(length(common)==0)
  stop("No rowname matches. 'rownames' need to match IlluminaHumanMethylation450k probe names.")

ind1 <- match(common,fData(gset)$Name)
ind2 <- match(common,locusNames)

preprocessing <- c(rg.norm=paste0('See GEO ',input1,' for details'))

what <- input3
if(what=="Beta"){
out <- GenomicRatioSet(gr=gr[ind2,],
                       Beta=gset[ind1,,drop=FALSE],
                       annotation=c(array = "IlluminaHumanMethylation450k"),
                       preprocessMethod=preprocessing)
} else {
  
out <- GenomicRatioSet(gr=gr[ind2,],
                         M=mat[ind1,,drop=FALSE],
                       annotation=c(array = "IlluminaHumanMethylation450k"),
                       preprocessMethod=preprocessing)
}
save(out, file = output)

  
