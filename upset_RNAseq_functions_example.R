
get_intersects = function(listINPUT,numberOfTopSets){
  # This function finds all combinations of size 1 to length(list)
  stopifnot(!is.null(listINPUT))
  stopifnot(is.list(listINPUT))
  stopifnot(length(listINPUT)>1)
  stopifnot(length(names(listINPUT))==length(listINPUT))
  if(missing(numberOfTopSets)){
    numberOfTopSets=10
  }
  stopifnot(is.numeric(numberOfTopSets))
  
  # First remove is.na() elements from the lists
  removeNAFromList = function(x){
    if(length(which(is.na(x)))){
      y=x[-which(is.na(x))];return(y)
    }else{
      return(x)
    }
  }
  listINPUT_filtered = lapply(listINPUT,FUN = removeNAFromList)
  
  outList=list() # initialize the final output list
  
  for(i in seq(length(listINPUT_filtered),1)){ # reverse order because genes get put into 1 group only & want to put them in the most inclusive group
    combos = combn(names(listINPUT_filtered),i)
    for(j in 1:ncol(combos)){
      listName=paste(combos[,j],collapse='_AND_')
      shared = Reduce(intersect,listINPUT_filtered[combos[,j]]) # get shared elements
      if(length(shared)>0){
        # we have some elements shared with this list
        shared.uniq=shared[which(! shared %in% unlist(outList))] # check that they are not in more inclusive lists
        if(length(shared.uniq)>0){
          outList[[listName]]=shared.uniq
          # print(cat('List',listName,'has',length(shared),'elements of which',length(shared.uniq),'are unique.'))
        }
      }
    }
  }
  
  if(length(outList)<numberOfTopSets){numberOfTopSets=length(outList)}
  names2include = names(unlist(lapply(outList,FUN=length))[order(unlist(lapply(outList,FUN=length)),decreasing = TRUE)])[1:numberOfTopSets] # get the ordered names
  
  outList.reduced = outList[names2include]
  
  return(outList.reduced)
}

tabMyList=function(listInput,tab){
  stopifnot(!missing(listInput))
  stopifnot(!missing(tab))
  stopifnot(length(listInput)>0)
  for(i in 1:length(listInput)){
    .list = listInput[i]
    l=length(.list[[1]])
    name=gsub("[\r\n]", "",names(.list))
    cat(paste('###',name,paste0('(N=',l,')')),'\n')
    .dataframe = data.frame(.list)
    colnames(.dataframe)=name
    print(knitr::kable(.dataframe),caption='')
  }
}  

# just like tab my list except returns a table
# could make this a toggle in tabMyList ...
tabMyList2table=function(listInput,tab,table2join=NULL,printTable=TRUE){
  require(dplyr)
  stopifnot(!missing(listInput))
  stopifnot(!missing(tab))
  stopifnot(length(listInput)>0)
  for(i in 1:length(listInput)){
    .list = listInput[i]
    l=length(.list[[1]])
    name=gsub("[\r\n]", "",names(.list))
    cat(paste('\n\n###',name,paste0('(N=',l,')')),' {.tabset}\n')
    cat(paste('\n\n####','TABLE','\n'))
    .dataframe = data.frame(.list)
    colnames(.dataframe)='ext_gene'
    if(!is.null(table2join)){
      .dataframe = dplyr::left_join(.dataframe,table2join,by='ext_gene') 
    }
    print(knitr::kable(.dataframe,caption='') %>% kable_styling(bootstrap_options = c("striped", "hover","condensed"),font_size = 12))
    
    print(cat("\n\n"))
    
    
    print(cat("\n\n"))
    
    print(cat(paste('\n\n####','ENRICHMENT','\n')))
    
    # Do enrichment testing
    GO.BP = enrichGO(gene=.dataframe$ext_gene,
                     OrgDb=config$org.db,
                     ont ="BP",
                     keyType = "SYMBOL",
                     pAdjustMethod = "BH"
    )
    if(nrow(GO.BP)>0){
      print(enrichplot::dotplot(GO.BP, showCategory=30) + ggtitle("DotPlot - GO:Biological Process"))
    }else{
      print("No GO Biological Processes Enriched\n")
    }
    
    print(cat(paste('\n\n####','ENRICHMENT TABLE','\n')))
    
    print(knitr::kable(data.frame(GO.BP)))
    
    cat("\n\n")
  }# Go to next set
}    



# Get the masterList
# masterList <- readRDS(file='masterList_in_vivo.Rds')
# names(masterList)
# masterList <- masterList[1:3]
# names(masterList)

# GETTING FROM DGE TSVs
# files2get <- list.files()[grep('d3.*DGE.tsv',list.files())]
# list2upset = NULL
# for(f in files2get){
#   x <- read.delim(f,sep="\t")
#   name <- gsub('.DGE.tsv','',f)
#   xFiltUp <- list(dplyr::filter(x,FDR<0.05 & logFC>0)$ext_gene)
#   names(xFiltUp) <- paste0(name,'_Up')
#   xFiltDown <- list(dplyr::filter(x,FDR<0.05 & logFC<0)$ext_gene)
#   names(xFiltDown) <- paste0(name,'_Down')
#   list2upset <- c(list2upset,xFiltUp,xFiltDown)
# }
# length(list2upset)
# names(list2upset)

# GETTING FROM masterList
# list2upset = NULL
# for(i in 1:length(masterList)){
#   name <- names(masterList)[i]
#   xFiltUp <- list(dplyr::filter(masterList[[i]],FDR<0.05 & logFC>0)$ext_gene)
#   names(xFiltUp) <- paste0(name,'_Up')
#   xFiltDown <- list(dplyr::filter(masterList[[i]],FDR<0.05 & logFC<0)$ext_gene)
#   names(xFiltDown) <- paste0(name,'_Down')
#   list2upset <- c(list2upset,xFiltUp,xFiltDown)
# }
# length(list2upset)
# names(list2upset)
# 
# lapply(list2upset,FUN=length)


# ```{r upset,fig.height=10,fig.width=10,eval=TRUE}
# 
# 
# library(UpSetR)
# 
# myOrder=rev(c(
#   names(list2upset)[grep('Up',names(list2upset))],
#   names(list2upset)[grep('Down',names(list2upset))]
# ))
# 
# UpSetR::upset(fromList(list2upset),
#               nsets=length(list2upset),
#               order.by = 'freq',
#               decreasing=FALSE,
#               nintersects = NA,
#               sets=myOrder,
#               sets.bar.color='blue',
#               number.angles = 0,
#               text.scale = 1,
#               main.bar.color='darkblue',
#               keep.order=FALSE
# )
# 
# UpSetR::upset(fromList(list2upset),
#               nsets=length(list2upset),
#               order.by = 'degree',
#               decreasing=FALSE,
#               nintersects = NA,
#               sets=myOrder,
#               sets.bar.color='blue',
#               number.angles = 0,
#               text.scale = 1,
#               main.bar.color='darkblue',
#               keep.order=TRUE
# )
# 
# ```

# ```{r tabListUpset,results='asis',fig.height=11,fig.width=7,eval=TRUE}
# 
# 
# cat('## Gene Sets {.tabset}\n')
# # UpSetR::upset(fromList(list2upset),nsets=50)
# .intersects = get_intersects(list2upset,numberOfTopSets = 100)
# 
# 
# # Filter for sets with >n genes
# n=1
# take <- which(as.numeric(unlist(lapply(.intersects,FUN=length)))>n)
# cat("Only showing tabs for intersects with greater than",n,"genes: N =",length(take),"sets\n")
# library(kableExtra)
# tabMyList2table(.intersects[take],'##')
# 
# 
# 
# cat('\n')
# 
# ```
