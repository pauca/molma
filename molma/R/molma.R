# naming: action + filetype + extra

unfactor<-function(f, type="n"){
  if(!is.factor(f)){
    return(f)
  }else{
    f <- levels(f)[f]
    if(type == "n"){
      return(as.numeric(f))
    }
    if(type == "c"){
      return(as.character(f))
    }
  }
}

readSmi <- function(file){
  d <- read.table(file,sep="\t",quote="\"", comment.char="~")
  colnames(d)<- c( "SMILES" , paste( "FIELD", 1:(ncol(d)-1),sep=""))
  return(as.data.frame(d))
}

readSdf <-function(file, addIndex = F){
  con    <- file(file, open = "r")
  fields <-  readSdfFieldNames(file)
  molCounter <- 0
  lineInMol <- 0
  struc <- ""
  
  if(addmolPositionOriginal){
    lineTemplate <- matrix(rep(NA,length(fields)+2),nrow=1)  
    colnames(lineTemplate) <- c("STRUCTURE",fields,"INDEX")
  }else{
    lineTemplate <- matrix(rep(NA,length(fields)+1),nrow=1)  
    colnames(lineTemplate) <- c("STRUCTURE",fields)
  }
  
  res <- list()
  toSave<- lineTemplate
  label <- ""
  while (length(line <- readLines(con, n = 1, warn = FALSE)) > 0) {
    if(grepl( "^(\\$\\$\\$\\$)", line )){
      molCounter <- molCounter + 1
      toSave[1,"STRUCTURE"] <- struc
      if(addmolPositionOriginal){
        toSave[1,"INDEX"] <- molCounter 
      }
      res[[length(res)+1]]<-toSave
      
      lineInMol <- 0
      struc <- ""
      
      toSave <- lineTemplate
    }else{
      lineInMol  <- lineInMol + 1
      if(lineInMol == 1 ){ readingStructure <- T }
      
      if(readingStructure){
        struc <- paste(struc ,  line , "\n",sep="")
      }
      
      if( grepl( "^(M\\s\\sEND)", line ) ){ readingStructure <- F }
      
      if( grepl( "^(>\\s\\s<)", line ) ){
        label = strsplit(line,"[<>]")[[1]][3]
      }else{
        if(label != ""){
          toSave[1,label]<- line
          label <- ""
        }
      }
    }
  }
  close(con)
  
  return( as.data.frame(do.call(rbind, res)))
}

isSdf <- function(file){
  return(grepl( "(\\.sdf)$", file ))
}

countSdfMols <- function(file){
  count <- 0
  con   <- file(file, open = "r")  
  while (length(line <- readLines(con, n = 1, warn = FALSE)) > 0) {
    if(grepl( "^(\\$\\$\\$\\$)", line )){
      count <- count + 1 
    }
  }
  return(count)
}

readSdfFieldNames <- function(inFile){
  con  <- file(inFile, open = "r")
  fields <- list()
  while (length(line <- readLines(con, n = 1, warn = FALSE)) > 0) {
    if(grepl( "^(>\\s\\s<)", line )){
      l <-  gsub("^(>\\s\\s<)","",line)
      l <-  gsub(">(\\s)*$","",l)
      fields[[length(fields)+1]] <- l
    }
  }
  close(con)
  unique(unlist(fields))
}

writeSdf<-function( data,file, addIndex = F){
  con  <- file(file, open = "w")  
  cols <- colnames(data)
  if(!addIndex){    
    if( "INDEX" %in% colnames(data)){
      data <- subset(data,select=-INDEX) 
    }
  }else{
    if( !("INDEX" %in% colnames(data))){
      warning("No Index Column Provided!")
    }
  }
  if( "STRUCTURE" %in% colnames(data)){
    data$STRUCTURE <- unfactor(data$STRUCTURE,"c")
  }else{
    warning("No Structure Column Provided!")
    data$STRUCTURE <- "\n NoName\n\n  0  0  0  0  0  0            999 V2000\nM  END\n"
  }
  cols <- cols[which(cols != "STRUCTURE")]
  for( i in 1:nrow(data)){    
    writeLines(data[i,"STRUCTURE"], con , sep = "", useBytes = FALSE)
    for(j in cols){
      writeLines( paste( ">  <",j,">\n",sep="" ), con , sep = "", useBytes = FALSE)
      writeLines( paste( data[i,j],"\n",sep="" ) , con , sep = "", useBytes = FALSE)
      writeLines(  "\n" , con , sep = "", useBytes = FALSE)      
    }
    writeLines(  "$$$$\n" , con , sep = "", useBytes = FALSE) 
    
  }
  close(con)
}

splitSdfSetRandom <- function( file, prob=c( "Train" = 0.8, "Test"=0.2), seed = NA ){
  if(is.numeric(seed)){    set.seed(seed)  }
  data <- readSdf(file)
  data$Set <- sample(names(prob), nrow(data),replace=T,prob=prob)
  res <- c()
  f <- gsub( "(\\.sdf)$","",file)
  for( n in names(prob)){
    ff <- paste(f,"_",n,".sdf",sep="")
    writeSdf( data[data$Set==n,], ff)
    res <- c(res,ff)
  }
  return( res) 
}


getCharge<-function(data){
  
  res <- sapply(unfactor(data$STRUCTURE,"c"), function(structure){
    ss <- unlist(strsplit(unlist(structure),"\\n"))
    indx <- sapply( ss,function(line){ 
      return( grepl( "^(M\\s\\sCHG)",line))
    })
    
    if( sum(indx)==0 ){
      return(NA)
    }else{
       line <- ss[indx]
       charge <- 0
       splline <- unlist(strsplit(line,"\\s\\s"))
       # get "M  CHG NR_CHARG_ATOMS ATOM_NR VALUE"
       for( i in 0:( as.numeric(splline[3])-1)){
         charge <- charge + as.numeric(splline[ 5 + 2*i ])
       }
       return(charge)
    }
  })
  names(res)<-1:length(res)
  return(res) 
  
  }
     
 
# from plyr 1.8
mapvalues <- function (x, from, to, warn_missing = TRUE) 
{
  if (length(from) != length(to)) {
    stop("`from` and `to` vectors are not the same length.")
  }
  if (!is.atomic(x)) {
    stop("`x` must be an atomic vector.")
  }
  if (is.factor(x)) {
    levels(x) <- mapvalues(levels(x), from, to)
    return(x)
  }
  mapidx <- match(x, from)
  mapidxNA <- is.na(mapidx)
  from_found <- sort(unique(mapidx))
  if (warn_missing && length(from_found) != length(from)) {
    message("The following `from` values were not present in `x`: ", 
            paste(from[!(1:length(from) %in% from_found)], collapse = ", "))
  }
  x[!mapidxNA] <- to[mapidx[!mapidxNA]]
  x
}

