 

drawSmile<-function( smi, width = 200, height = 200){
  require(rcdk)
  require(ggplot2)
  require(reshape2)
  #smile<- "c1ccccc1C(Cl)(Br)c1ccccc1"
  mol <- parse.smiles(smi)
  atoms <- get.atoms(mol[[1]])
  get.point2d(atoms[[1]])
  
  
  img <- view.image.2d(mol[[1]],width, height )
  
  d1 <- melt(img[,,1]) 
  d2 <- melt(img[,,2]) 
  d3 <- melt(img[,,3]) 
  
  colnames(d1)[3]<- "r"
  colnames(d2)[3]<- "g"
  colnames(d3)[3]<- "b"
  
  dd <- cbind(cbind(d1,d2$g),d3$b)

  dd$hex <- apply( dd , 1, function(l){return(rgb(l[3],l[4],l[5]))})
  
  g <- ggplot(dd)+
    theme_bw()+
    geom_raster(aes(Var1, Var2, fill=hex))+coord_equal()+
    xlab("")+ylab("")+scale_fill_identity()+
    theme(axis.ticks = element_blank(), axis.text = element_blank())
  
  return(g)
}