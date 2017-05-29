
#-------------asdfasfsdf----------------------------------------------------------------------------------------------------------
#----------------------------function to perform cut for each data point of a dimention---------------------------------
#-----------------------------------------------------------------------------------------------------------------------
cut <- function(cutdata){
          Range<- max(cutdata) - min(cutdata)
          N <- length(cutdata)
          Y <- length(cutdata)
          D <- Y + N
          # first split
          info.gain <- as.numeric()
          for (i in 1:length(cutdata)){
                    # initialize
                    left_Y <- 0 
                    left_N <- 0
                    right_Y <- 0
                    right_N <- 0
                    
                    x <- as.numeric(cutdata[i])
                    #update for root
                    left_Y <- sum(cutdata < x) 
                    left_N <- (x -min(cutdata)) /Range * N
                    right_Y <- length(cutdata) - left_Y 
                    right_N <- N - left_N
                    
                    if(right_Y > right_N){right_N <- right_Y}
                    if(left_Y > left_N){left_N <- left_Y}
                    
                    
                    # create variables for info calc
                    dl <- left_Y + left_N
                    dr <- right_Y + right_N
                    ydl <- left_Y
                    ydr <- right_Y
                    ndl <- left_N
                    ndr <- right_N
                    # calc information gain -- info(x)D
                    
                    l.info <- - dl/D *(- ydl/dl * log2(ydl/dl) - ndl/dl * log2(ndl/dl))
                    r.info <- - dr/D *(- ydr/dr * log2(ydr/dr) - ndr/dr * log2(ndr/dr))
                    
                    infox <- l.info + r.info
                    info.gain[i] <- 1 - infox 
          }
          # cut 1
          d_cut <- as.numeric(cutdata[which(info.gain == max(info.gain,na.rm = T))])
          return(d_cut)
}


#-----------------------------------------------------------------------------------------------------------------------
#----------------------------function to return the best cut for a dimention---------------------------------
#-----------------------------------------------------------------------------------------------------------------------

getCut <- function(data,feature){
          #t_data <- clusterdata[,"X2"]
          t_data <- data[,feature]
          total.num<- length(t_data)
          #for know dimension each value calc the proportion of data by split
          Range <- max(t_data) - min(t_data)
          #get best cut 1
          dcut1 <- cut(t_data)
          # right of cut1 region get next best cut 
          cutdir <- density_dir(t_data,dcut1)
          
          if(cutdir == "Left"){
                    temp <- t_data[t_data < as.numeric(dcut1)]}else{
                              temp <- t_data[t_data > as.numeric(dcut1)]
                    }
          
          dcut2 <- cut(temp)
          # calc density measure
          
          dcut1_dcut2_Y <-  length(temp[temp < dcut2]) 
          dcut1_dcut2_N <-  total.num*(dcut2-dcut1)/Range  
          
          desity_cut12 <- dcut1_dcut2_Y / dcut1_dcut2_N
          
          dcut2_r_Y <-  length(temp[temp > dcut2]) 
          dcut2_r_N <-  total.num*(max(t_data)-dcut2)/Range  
          
          desity_cut2r <- dcut2_r_Y / dcut2_r_N
          
          if(desity_cut2r < desity_cut12){ 
                    # if the region between cut1 and cut2 denser than that between cut2 and right boundary then cut2 is final
                    cut <- dcut2} else{
                              # if the region between cut1 and cut2 less dense than that between cut2 and right boundary then perform cut3 and cut3 is final
                              temp1 <- temp[temp < dcut2]
                              cut <- cut(temp1)
                    }
         return(cut)
          
}

#------------------------------------------------------------------------------
# get the best cut from variable
#------------------------------------------------------------------------------




#-----------------------------------------------------------------------------------------
#------------------------------ Density test
#-----------------------------------------------------------------------------------------

density_dir <- function(data,numcut){
         Yr <- sum(numcut < data) 
         Range <- max(data) - min(data)
         total.num <- length(data)
         Nr <- total.num*(max(data) - numcut)/Range
         rdensity <- Yr/Nr
         Yl <- sum(numcut > data) 
         Nl <- total.num*(numcut - min(data))/Range
         ldensity <- Yl/Nl
         if(rdensity > ldensity){
                   direction <- "Left"
         }else{direction <- "Right"}
         return(direction)
}








clusterdata <- data.frame(A$X)
t_data <- clusterdata[,"X2"]
names(clusterdata) <- c("X1","X2")
x <- getCut(clusterdata,"X1")
y <- getCut(clusterdata,"X2")



ggplot(clusterdata, aes(x = X1 , y =X2 )) + 
          geom_point()+
          geom_vline(xintercept =  as.numeric(x),col ="blue") +
          geom_hline(yintercept =  as.numeric(y),col = "purple")


rm(dcut1)
