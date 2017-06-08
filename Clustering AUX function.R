cutdata <- tempdata[,1]
#-----------------------------------------------------------------------------------------------------------------------
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
          d_cut <- as.numeric(cutdata[which(info.gain == max(info.gain,na.rm = T))])[1]
          return(d_cut)
}



#-----------------------------------------------------------------------------------------------------------------------
#----------------------------function to return the best cut for a dimention and the 2 step ahead cut-------------------
#-----------------------------------------------------------------------------------------------------------------------

getCut <- function(data){
          #data <- cutdata
          if(length(data) < 20){
                    print("no cut")
          }else{
                    total.num<- length(data)
                    #for know dimension each value calc the proportion of data by split
                    Range <- max(data) - min(data)
                    #get best cut 1
                    dcut1 <- cut(data)
                    # right of cut1 region get next best cut 
                    cutdir <- density_dir(data,dcut1)
                    
                    if(cutdir == "Left"){
                              temp <- data[data < as.numeric(dcut1)]}else{
                                        temp <- data[data > as.numeric(dcut1)]
                              }
                    
                    dcut2 <- cut(temp)
                    # calc density measure
                    if(cutdir == "Left"){
                              
                              dcut1_dcut2_Y <-  length(temp[temp > dcut2]) 
                              dcut1_dcut2_N <-  total.num*(abs(dcut2-dcut1))/Range  
                              
                              desity_cut12 <- dcut1_dcut2_Y / dcut1_dcut2_N
                              
                              dcut2_r_Y <-  length(temp[temp < dcut2]) 
                              dcut2_r_N <-  total.num*(abs(dcut2 - min(data)))/Range  
                              
                              desity_cut2r <- dcut2_r_Y / dcut2_r_N
                              
                              if(desity_cut2r < desity_cut12){ 
                                        # if the region between cut1 and cut2 denser than that between cut2 and right boundary then cut2 is final
                                        bestcut <- dcut2
                                        final_N <- total.num*abs(bestcut - min(data))/Range
                                        final_Y <- length(data[data < bestcut])
                                        final.density <- final_Y/final_N} else{
                                                  # if the region between cut1 and cut2 less dense than that between cut2 and right boundary then perform cut3 and cut3 is final
                                                  temp1 <- temp[temp > dcut2]
                                                  bestcut <- cut(temp1)
                                                  final_N <- total.num*abs(bestcut - dcut1)/Range
                                                  final_Y <- length(data[data < dcut1 & data > bestcut])
                                                  final.density <- final_Y/final_N
                                        }
                    } else {
                              dcut1_dcut2_Y <-  length(temp[temp < dcut2]) 
                              dcut1_dcut2_N <-  total.num*(abs(dcut2-dcut1))/Range  
                              
                              desity_cut12 <- dcut1_dcut2_Y / dcut1_dcut2_N
                              
                              dcut2_r_Y <-  length(temp[temp > dcut2]) 
                              dcut2_r_N <-  total.num*(abs(dcut2 - max(data)))/Range  
                              
                              desity_cut2r <- dcut2_r_Y / dcut2_r_N
                              
                              if(desity_cut2r < desity_cut12){ 
                                        # if the region between cut1 and cut2 denser than that between cut2 and right boundary then cut2 is final
                                        bestcut <- dcut2
                                        final_N <- total.num*abs(bestcut - max(data))/Range
                                        final_Y <- length(data[data > bestcut])
                                        final.density <- final_Y/final_N
                              } else{
                                        # if the region between cut1 and cut2 less dense than that between cut2 and right boundary then perform cut3 and cut3 is final
                                        temp1 <- temp[temp < dcut2]
                                        bestcut <- cut(temp1)
                                        final_N <- total.num*abs(bestcut - dcut1)/Range
                                        final_Y <- length(data[data > dcut1 & data < bestcut])
                                        final.density <- final_Y/final_N
                              }
                              
                    }
                    
                    # store result in a list -- bestcut, lhsdata, rhsdata
                    final.cut <- data.frame(bestcut,final.density)
                    return(final.cut)
          }
}


#-----------------------------------------------------------------------------------------
#------------------------------ Density test make sure ahead looking cut only in low density regions
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

#-----------------------------------------------------------------
# get range
#-----------------------------------------------------------------

getRange <- function(data){
          temp <- range(data)
          temp <- temp[2] - temp[1]
          return(temp)
}

#-----------------------------------------------------------------
# get overall best cut
#-----------------------------------------------------------------

getFeatureCut <- function(data){
          
          features <- names(data)
          
          result <- list()
          N <- nrow(data)
          for (c in features){
                    
                    if (is.numeric(data[,c])){
                              featuredata <- data[,c]
                              result[[c]] <- getCut(featuredata)}
                    else{
                              next
                    }
          }
          result <- do.call(rbind,result)
          return(result)
}

cutResult <- getFeatureCut(test)

#------------------------------------------------------------------------------
# get the best cut from variable and apply to data, split data frame
#------------------------------------------------------------------------------

getBestCut<- function(data,cutResult){
          cutResult$var<- rownames(cutResult)
          out <- cutResult[cutResult$final.density == min(cutResult$final.density),]
          return(out)
} 

runCut(data)


runCut <- function(data){
          
          cutResult<- getFeatureCut(data)
          output <- getBestCut(data,cutResult)
          
          return(output)
}

#-------------------------------------------------------------------------------
#  function to grow tree
#-------------------------------------------------------------------------------

data <- test
depth <- 3
minobs <- 10
growClust <- function(data ,depth, minobs){
          numnodes <- 2^(depth + 1) - 1
          # initialize trace matrix, result table
          idx <- matrix(nrow = nrow(data),ncol = numnodes )
          rtable <- data.frame(var = NA,bestcut=NA)
          # root nodes
          idx[,1] <-  rep(x = 1,nrow(data))
          
          # first cut
          rtable$var[1] <- runCut(data)$var
          rtable$bestcut[1]<- runCut(data)$bestcut
          indexL <- 1*(data[,rtable[1,1]] < rtable[1,2]+ 0.00001 & as.logical(idx[,1]))
          indexR <- 1*((data[,rtable[1,1]] > rtable[1,2]) & as.logical(idx[,1]))
          idx[,2] <- indexL
          idx[,3] <- indexR
          
          # build out tree
          for(i in 2:(2^(depth)-1)){

                    
                    # add min obs stoping 
                    if (sum(idx[,i]) > minobs){
                             print(i)
                              #i = 7
                              tempdata <- data[as.logical(idx[,i]),]
                              rtable[i,1] <- runCut(tempdata)$var
                              rtable[i,2] <- runCut(tempdata)$bestcut
                              indexL <- 1*(data[,rtable[i,1]] < rtable[i,2] & as.logical(idx[,i]))
                              indexR <- 1*(!(data[,rtable[i,1]] < rtable[i,2]) & as.logical(idx[,i]))
                              idx[,2*i] <- indexL
                              idx[,2*i + 1] <- indexR          
                    }else{
                              rtable[i,1] <- "no cut"
                              rtable[i,2] <- NA
                              idx[,2*i] <- 0
                              idx[,2*i + 1] <- 0
                    }
                    
          }
          return(rtable)
}

ggplot(clusterdata, aes(x = X1 , y =X2 )) + 
          geom_point()+
          geom_vline(xintercept =  as.numeric(0.2263293),col ="blue") +
          geom_segment(x = -0.2 ,xend =0.2263293 ,y = 0.56587511,yend = 0.56587511 , col = "purple")+
          geom_segment(x = 0.2263293 ,xend =1.2263293 ,y = 0.55646372,yend = 0.55646372 , col = "purple")+
          geom_segment(x = -0.2 ,xend =0.2263293 ,y = 0.53092844,yend = 0.53092844 , col = "purple") +
          geom_segment(x = -0.2 ,xend =0.2263293 ,y = 0.68425761,yend = 0.68425761 , col = "purple") +
          geom_segment(x = 0.2263293 ,xend =1.2263293 ,y = 0.50764796,yend = 0.50764796 , col = "purple") +
          geom_segment(x = 0.55860068 ,xend =0.55860068 ,y = 0.55646372,yend = 1 , col = "purple") +
          
          geom_segment(x = 0.2779978 ,xend =0.2779978 ,y = 0,yend = 0.5149316 , col = "purple")+
          geom_segment(x = 0.3518865 ,xend =0.3518865 ,y = 0.5149316,yend = 1 , col = "purple")+
          geom_vline(xintercept =  as.numeric(0.7411899),col ="blue")

geom_vline(xintercept =  as.numeric(0.1582618),col ="blue") +
          geom_vline(xintercept =  as.numeric(-0.03242239),col ="blue") +
          geom_segment(x = 0.4853177,xend =1.1 ,y = 0.4658396,yend = 0.4658396 , col = "purple") + 
          geom_segment(x = -0.03242239 ,xend =0.1821867 ,y = 0.5149316,yend = 0.5149316 , col = "purple")




