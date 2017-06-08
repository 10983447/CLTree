testdata <- read.csv("K:/R&D/Personal Lines/External Data/ABS/ABS 2011 Census/Postal Data/Census Processed Data Postcode.csv")

names(testdata)

test <- testdata[,c("c_Aus_Birthplace","c1_Popn_Priv_Dwellings","c_T_Own_with_Mortgage","c_Occ_Manager_Professional")]

testdata <- testdata[testdata$c_Aus_Birthplace > 0.61684 &
                       testdata$c1_Popn_Priv_Dwellings < 19174.00000,
             ]


test[is.na(test)] <- 0

summary(test)

growClust(test,depth = 3,minobs = 5)


names(test)

ggplot(testdata, aes(x = c_Aus_Birthplace , y =c_Occ_Manager_Professional , col = Region_Classification)) + 
          geom_point(alpha = 0.6) +
          geom_hline(yintercept =  as.numeric(19174.00000),col ="blue")  +
          geom_segment(x = 0.61684, xend =  as.numeric(0.61684), y =0, yend = 19174, col ="blue") 