
featureAll <- read.csv("processed_data/attributeName/varName.csv",row.names = 1) %>% unlist
featureMes <- read.csv("processed_data/attributeName/varName_mes.csv",row.names = 1) %>% unlist
featureExclude <- setdiff(featureAll,featureMes)
write
