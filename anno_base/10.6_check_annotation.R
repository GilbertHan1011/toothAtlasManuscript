#== check
library(dplyr)
library(purrr)
library(ggplot2)
fileName = list.files("process/annotation/first_round_base/anno/",pattern = "*.csv",full.names = T)
fileNameShort <- fileName %>% gsub(".*?/(\\w+).csv","\\1",.)
files <- lapply(fileName,read.csv,stringsAsFactors = F)
# check whether column number of files equal 2
lapply(files,ncol)
annoBind <- do.call(rbind,files)
annoBind[,2] %>% unique()

# find which files have error 
names(files) <- fileNameShort

files <- map2(files, fileNameShort, function(x, y) {
  x[["project"]] <- y
  return(x)
})
annoBind <- do.call(rbind,files)
annoBind[annoBind$Coarse_Label_1 == "Mesenchme",]

annoBind$Coarse_Label_1 <- gsub("Mesenchme","Mesenchyme",annoBind$Coarse_Label_1)
annoBind$Coarse_Label_1 <- gsub("Perivasular","Perivasular",annoBind$Coarse_Label_1)
annoBind$Coarse_Label_1 <- gsub("Imuune","Immune",annoBind$Coarse_Label_1)
annoBind$Coarse_Label_1 <- gsub("Plasma","Immune",annoBind$Coarse_Label_1)
annoBind$Coarse_Label_1 <- gsub("Perivasular","Perivascular",annoBind$Coarse_Label_1)
annoBind$Coarse_Label_1 <- gsub("14","Unknown",annoBind$Coarse_Label_1)
annoBind$Coarse_Label_1 <- gsub("15","Unknown",annoBind$Coarse_Label_1)
annoBind$Coarse_Label_1 <- gsub("Unkown","Unknown",annoBind$Coarse_Label_1)
annoBind$Coarse_Label_1 <- gsub("muscle","Muscle",annoBind$Coarse_Label_1)
annoBind$Coarse_Label_1 <- gsub("muscle and RBC","Muscle",annoBind$Coarse_Label_1)
annoBind$Coarse_Label_1 <- gsub("Muscle and RBC","Muscle",annoBind$Coarse_Label_1)
annoBind$Coarse_Label_1 <- gsub("Mesenchme","Mesenchyme",annoBind$Coarse_Label_1)
annoBind$Coarse_Label_1 <- gsub("perivascular and neuron","Perivascular",annoBind$Coarse_Label_1)
annoBind$Coarse_Label_1 <- gsub("perivascular","Perivascular",annoBind$Coarse_Label_1)
annoBind$Coarse_Label_1 <- gsub("neuron","Neuron",annoBind$Coarse_Label_1)
#annoBind[annoBind$Coarse_Label_1 == "14",]

# plot percentage of annoBind
annoBind$Coarse_Label_1 %>% 
  table() %>% 
  prop.table() %>% 
  as.data.frame() %>% 
  arrange(desc(Freq)) %>% 
  ggplot(aes(x = reorder(., Freq), y = Freq)) +
  geom_bar(stat = "identity", fill = "skyblue") +theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
        axis.title.x = element_text(size = 16, face = "bold"),
        axis.title.y = element_text(size = 16)) +
  labs(x = "Cell type", y = "Percentage")
ggsave("process/annotation/first_round_base/plot_summary//20241007_percentage.pdf",height = 6,width = 8)

uniqueAnno <- annoBind[c("Coarse_Label_1","project")]%>%unique
uniqueAnno$Coarse_Label_1 %>% 
  table() %>% 
  as.data.frame() %>% 
  arrange(desc(Freq)) %>% 
  ggplot(aes(x = reorder(., Freq), y = Freq)) +
  geom_bar(stat = "identity", fill = "skyblue") +theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
        axis.title.x = element_text(size = 16, face = "bold"),
        axis.title.y = element_text(size = 16)) +
  labs(x = "Cell type", y = "Freq")+
  ggtitle("Frequency of cell type across project")
ggsave("process/annotation/first_round_base/plot_summary//20241007_frequency.pdf",height = 6,width = 8)

write.csv(annoBind,"process/annotation/first_round_base/anno/Anno_summary.csv",quote = F,row.names = F)

#anno_inspect <- read.csv("process/annotation/first_round_base/anno/Anno_summary.csv")
#head(anno_inspect)
