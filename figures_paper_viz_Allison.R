####make the working directory, upload data####
setwd("~/Project_Managment/01_PQuiescence/02_paper/04_Results/01_FlowCytometry_Data")

####Prepare Data of EGFP and Viability for Supp_Table2 and Fig2A####
data1 <- read_csv("01_Expe_EGFP_PAT_BPK/03_FCS_Analysis/2024_BPK_01_PAT_D1-D3.csv")
data2 <- read_csv("01_Expe_EGFP_PAT_BPK/03_FCS_Analysis/2024_BPK_02_PAT_D4-D8.csv")
data3 <- read_csv("01_Expe_EGFP_PAT_BPK/03_FCS_Analysis/2024_BPK_03_PAT_D10.csv")
data4 <- read_csv("01_Expe_EGFP_PAT_BPK/03_FCS_Analysis/2024_BPK_04_PAT_D1-D4_B275.csv")
data5 <- read_csv("01_Expe_EGFP_PAT_BPK/03_FCS_Analysis/2024_BPK_05_PAT_D6-D10_B275.csv")

L1<-list(data1=data1, data2=data2, data3=data3, data4=data4, data5=data5)
L1b<-list()
for (data_name in names(L1)) {
  LX<-L1[[data_name]][-((nrow(L1[[data_name]])-1):nrow(L1[[data_name]])),]
  LX<-LX[,-5]
  colnames(LX)<- c("name","per_TrueAlive","median_EGFP","count")
  LX[[3]]<-as.numeric(LX[[3]])
  LX<-LX[!grepl("_004|_008|_012|_016|_020|_024|_028|_032|_036|_040|_044|_048|_052", LX$name),]
  LX<-LX[!grepl("_056|_060|_064|_068|_072|_076", LX$name),]
  if (grepl(c("data1|data2"), data_name)){
    LX$BPK<-rep(c("BPK026","BPK031","BPK156","BPK080","BPK190","BPK087","BPK085",
                  "BPK282","BPK294","LdBOB"),each=6,times=3)
  } else if (grepl("data3", data_name)){
    LX$BPK<-rep(c("BPK026","BPK031","BPK156","BPK080","BPK190","BPK087","BPK085",
                  "BPK282","BPK294","LdBOB"),each=6,times=1)
  } else if (grepl("data4", data_name)){
    LX$BPK<-rep("BPK275", each=6, times=4)
  } else if (grepl("data5", data_name)){
    LX$BPK<-rep("BPK275", each=6, times=3)
  } 
  if (grepl("data1", data_name)){
    LX$day<-rep(c(1,2,3),each=60)
  } else if (grepl("data2", data_name)){
    LX$day<-rep(c(4,6,8), each=60)
  } else if (grepl("data3", data_name)){
    LX$day<-rep(c(10), each=60)
  } else if (grepl("data4", data_name)) {
    LX$day<-rep(c(1,2,3,4), each=6)
  } else if (grepl("data5", data_name)) {
    LX$day<-rep(c(6,8,10), each=6)
  }
  L1b[[data_name]]<-LX
  rm(LX)
}

data6<-rbind(L1b[["data1"]],L1b[["data2"]],L1b[["data3"]],L1b[["data4"]],L1b[["data5"]])

L2<-split(data6, data6$BPK)
L2b<-list()
for (BPK_name in names(L2)) {
  LX<-L2[[BPK_name]]
  LX$drug<-rep(c("D","ND"),each=3, times=7)
  LX$rep<-rep(c("R1","R2","R3"), times=14)
  if (grepl("BPK026|BPK031|BPK156", BPK_name)){
    LX$ISC<-"ISC1"
  } else if (grepl("BPK080|BPK190|BPK087", BPK_name)){
    LX$ISC<-"ISC4"
  } else if (grepl("BPK085|BPK282", BPK_name)){
    LX$ISC<-"ISC6"
  } else if (grepl("BPK294", BPK_name)){
    LX$ISC<-"ISC9"
  } else if (grepl("BPK275", BPK_name)){
    LX$ISC<-"ISC5"
  } else if (grepl("LdBOB", BPK_name)){
    LX$ISC<-"Other"
  }
  L2b[[BPK_name]]<-LX
  rm(LX)
}

data7<-rbind(L2b[["BPK026"]],L2b[["BPK031"]],L2b[["BPK156"]],L2b[["BPK080"]],
             L2b[["BPK190"]],L2b[["BPK087"]],L2b[["BPK085"]],L2b[["BPK282"]],
             L2b[["BPK294"]],L2b[["BPK275"]],L2b[["LdBOB"]])

rm(data1,data2,data3,data4,data5,data6,L1,L1b,L2,L2b)

write_csv(data7,"01_Expe_EGFP_PAT_BPK/02_R_files/data_f_1.csv")


####Supp_Table2_viability Analysis under PAT pressure####
#Loading the data#
data8<-read.csv("01_Expe_EGFP_PAT_BPK/02_R_files/data_f_1.csv")
head(data8)

#adding viability for day 0 at 100% and filtrate the data under PAT pressure only# 
L1<-split(data8, data8$day)
L1$"0"<-L1$"1"
L1$"0"$day<-0
L1$"0"$per_TrueAlive<-100

data8b<-do.call(rbind,L1)

L2<-split(data8b,data8b$drug)
data9<-L2$D
head(data9)

write_csv(data9,"01_Expe_EGFP_PAT_BPK/02_R_files/data_f_2.csv")

#Descriptive Statistics#
data10<-data9%>%
  dplyr::group_by(BPK,day,ISC)%>%
  dplyr::summarise_at(vars(c(per_TrueAlive)),list(mean=mean,sd=sd,se=std.error))
data10

write_xlsx(data10,"01_Expe_EGFP_PAT_BPK/02_R_files/viability_results.xlsx")
rm(data8,data8b,data9,data10,data7,L1,L2)


####Supp_Fig2_Number of Parasites####
##preparing the data##
data1 <- read_csv("01_Expe_EGFP_PAT_BPK/04_NumberParasite/2024_BPK_01_PAT_D1-D3_cells.csv")
data2 <- read_csv("01_Expe_EGFP_PAT_BPK/04_NumberParasite/2024_BPK_02_PAT_D4-D8_cells.csv")
data3 <- read_csv("01_Expe_EGFP_PAT_BPK/04_NumberParasite/2024_BPK_03_PAT_D10_cells.csv")
data4 <- read_csv("01_Expe_EGFP_PAT_BPK/04_NumberParasite/2024_BPK_04_PAT_D1-D4_B275_cells.csv")
data5 <- read_csv("01_Expe_EGFP_PAT_BPK/04_NumberParasite/2024_BPK_05_PAT_D6-D10_B275_cells.csv")

L1<-list(data1=data1, data2=data2, data3=data3, data4=data4, data5=data5)
L1b<-list()
for (data_name in names(L1)) {
  LX<-L1[[data_name]]
  colnames(LX)<- c("name","events","para")
  LX[[3]]<-as.numeric(LX[[3]])
  LX<-LX[!grepl("_004|_008|_012|_016|_020|_024|_028|_032|_036|_040|_044|_048|_052", LX$name),]
  LX<-LX[!grepl("_056|_060|_064|_068|_072|_076", LX$name),]
  if (grepl(c("data1|data2"), data_name)){
    LX$BPK<-rep(c("BPK026","BPK031","BPK156","BPK080","BPK190","BPK087","BPK085",
                  "BPK282","BPK294","LdBOB"),each=6,times=3)
  } else if (grepl("data3", data_name)){
    LX$BPK<-rep(c("BPK026","BPK031","BPK156","BPK080","BPK190","BPK087","BPK085",
                  "BPK282","BPK294","LdBOB"),each=6,times=1)
  } else if (grepl("data4", data_name)){
    LX$BPK<-rep("BPK275", each=6, times=4)
  } else if (grepl("data5", data_name)){
    LX$BPK<-rep("BPK275", each=6, times=3)
  } 
  if (grepl("data1", data_name)){
    LX$day<-rep(c(1,2,3),each=60)
  } else if (grepl("data2", data_name)){
    LX$day<-rep(c(4,6,8), each=60)
  } else if (grepl("data3", data_name)){
    LX$day<-rep(c(10), each=60)
  } else if (grepl("data4", data_name)) {
    LX$day<-rep(c(1,2,3,4), each=6)
  } else if (grepl("data5", data_name)) {
    LX$day<-rep(c(6,8,10), each=6)
  }
  L1b[[data_name]]<-LX
  rm(LX)
}

data6<-rbind(L1b[["data1"]],L1b[["data2"]],L1b[["data3"]],L1b[["data4"]],L1b[["data5"]])

L2<-split(data6, data6$BPK)
L2b<-list()
for (BPK_name in names(L2)) {
  LX<-L2[[BPK_name]]
  LX$drug<-rep(c("D","ND"),each=3, times=7)
  LX$rep<-rep(c("R1","R2","R3"), times=14)
  if (grepl("BPK026|BPK031|BPK156", BPK_name)){
    LX$ISC<-"ISC1"
  } else if (grepl("BPK080|BPK190|BPK087", BPK_name)){
    LX$ISC<-"ISC4"
  } else if (grepl("BPK085|BPK282", BPK_name)){
    LX$ISC<-"ISC6"
  } else if (grepl("BPK294", BPK_name)){
    LX$ISC<-"ISC9"
  } else if (grepl("BPK275", BPK_name)){
    LX$ISC<-"ISC5"
  } else if (grepl("LdBOB", BPK_name)){
    LX$ISC<-"Other"
  }
  L2b[[BPK_name]]<-LX
  rm(LX)
}

data7<-rbind(L2b[["BPK026"]],L2b[["BPK031"]],L2b[["BPK156"]],L2b[["BPK080"]],
             L2b[["BPK190"]],L2b[["BPK087"]],L2b[["BPK085"]],L2b[["BPK282"]],
             L2b[["BPK294"]],L2b[["BPK275"]],L2b[["LdBOB"]])

rm(data1,data2,data3,data4,data5,data6,L1,L1b,L2,L2b)

write_csv(data7,"01_Expe_EGFP_PAT_BPK/04_NumberParasite/02_R_files/data_f_4_cells.csv")

##preparing plot for number of parasites##
#Loading the data#
data8<-read.csv("01_Expe_EGFP_PAT_BPK/04_NumberParasite/02_R_files/data_f_4_cells.csv")

#adding viability for day 0 at 100% and filtrate the data under PAT pressure only# 
L1<-split(data8, data8$day)
L1$"0"<-L1$"1"
L1$"0"$day<-0
L1$"0"$para<-500000

data8b<-do.call(rbind,L1)
tail(data8b)  

L2<-split(data8b,data8b$drug)
data9<-L2$D

data9b<-L2$ND[L2$ND$day %in% c(0,1,2,3,4),]

data9<-rbind(data9,data9b)

#adding grouping strategy
data9[data9$BPK %in% c("BPK026","BPK031","BPK156"),"profile"]<-"P1"
data9[data9$BPK %in% c("BPK087","BPK190","LdBOB"),"profile"]<-"P2"
data9[data9$BPK %in% c("BPK085","BPK282","BPK294"),"profile"]<-"P3"
data9[data9$BPK %in% c("BPK275","BPK080"),"profile"]<-"P4"

write_csv(data9,"01_Expe_EGFP_PAT_BPK/04_NumberParasite/02_R_files/data_f_5_cells_filter.csv")

#Descriptive statistics and plot for number of cells #
data9[data9$events <10,"para"]<-NA
data10<-data9%>%
  dplyr::group_by(BPK,day,ISC,profile,drug)%>%
  dplyr::summarise_at(vars(c(para)),list(mean=~mean(.x,na.rm=TRUE),
                                         sd=~sd(.x,na.rm=TRUE),
                                         se=~std.error(.x, na.rm=TRUE)))

#factorizing data for plots#
data10$ISC<-factor(data10$ISC, levels = c("ISC1","ISC4","Other","ISC6","ISC9","ISC5"), ordered = TRUE)
data10$BPK<-factor(data10$BPK, levels = c("BPK026","BPK031","BPK156","BPK190","BPK087","LdBOB","BPK080","BPK085",
                                          "BPK282","BPK294","BPK275"))

#plot for number of parasite#BLACK AND WHITE WITH NO PAT
p10<-ggplot(data10,aes(day,mean))+
  geom_point(aes(shape = drug),size=2,colour="black")+
  geom_smooth(aes(linetype = drug), color="black", size=0.4, se = FALSE) +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se),width=0.2)+
  scale_shape_manual(name="",values = c("D"=16,"ND"=0), label = c("D"="PAT", "ND"="No PAT"))+
  scale_x_continuous(breaks=seq(0,10,by=1))+
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x)),
                limits = c(1e0,1e8))+
  facet_wrap(~BPK, scales = "free", nrow = 4)+
  theme(axis.line=element_line(), strip.text = element_text(size=12),
        strip.background = element_blank(),
        axis.text.x = element_text(face="bold", size=10),
        axis.text.y = element_text(face="bold", size=10),
        legend.title = element_text(size=9, face="bold"),
        legend.position = "right",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill="white"),
        axis.title = element_text(size = 11, face = "bold"),
        plot.title = element_text(size = 16, face = "bold"))+
  guides(x=guide_axis(cap = "both",minor.ticks = TRUE), 
         y=guide_axis(cap="both",minor.ticks = TRUE),
         fill=guide_legend(nrow = 11))+
  labs(title= " ",y="Number of viable parasites / mL", x="Time (day)")

p10

#saving the plot#
ggsave(file="01_Expe_EGFP_PAT_BPK/04_NumberParasite/02_R_files/Rplot20_NumberPara_D10.pdf", width = 6.5, height = 7)
ggsave(file="01_Expe_EGFP_PAT_BPK/04_NumberParasite/02_R_files/Rplot20_NumberPara_D10.tif", width = 4, height = 4)

rm(data9,data9b,data8,data8b,data7,data10,L1,L2,p10)


####Fig1_Delta log10 values####
data9<-read.csv("01_Expe_EGFP_PAT_BPK/04_NumberParasite/02_R_files/data_f_5_cells_filter.csv")
data11<-data9
data11<-data9[data9$drug=="D",]

L1<-split(data11,data11$day)
L1b<-L1[['0']]$para

L12<- lapply(L1, function(x) {
  x$delta <-log10(x$para/L1b)
  return(x)
})

data12<-do.call(rbind, L12)

#to decrease events lo less to 10 
data12[data12$events <10,"delta"]<-NA

#to save the files
write_xlsx(data12,"01_Expe_EGFP_PAT_BPK/04_NumberParasite/02_R_files/Delta_NumberParasites.xlsx")

#Descriptive statistics# 
data13<-data12%>%
  dplyr::group_by(BPK,day,ISC,profile)%>%
  dplyr::summarise_at(vars(c(para,delta)),list(mean=~mean(.x,na.rm=TRUE),
                                               sd=~sd(.x,na.rm=TRUE),
                                               se=~std.error(.x, na.rm=TRUE)))
data13


#factorizing data for plots#
data13$ISC<-factor(data13$ISC, levels = c("ISC1","ISC4","Other","ISC6","ISC9","ISC5"), ordered = TRUE)
data13$BPK<-factor(data13$BPK, levels = c("BPK026","BPK031","BPK156","BPK190","BPK087","LdBOB","BPK080","BPK085",
                                          "BPK282","BPK294","BPK275"))

#plot based on profile#
p11<-ggplot(data13,aes(day,delta_mean,fill=BPK))+
  geom_point(aes(colour=BPK,fill=BPK),shape=22,size=3,colour="black")+
  geom_errorbar(aes(ymin=delta_mean-delta_se, ymax=delta_mean+delta_se),width=0.2)+
  scale_x_continuous(breaks=seq(0,10,by=1))+
  scale_y_continuous(breaks =seq(-4,5, by=1),limits = c(-4.2,1))+  
  geom_smooth(method = "loess", se = FALSE, color="#636363", size=0.25)+
  scale_fill_manual(name="Ldo strain ",values = c("#a50026","#d73027","#f46d43","#053061","#2166ac",
                                                  "#f7f7f7","#4393c3","#8073ac", "#542788","#000000","#feb24c"))+
  facet_wrap(~profile, scales = "free")+
  theme(axis.line=element_line(), strip.text = element_text(size=12, face ="bold"),
        strip.background = element_blank(),
        axis.text.x = element_text(face="bold", size=6),
        axis.text.y = element_text(face="bold", size=10),
        legend.title = element_text(size=9, face="bold"),
        legend.position = "right",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill="white"),
        axis.title = element_text(size = 11, face = "bold"),
        plot.title = element_text(size = 16, face = "bold"))+
  guides(x=guide_axis(cap = "both",minor.ticks = TRUE), 
         y=guide_axis(cap="both",minor.ticks = TRUE),
         fill=guide_legend(nrow = 11))+
  labs(title= " ",y="Δ log 10 from initial concentration", x="Time (day)")

p11

#saving the plot#
ggsave(file="01_Expe_EGFP_PAT_BPK/04_NumberParasite/02_R_files/Rplot18_Delta10.pdf", width = 6.5, height = 5)
ggsave(file="01_Expe_EGFP_PAT_BPK/04_NumberParasite/02_R_files/Rplot18_Delta10.tif", width = 6.5, height = 5)


#plot based on BPK line#
p12<-ggplot(data13,aes(day,delta_mean,fill=BPK))+
  geom_point(aes(colour=BPK,fill=BPK),shape=22,size=3,colour="black")+
  geom_errorbar(aes(ymin=delta_mean-delta_se, ymax=delta_mean+delta_se),width=0.2)+
  scale_x_continuous(breaks=seq(0,10,by=1))+
  scale_y_continuous(breaks =seq(-4,5, by=1),limits = c(-4.2,1))+  
  geom_smooth(method = "loess", se = FALSE, color="#636363", size=0.25)+
  scale_fill_manual(name="Ldo strain ",values = c("#a50026","#d73027","#f46d43","#053061","#2166ac",
                                                  "#f7f7f7","#4393c3","#8073ac", "#542788","#000000","#feb24c"))+
  facet_wrap(~BPK, scales = "free")+
  theme(axis.line=element_line(), strip.text = element_text(size=12, face ="bold"),
        strip.background = element_blank(),
        axis.text.x = element_text(face="bold", size=6),
        axis.text.y = element_text(face="bold", size=10),
        legend.title = element_text(size=9, face="bold"),
        legend.position = "right",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill="white"),
        axis.title = element_text(size = 11, face = "bold"),
        plot.title = element_text(size = 16, face = "bold"))+
  guides(x=guide_axis(cap = "both",minor.ticks = TRUE), 
         y=guide_axis(cap="both",minor.ticks = TRUE),
         fill=guide_legend(nrow = 11))+
  labs(title= " ",y="Δ log 10 from initial concentration", x="Time (day)")

p12

#saving the plot#
ggsave(file="01_Expe_EGFP_PAT_BPK/04_NumberParasite/02_R_files/Rplot19_Delta10_BPK.pdf", width = 7.5, height = 5)
ggsave(file="01_Expe_EGFP_PAT_BPK/04_NumberParasite/02_R_files/Rplot19_Delta10_BPK.tif", width = 7.5, height = 5)




####Fig2A_Percentage of reduction of EGFP####
#load the data#
data8<-read.csv("01_Expe_EGFP_PAT_BPK/02_R_files/data_f_1.csv")
head(data8)

#prepare the data#
L2<-split(data8,data8$drug)
data14<-L2$ND[L2$ND$day %in% c(2,4),]
data14b<-L2$D[L2$D$day %in% c(2,4,6,8,10),]
head(data14)

data14$stage<-ifelse(data14$day==2, "log",
                     ifelse(data14$day==4, "sta", NA))

data14b$stage<-ifelse(data14b$day==2, "PAT_2",
                      ifelse(data14b$day==4, "PAT_4",
                             ifelse(data14b$day==6, "PAT_6",
                                    ifelse(data14b$day==8, "PAT_8",
                                           ifelse(data14b$day==10,"PAT_10", NA)))))


data15<-rbind(data14,data14b)

head(data15)
write_csv(data15,"01_Expe_EGFP_PAT_BPK/02_R_files/data_f_3_EGFP.csv")
rm(data14, data14b,data8,L2)

#calculate the ratio#
data8<-data15
rm(data15)
L2<-split(data8,data8$drug)
data14<-L2$ND[L2$ND$day==2,]
data14b<-L2$D 
data14b<-rbind(L2$D,L2$ND[L2$ND$day==4,])
data14c<-split(data14b,data14b$day)

L1b<-list()

for (day_name in names(data14c)) {
  LX<-data14c[[day_name]]
  LX$ratio<-LX$median_EGFP/data14$median_EGFP
  LX$red<-100-(LX$ratio * 100)
  L1b[[day_name]]<-LX
}

data15<-do.call(rbind,L1b)

#calculate the Percentage of reduction, lower than 40 events##
data15b<-data15[data15$count>=40,]
data15b<-data15b[data15b$day %in% c(2,4),]

write_csv(data15b,"01_Expe_EGFP_PAT_BPK/02_R_files/data_f_4_PerEGFP.csv")

#Descriptive statistics#
data16<-data15b%>%
  dplyr::group_by(BPK,day,ISC,stage,drug)%>%
  dplyr::summarise_at(vars(c(per_TrueAlive,ratio,red)), 
                      list(mean=mean, sd=sd, se=std.error))
data16
write_csv(data16,"01_Expe_EGFP_PAT_BPK/02_R_files/data_f_4_Statistics.csv")

#adding the NA data to the plot#
b<- data.frame(
  BPK = c("BPK190","BPK190","BPK087","BPK085","BPK282","BPK294","LdBOB","LdBOB"),
  day = c(2,4,4,4,4,4,2,4),
  ISC = c("ISC4","ISC4","ISC4","ISC6","ISC6","ISC9", "Other","Other"),
  stage = c("PAT_2", "PAT_4","PAT_4", "PAT_4", "PAT_4","PAT_4","PAT_2","PAT_4"),
  drug = c("D", "D","D", "D", "D","D","D","D"),
  per_TrueAlive_mean = c(0,0,0,0,0,0,0,0),
  ratio_mean = c(0,0,0,0,0,0,0,0),
  red_mean = c(0,0,0,0,0,0,0,0),
  per_TrueAlive_sd = c(0,0,0,0,0,0,0,0),
  ratio_sd = c(0,0,0,0,0,0,0,0),
  red_sd = c(0,0,0,0,0,0,0,0),
  per_TrueAlive_se = c(0,0,0,0,0,0,0,0),
  ratio_se = c(0,0,0,0,0,0,0,0),
  red_se = c(0,0,0,0,0,0,0,0)
)

b$day<-as.integer(b$day)
data16<-rbind(data16,b)

#plot the ratio of EGFP_median#
data16[data16$BPK %in% c("BPK026","BPK031","BPK156"),"profile"]<-"P1"
data16[data16$BPK %in% c("BPK190","BPK087","LdBOB"),"profile"]<-"P2"
data16[data16$BPK %in% c("BPK085","BPK282","BPK294"),"profile"]<-"P3"
data16[data16$BPK %in% c("BPK275","BPK080"),"profile"]<-"P4"

data16$BPK<-factor(data16$BPK, levels = c("BPK026","BPK031","BPK156","BPK190","BPK087","BPK080","LdBOB","BPK085",
                                          "BPK282","BPK294","BPK275"))
data16$stage<-factor(data16$stage, levels = c("log","sta","PAT_2","PAT_4"))

data16a<-data16[data16$profile == "P1",]
data16b<-data16[data16$profile == "P2",]
data16c<-data16[data16$profile == "P3",]
data16d<-data16[data16$profile == "P4",]


p9a<-ggplot(data16a,aes(BPK,red_mean,fill=stage))+
  geom_col(position = 'dodge', 
           color='black')+
  geom_errorbar(aes(ymin=red_mean-red_se, ymax=red_mean+red_se),width=0.3,
                position=position_dodge(0.9))+
  scale_fill_manual(name=" ",values = c("#fc8d59","#91bfdb","#f7f7f7"),labels=c("No PAT D4","PAT D2","PAT D4"))+
  scale_x_discrete(guide = guide_axis(angle=90))+
  scale_y_continuous(breaks = seq(-100,100,by=20), limits = c(-100,100))+
  facet_wrap(~profile)+
  theme_classic()+
  theme(axis.line = element_line(), panel.background = element_rect(fill = "white"), 
        axis.text.x = element_text(face = "bold", size = 10), 
        axis.text.y = element_text(face = "bold", size = 10),
        axis.title.y = element_text(face = "bold", size = 10),
        legend.text = element_text(face = "bold", size = 9),
        legend.position = "none",
        strip.background = element_blank(),  
        strip.text = element_text(face = "bold", size = 14),
        plot.title = element_text(size = 19, face = "bold", hjust = 0.5)) +
  guides(y = guide_axis(cap = "both", minor.ticks = TRUE)) + 
  labs(title = " ", x = " ", y = "rEGFP reduction (%)")

p9a

p9b<-ggplot(data16b,aes(BPK,red_mean,fill=stage))+
  geom_col(position = 'dodge', 
           color='black')+
  geom_errorbar(aes(ymin=red_mean-red_se, ymax=red_mean+red_se),width=0.3,
                position=position_dodge(0.9))+
  scale_fill_manual(name=" ",values = c("#fc8d59","#91bfdb","#f7f7f7"),labels=c("No PAT D4","PAT D2","PAT D4"))+
  scale_x_discrete(guide = guide_axis(angle=90))+
  scale_y_continuous(breaks = seq(-100,100,by=20), limits = c(-100,100))+
  facet_wrap(~profile)+
  theme_classic()+
  theme(axis.line = element_line(), panel.background = element_rect(fill = "white"), 
        axis.text.x = element_text(face = "bold", size = 10), 
        axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        legend.text = element_text(face = "bold", size = 9),
        legend.position = "none",
        strip.background = element_blank(),  
        strip.text = element_text(face = "bold", size = 14),
        plot.title = element_text(size = 19, face = "bold", hjust = 0.5)) +
  guides(y = guide_axis(cap = "both", minor.ticks = TRUE)) + 
  labs(title = " ", x = " ", y = "rEGFP reduction (%)")

p9b
p9b <- p9b +
  annotate("text", x=1.05, y = 10.5, label="n.d.", size=2)+
  annotate("text", x=1.32, y = 10.5, label="n.d.", size=2)+
  annotate("text", x=2.30, y = 10.5, label="n.d.", size=2)+
  annotate("text", x=3.05, y = 10.5, label="n.d.", size=2)+
  annotate("text", x=3.30, y = 10.5, label="n.d.", size=2)
p9b

p9c<-ggplot(data16c,aes(BPK,red_mean,fill=stage))+
  geom_col(position = 'dodge', 
           color='black')+
  geom_errorbar(aes(ymin=red_mean-red_se, ymax=red_mean+red_se),width=0.3,
                position=position_dodge(0.9))+
  scale_fill_manual(name=" ",values = c("#fc8d59","#91bfdb","#f7f7f7"),labels=c("No PAT D4","PAT D2","PAT D4"))+
  scale_x_discrete(guide = guide_axis(angle=90))+
  scale_y_continuous(breaks = seq(-100,100,by=20), limits = c(-100,100))+
  facet_wrap(~profile)+
  theme_classic()+
  theme(axis.line = element_line(), panel.background = element_rect(fill = "white"), 
        axis.text.x = element_text(face = "bold", size = 10), 
        axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        legend.text = element_text(face = "bold", size = 9),
        legend.position = "none",
        strip.background = element_blank(),  
        strip.text = element_text(face = "bold", size = 14),
        plot.title = element_text(size = 19, face = "bold", hjust = 0.5)) +
  guides(y = guide_axis(cap = "both", minor.ticks = TRUE)) + 
  labs(title = " ", x = " ", y = "rEGFP reduction (%)")
p9c

p9c <- p9c +
  annotate("text", x=1.30, y = 10.5, label="n.d.", size=2)+
  annotate("text", x=2.31, y = 10.5, label="n.d.", size=2)+
  annotate("text", x=3.31, y = 10.5, label="n.d.", size=2)
p9c


p9d<-ggplot(data16d,aes(BPK,red_mean,fill=stage))+
  geom_col(position = 'dodge', 
           color='black')+
  geom_errorbar(aes(ymin=red_mean-red_se, ymax=red_mean+red_se),width=0.3,
                position=position_dodge(0.9))+
  scale_fill_manual(name=" ",values = c("#fc8d59","#91bfdb","#f7f7f7"),labels=c("No PAT D4","PAT D2","PAT D4"))+
  scale_x_discrete(guide = guide_axis(angle=90))+
  scale_y_continuous(breaks = seq(-100,100,by=20), limits = c(-100,100))+
  facet_wrap(~profile)+
  theme_classic()+
  theme(axis.line = element_line(), panel.background = element_rect(fill = "white"), 
        axis.text.x = element_text(face = "bold", size = 10), 
        axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        legend.text = element_text(face = "bold", size = 9),
        legend.position = "right",
        strip.background = element_blank(),  
        strip.text = element_text(face = "bold", size = 14),
        plot.title = element_text(size = 19, face = "bold", hjust = 0.5)) +
  guides(y = guide_axis(cap = "both", minor.ticks = TRUE)) + 
  labs(title = " ", x = " ", y = "rEGFP reduction (%)")

p9d

#combine the plots
combi_plot<-ggdraw()+
  draw_plot(p9a, x = 0.0, y = 0.2, width = 0.28, height = 0.75)+
  draw_plot(p9b, x = 0.28, y = 0.2, width = 0.22, height = 0.75)+
  draw_plot(p9c, x = 0.50, y = 0.2, width = 0.22, height = 0.75)+
  draw_plot(p9d, x = 0.71, y = 0.2, width = 0.30, height = 0.75)

combi_plot

#save the pot#
ggsave(file="01_Expe_EGFP_PAT_BPK/02_R_files/Rplot04_EGFPreduction_figure02a.pdf", width = 11, height = 5.5)
ggsave(file="01_Expe_EGFP_PAT_BPK/02_R_files/Rplot04_EGFPreduction_figure02a.tif", width = 11, height = 5.5)

rm(b,combi_plot,data14,data14b,data14c,data15,data15b,data16,data16a,data16c,data16d,data8)
rm(data16b,L1b,L2,LX,p9a,p9b,p9c,p9d)



####Figu2a
####Prepare the data for Figure 2B and Supp_Fig4b####
data1 <- read_csv("02_Expe_GC_PAT_B026_B275/03_FCS_Analysis/2024_BPK_GC.csv")
head(data1)

L1<-list(data1=data1)
L1b<-list()
for (data_name in names(L1)) {
  LX<-L1[[data_name]][-((nrow(L1[[data_name]])-1):nrow(L1[[data_name]])),]
  LX<-LX[,-5]
  colnames(LX)<- c("name","per_TrueAlive","median_EGFP","count")
  LX[[3]]<-as.numeric(LX[[3]])
  LX$BPK <- ifelse(grepl("B026", LX$name), "BPK026",
                   ifelse(grepl("B275", LX$name), "BPK275", NA))
  LX$day <- ifelse(grepl("-D0", LX$name), "0",
                   ifelse(grepl("D1.fcs", LX$name), "1",
                          ifelse(grepl("D2.fcs", LX$name), "2",
                                 ifelse(grepl("D3.fcs", LX$name), "3",
                                        ifelse(grepl("D4.fcs", LX$name), "4",
                                               ifelse(grepl("D5.fcs", LX$name), "5",
                                                      ifelse(grepl("D6.fcs", LX$name), "6",
                                                             ifelse(grepl("D7.fcs", LX$name),"7",
                                                                    ifelse(grepl("D10.fcs", LX$name),"10",NA)))))))))
  LX$drug<-ifelse(grepl("-ND", LX$name), "No_PAT",
                  ifelse(grepl("-D", LX$name),"PAT", NA))
  L1b[[data_name]]<-LX
  rm(LX)
}
L1b

data2<-do.call(rbind, L1b)

grep("D7_2", data2$name)
data2[c(79,80,81),"day"] <- "7"

write_csv(data2,"02_R_files/data_GC_1.csv")


####Fig2B_EGFPreduction####
#prepare the data#
data3<-read.csv("02_R_files/data_GC_1.csv")

data3[is.na(data3)]<-7

data3[data3$drug == "No_PAT","drug"]<-"ND"
data3[data3$drug == "PAT","drug"]<-"D"

L2<-split(data3,data3$drug)
data14<-L2$ND[L2$ND$day==2,]
data14b<-L2$ND[L2$ND$day==0,]
data14c<-rbind(L2$D,data14b) 

data14d<-split(data14c,data14c$day)

L1b<-list()

for (day_name in names(data14d)) {
  LX<-data14d[[day_name]]
  LX$ratio<-LX$median_EGFP/data14$median_EGFP
  LX$red<-100-(LX$ratio * 100)
  L1b[[day_name]]<-LX
}

data15<-do.call(rbind,L1b)
data4<-data15
head(data15)

write_csv(data15,"02_R_files/data_GC_2.csv")


#some statistics
data4a<-data4%>%
  dplyr::group_by(BPK,day,drug)%>%
  dplyr::summarise_at(vars(c(per_TrueAlive,median_EGFP,red)),
                      list(mean=mean,sd=sd,se=std.error))

data4a

#generating the plot
p1<-ggplot(data4a,aes(day,red_mean))+
  geom_point(aes(shape = BPK),size=3, color="black")+
  geom_line(aes(linetype=BPK), color="black", size=0.4) +
  geom_errorbar(aes(ymin = red_mean-red_se, 
                    ymax=red_mean+red_se), width=0.3)+
  scale_x_continuous(limits = c(0,10),breaks = seq(0,10, by=1))+
  scale_y_continuous(limits = c(-2,105), breaks = seq(0,105,by=10))+
  scale_shape_manual(name="Ldo strain",values = c("BPK026"=16,"BPK275"=0))+
  scale_linetype_manual(values = c("BPK026"="solid","BPK275"="dashed"))+
  theme(axis.line = element_line(),
        panel.background = element_rect(fill = "white"),
        axis.text.x = element_text(face = "bold", size=10),
        axis.text.y = element_text(face = "bold", size=10),
        axis.title.x = element_text(face="bold",size=12),
        axis.title.y = element_text(face="bold",size=12),
        strip.background = element_rect(fill = "white",color = "black"),  
        legend.title = element_text(face = "bold", size=10),
        legend.text = element_text(size=9),
        legend.position = "right")+
  guides(x=guide_axis(cap="both",minor.ticks = TRUE),
         y=guide_axis(cap="both", minor.ticks = TRUE))+
  labs(title = " ",x="Time (day)", 
       y="rEGFP reduction (%)")

p1


#get the legend#
legendC<-get_legend(p1 + theme(legend.position = "right",
                               legend.key.size = unit(0.5,"cm"),
                               legend.text = element_text(size=8),
                               legend.title = element_text(size=10)))


#saving the plot
ggsave(file="02_R_files/Rplot02_EGFP_red_fig2b.pdf", width = 5.5, height = 3.7)
ggsave(file="02_R_files/Rplot02_EGFP_red_fig2b.tif", width = 5.5, height = 3.7)

####Supp_Fig4b_plot Viability####
data3<-read.csv("02_R_files/data_GC_1.csv")
head(data3)

data3<-data3[data3$count>40,]
data3b<-data3[data3$day==0 & data3$drug=="No_PAT",]
data3b$drug<-"PAT"

data4<-rbind(data3,data3b)
head(data4)
data4[is.na(data4)]<-7

write_csv(data4,"02_R_files/data_GC_via_3.csv")

#some statistics
data4<-data4%>%
  dplyr::group_by(BPK,day,drug)%>%
  dplyr::summarise_at(vars(c(per_TrueAlive,median_EGFP)),
                      list(mean=mean,sd=sd,se=std.error))

data4

#generating the plot
plot1<-ggplot(data4,aes(day,per_TrueAlive_mean))+
  geom_point(aes(shape = BPK),size=2, color="black")+
  geom_line(aes(linetype=drug,group=interaction(drug,BPK)), color="black", size=0.4) +
  geom_errorbar(aes(ymin = per_TrueAlive_mean-per_TrueAlive_se, 
                    ymax=per_TrueAlive_mean+per_TrueAlive_se), width=0.3)+
  scale_x_continuous(limits = c(0,10),breaks = seq(0,10, by=1))+
  scale_y_continuous(limits = c(-2,110), breaks = seq(0,110,by=20))+
  scale_shape_manual(name="Ldo strain",values = c("BPK026"=16,"BPK275"=0))+
  scale_linetype_manual(name="Drug",values = c("No_PAT"="dashed","PAT"="solid"))+
  theme(axis.line = element_line(),
        panel.background = element_rect(fill = "white"),
        axis.text.x = element_text(face = "bold", size=9),
        axis.text.y = element_text(face = "bold", size=9),
        axis.title.x = element_text(face="bold",size=9),
        axis.title.y = element_text(face="bold",size=9),
        strip.background = element_rect(fill = "white",color = "black"),  
        legend.title = element_text(face = "bold", size=10) )+
  guides(x=guide_axis(cap="both",minor.ticks = TRUE),
         y=guide_axis(cap="both", minor.ticks = TRUE))+
  labs(title = " ",x="Time (day)", 
       y="Viable promastigotes (%)")

plot1

#saving the plot
ggsave(file="02_R_files/Rplot02_Via_B026_B275.pdf", width = 4, height = 2.5)
ggsave(file="02_R_files/Rplot02_Via_B026_B275.tif", width = 4, height = 2.5)

####Supp_Fig4a_Number of Parasites####
#loading the data#
data5<-read.csv("02_Expe_GC_PAT_B026_B275/02_R_files/2024_BPK_Nparasites.csv")
head(data5)

data5b<-data5[data5$day==0 & data5$drug=="No_PAT",]
data5b$drug<-"PAT"

data6<-rbind(data5,data5b)

#some statistics
data7<-data6%>%
  dplyr::group_by(BPK,day,drug)%>%
  dplyr::summarise_at(vars(c(para)),
                      list(mean=mean,sd=sd,se=std.error))

data7

#generating the plot#
plot1<-ggplot(data7, aes(day,mean)) +
  geom_point(aes(shape = BPK),size=2,color="black")+
  geom_line(aes(linetype=drug,group=interaction(drug,BPK)), color="black", size=0.4) +
  geom_errorbar(aes(ymin = mean-se, 
                    ymax = mean+se), width=0.3)+
  scale_shape_manual(name="Ldo strain",values = c("BPK026"=16,"BPK275"=0))+
  scale_x_continuous(limits = c(0,10),breaks = seq(0,10, by=1))+
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x)),
                limits = c(1e2,1e8))+
  scale_linetype_manual(name="Drug",values = c("No_PAT"="dashed","PAT"="solid"))+
  theme(axis.line = element_line(),
        panel.background = element_rect(fill = "white"),
        axis.text.x = element_text(face = "bold", size=9),
        axis.text.y = element_text(face = "bold", size=9),
        axis.title.x = element_text(face="bold",size=9),
        axis.title.y = element_text(face="bold",size=9),
        legend.title = element_text(face = "bold", size=10))+
  guides(x=guide_axis(cap="both",minor.ticks = TRUE),
         y=guide_axis(cap="both", minor.ticks = TRUE))+
  labs(title= " ",y="Number of viable parasites / mL", x="Time (day)")

plot1

#saving the plot#
ggsave(file="02_R_files/Rplot03_Para_B026_B275.pdf", width = 4, height = 2.5)
ggsave(file="02_R_files/Rplot03_Para_B026_B275.tif", width = 4, height = 2.5)



####Fig2c_Proliferation_CELL_TRACKER####
##Preparing the data##
#loading the data
data9<-read.csv("04_Expe_CellTrack_PAT_B026_B275/03_FCS_Files/2024_BPK_CellTra.csv")
head(data9)

#preparing the data
L9<-list(data1=data9)
L9b<-list()
for (data_name in names(L9)) {
  LX<-L9[[data_name]][-((nrow(L9[[data_name]])-1):nrow(L9[[data_name]])),]
  LX<-LX[,-6]
  colnames(LX)<-c("name","per_TrueAlive","median_EGFP","count","median_CellTra")
  LX$BPK<- ifelse(grepl("B026_|BPK026",LX$name),"BPK026",
                  ifelse(grepl("B275_|BPK275", LX$name),"BPK275",NA))
  LX$day<-ifelse(grepl("0119",LX$name),"0",
                 ifelse(grepl("0423",LX$name), "4",
                        ifelse(grepl("0429",LX$name), "10",NA)))
  LX$drug<-ifelse(grepl("_T1_|_T2_|_T3_|_T4_|_T5_|_T6_",LX$name),"NoPAT",
                  ifelse(grepl("_T7_|_T8_|_T9_|_T10_|_T11_|_T12_",LX$name), "PAT",
                         ifelse(grepl("_1_",LX$name),"NoPAT", NA)))
  L9b[[data_name]]<-LX
}

data9<-do.call(rbind, L9b)
head(data9)

#preparing the data adding the stage column
data9<-data9[data9$count>40,]
data10<-data9
data10[data10$day==0 & data10$drug=="NoPAT","stage"]<-"D0"
data10[data10$day==4 & data10$drug=="NoPAT","stage"]<-"D4"
data10[data10$day==10 & data10$drug=="NoPAT","stage"]<-"D10"
data10[data10$day==4 & data10$drug=="PAT","stage"]<-"PAT_4"
data10[data10$day==10 & data10$drug=="PAT","stage"]<-"PAT_10"

head(data10)

#saving the data
write.csv(data10,"02_R_Files/data_f_CT_1.csv")

##statistical analysis##
#01#analysis by BPK using anova (analysis of D0,D4,D10,PAT_4,PAT_10
#From previous analysis, the data was homogeneous but not normal
#A log transformation was performed
#some values are 0, I converted to 1 first before log transformation

#loading the data
data10<-read.csv("02_R_Files/data_f_CT_1.csv")

#preparing the data and log transformation
d2<-data10
d2[d2$median_CellTra==0,"median_CellTra"] <- 1
d2$log_CellTra<-log10(d2$median_CellTra)

#saving the log data
write.csv(d2,"02_R_Files/data_f_CT_1_log.csv")

#be sure that data is numeric and factoring data 
d2$log_CellTra<-as.numeric(d2$log_CellTra)
d2$BPK<-factor(d2$BPK, levels = c("BPK026","BPK275"))
d2$stage<-factor(d2$stage, levels = c("D0","D4","D10","PAT_4","PAT_10"))

#statistical analysis
L4<-split(d2,d2$BPK)
results <- list()
for (BPK_name in names(L4)){
  LX<-L4[[BPK_name]]
  fit2<-lm(log_CellTra~stage,LX)
  anova_results <- Anova(fit2)
  residuals <- augment(fit2)$.resid
  normality_test <- ols_test_normality(fit2)
  shapirox<-shapiro.test(fit2$residuals)
  plot_residx<- ols_plot_resid_fit(fit2)
  plot_qqx<-ols_plot_resid_qq(fit2)
  bartx<-bartlett.test(log_CellTra~stage,data=LX)
  levex<-leveneTest(log_CellTra~stage,data=LX)
  Tukeyx<- glht(fit2,linfct = mcp(stage="Tukey"))
  Tukeyy<-summary(Tukeyx)
  results[[BPK_name]] <- list(
    fit=fit2,
    anova= anova_results,
    residuals = residuals,
    normality=normality_test,
    shapiro=shapirox,
    plot_resid=plot_residx,
    plot_qq=plot_qqx,
    barttest=bartx,
    levetest=levex,
    Tukey=Tukeyy)
}
results$BPK026$Tukey

#code for siginificance 
significance_codes <- function(p_value){
  if (p_value < 0.001) {
    return("***")
  } else if (p_value < 0.01) {
    return("**")
  } else if (p_value < 0.05) {
    return("*")
  } else if (p_value < 0.1) {
    return(".")
  } else {
    return("")
  } 
}

#loop for statistical analysis for excell using anova
results_list <- list()
for (BPK_name in names(L4)){
  LX<-L4[[BPK_name]]
  fit2<-lm(log_CellTra~stage,LX)
  anova_results <- Anova(fit2)
  residuals <- augment(fit2)$.resid
  normality_test <- ols_test_normality(fit2)
  shapirox<-shapiro.test(fit2$residuals)
  plot_residx<- ols_plot_resid_fit(fit2)
  plot_qqx<-ols_plot_resid_qq(fit2)
  bartx<-bartlett.test(log_CellTra~stage,data=LX)
  levex<-leveneTest(log_CellTra~stage,data=LX)
  Tukeyx<- glht(fit2,linfct = mcp(stage="Tukey"))
  Tukeyy<-summary(Tukeyx)
  Tukey_p_values <-Tukeyy$test$pvalues
  Tukey_significance <-sapply(Tukey_p_values,significance_codes)
  Tukey_results <-paste0(round(Tukey_p_values, 5), " ",Tukey_significance)
  results_df <- data.frame(
    BPK = BPK_name,
    anovar_p_value= anova_results$`Pr(>F)`[1],
    residuals=paste(round(residuals, 3), collapse = ","),
    shapiro_W= shapirox$statistic,
    shapiro_p_value= shapirox$p.value,
    Bartlett_p_value= bartx$p.value, 
    levene_p_value= levex$`Pr(>F)`[1],
    D4_D0= Tukey_results[1],
    D10_D0= Tukey_results[2],
    PAT4_D0= Tukey_results[3],
    PAT10_D0= Tukey_results[4],
    D10_D4= Tukey_results[5],
    PAT4_D4= Tukey_results[6],
    PAT10_D4= Tukey_results[7],
    PAT4_D10= Tukey_results[8],
    PAT10_D10= Tukey_results[9],
    PAT10_PAT4= Tukey_results[10]
  )
  results_list[[BPK_name]]<-results_df
}

#exporting the data of the analysis
final_results<- do.call(rbind, results_list)
write_xlsx(final_results,"02_R_Files/statistical_results_Ctra_TransData_Tukey.xlsx")

##descriptive statistics on the log data##
##loading the data for mean and plots##
d2<-read.csv("02_R_Files/data_f_CT_1_log.csv")
#mean,sd,se
data11<-d2%>%
  dplyr::group_by(BPK,day,drug,stage)%>%
  dplyr::summarise_at(vars(c(per_TrueAlive,median_EGFP,median_CellTra,log_CellTra)),
                      list(mean=mean,sd=sd,se=std.error))

data11$stage<-factor(data11$stage, 
                     levels = c("D0","D4","D10","PAT_4","PAT_10"), 
                     ordered = TRUE)
data11a<-data11[data11$BPK == "BPK026",]
data11b<-data11[data11$BPK == "BPK275",]

##generating BAR plot for BPK026##
plot10a<-ggplot(data = data11a, aes(x = stage, y = log_CellTra_mean, fill = drug))+
  geom_col(position = "dodge", colour = "black",width = 0.8)+
  geom_errorbar(aes(ymin = log_CellTra_mean - log_CellTra_se, 
                    ymax = log_CellTra_mean + log_CellTra_se), 
                width = 0.3, position = position_dodge(0.9), color = "black") +
  scale_fill_manual(name = " ", values = c("#525252", "#d9d9d9"),
                    labels = c("No PAT","PAT")) +
  scale_x_discrete(labels = c("D0" = "No PAT\nD0", 
                              "D4" = "No PAT\nD4", 
                              "D10" = "No PAT\nD10",
                              "PAT_4" = "PAT\nD4",
                              "PAT_10" = "PAT\nD10")) +
  scale_y_continuous(limits=c(0,5), breaks=c(0,1,2,3,4,5),
                     labels = scales::label_math(expr="1"~x10^.x))+
  theme(axis.line = element_line(), panel.background = element_rect(fill = "white"), 
        axis.text.x = element_text(face = "bold", size = 8), 
        axis.text.y = element_text(face = "bold", size = 10),
        axis.title.y = element_text(face="bold",size=10),
        legend.text = element_text(face = "bold", size = 9),
        legend.position = "bottom",
        strip.background = element_rect(fill = "white"),  
        strip.text = element_text(face = "bold", size = 10),
        plot.title = element_text(size = 12, face = "bold", hjust = 0.5)) +
  guides(y = guide_axis(cap = "both", minor.ticks = TRUE)) + 
  labs(title = "BPK026", x = " ", y = "Cell Tracker (MFI)")

plot10a

plot11a <- plot10a +
  annotate("text", x=2.00, y = 2.3, label="***", size=4)+
  annotate("text", x=3.00, y = 2.0, label="***", size=4)+
  annotate("text", x=4.00, y = 4.2, label="***", size=4)+
  annotate("text", x=5.00, y = 4.0, label="***", size=4)

plot11a

##generating BAR plot for BPK275##
plot10b<-ggplot(data = data11b, aes(x = stage, y = log_CellTra_mean, fill = drug))+
  geom_col(position = "dodge", colour = "black",width = 0.8)+
  geom_errorbar(aes(ymin = log_CellTra_mean - log_CellTra_se, 
                    ymax = log_CellTra_mean + log_CellTra_se), 
                width = 0.3, position = position_dodge(0.9), color = "black") +
  scale_fill_manual(name = " ", values = c("#525252", "#d9d9d9"),
                    labels = c("No PAT", "PAT")) +
  scale_x_discrete(labels = c("D0" = "No PAT\nD0", 
                              "D4" = "No PAT\nD4", 
                              "D10" = "No PAT\nD10",
                              "PAT_4" = "PAT\nD4",
                              "PAT_10" = "PAT\nD10")) +
  scale_y_continuous(limits=c(0,5), breaks=c(0,1,2,3,4,5),
                     labels = scales::label_math(expr="1"~x10^.x))+
  theme(axis.line = element_line(), panel.background = element_rect(fill = "white"), 
        axis.text.x = element_text(face = "bold", size = 8), 
        axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        legend.text = element_text(face = "bold", size = 9),
        legend.position = "bottom",
        strip.background = element_rect(fill = "white"),  
        strip.text = element_text(face = "bold", size = 10),
        plot.title = element_text(size = 12, face = "bold", hjust = 0.5)) +
  guides(y = guide_axis(cap = "both", minor.ticks = TRUE)) + 
  labs(title = "BPK275", x = " ", y = "Cell Tracker (MFI)")

plot10b

plot11b <- plot10b +
  annotate("text", x=2.00, y = 1.9, label="***", size=4)+
  annotate("text", x=3.00, y = 1.4, label="***", size=4)+
  annotate("text", x=4.00, y = 2.8, label="***", size=4)+
  annotate("text", x=5.00, y = 0.5, label="***", size=4)
plot11b

#combine the plots
combin_plot1<-ggdraw() +
  draw_plot(plot11a, x = 0, y = 0.2, width = 0.52, height = 0.8) +
  draw_plot(plot11b, x = 0.52, y = 0.2, width = 0.45, height = 0.8)

combin_plot1


#get the legend#
legendB<-get_legend(plot10b + theme(legend.position = "right",
                                    legend.key.size = unit(0.5,"cm"),
                                    legend.text = element_text(size=8),
                                    legend.title = element_text(size=10)))

#export figures
ggsave(file="02_R_files/Rplot01_CellTra_fig2c.pdf", width = 6, height = 5)
ggsave(file="02_R_files/Rplot01_CellTra_fig2c.tif", width = 6, height = 5)




####Fig2D_MitochondrialMembranePotential####
data6<-read.csv("03_Expe_MiPo_PAT_B026_B275/03_FCS_files/2024_BPK_MiPo_A.csv")
data7<-read.csv("03_Expe_MiPo_PAT_B026_B275/03_FCS_files/2024_BPK_MiPo_B.csv")

head(data6)

L8<-list(data1=data6, data2=data7)
L8b<-list()
for (data_name in names(L8)) {
  LX<-L8[[data_name]][-((nrow(L8[[data_name]])-1):nrow(L8[[data_name]])),]
  LX<-LX[,-6]
  colnames(LX)<-c("name","per_TrueAlive","median_EGFP","count","median_MitoPo")
  LX$BPK<- ifelse(grepl("1109_",LX$name),"BPK026",
                  ifelse(grepl("1110_", LX$name),"BPK275",NA))
  LX$stage<-ifelse(grepl("T1_G|T2_G|T3_G",LX$name),"Pro_Log",
                   ifelse(grepl("T4_G|T5_G|T6_G",LX$name), "Pro_Sta_4",
                          ifelse(grepl("T7_G|T8_G|T9_G",LX$name), "Pro_Sta_10",
                                 ifelse(grepl("T10_G|T11_G|T12_G",LX$name),"Pro_PAT 2",
                                        ifelse(grepl("T13_G|T14_G|T15_G",LX$name), "Pro_PAT 4",
                                               ifelse(grepl("T16_G|T17_G|T18_G", LX$name), "Pro_PAT 10",NA))))))
  L8b[[data_name]]<-LX
}

data8<-do.call(rbind, L8b)
head(data8)

#save the data# 
write.csv(data8,"02_R_Files/data_f_Mi.csv")


#load the data#
data8<-read.csv("02_R_Files/data_f_Mi.csv")
d2<-data8
head(data8)
unique(d2$stage)

#be sure that data is numeric and factoring data 
d2$median_MitoPo<-as.numeric(d2$median_MitoPo)
d2$BPK<-factor(d2$BPK, levels = c("BPK026","BPK275"))
d2$stage<-factor(d2$stage, levels = c("Pro_Log","Pro_Sta_4","Pro_Sta_10","Pro_PAT 2","Pro_PAT 4","Pro_PAT 10"))

#statistical analysis
L4<-split(d2,d2$BPK)
results <- list()
for (BPK_name in names(L4)){
  LX<-L4[[BPK_name]]
  fit2<-lm(median_MitoPo~stage,LX)
  anova_results <- Anova(fit2)
  residuals <- augment(fit2)$.resid
  normality_test <- ols_test_normality(fit2)
  shapirox<-shapiro.test(fit2$residuals)
  plot_residx<- ols_plot_resid_fit(fit2)
  plot_qqx<-ols_plot_resid_qq(fit2)
  bartx<-bartlett.test(median_MitoPo~stage,data=LX)
  levex<-leveneTest(median_MitoPo~stage,data=LX)
  Tukeyx<- glht(fit2,linfct = mcp(stage="Tukey"))
  Tukeyy<-summary(Tukeyx)
  results[[BPK_name]] <- list(
    fit=fit2,
    anova= anova_results,
    residuals = residuals,
    normality=normality_test,
    shapiro=shapirox,
    plot_resid=plot_residx,
    plot_qq=plot_qqx,
    barttest=bartx,
    levetest=levex,
    Tukey=Tukeyy)
}
results$BPK275$Tukey

#code for siginificance 
significance_codes <- function(p_value){
  if (p_value < 0.001) {
    return("***")
  } else if (p_value < 0.01) {
    return("**")
  } else if (p_value < 0.05) {
    return("*")
  } else if (p_value < 0.1) {
    return(".")
  } else {
    return("")
  } 
}

#loop for statistical analysis for excell using anova
results_list <- list()
for (BPK_name in names(L4)){
  LX<-L4[[BPK_name]]
  fit2<-lm(median_MitoPo~stage,LX)
  anova_results <- Anova(fit2)
  residuals <- augment(fit2)$.resid
  normality_test <- ols_test_normality(fit2)
  shapirox<-shapiro.test(fit2$residuals)
  plot_residx<- ols_plot_resid_fit(fit2)
  plot_qqx<-ols_plot_resid_qq(fit2)
  bartx<-bartlett.test(median_MitoPo~stage,data=LX)
  levex<-leveneTest(median_MitoPo~stage,data=LX)
  Tukeyx<- glht(fit2,linfct = mcp(stage="Tukey"))
  Tukeyy<-summary(Tukeyx)
  Tukey_p_values <-Tukeyy$test$pvalues
  Tukey_significance <-sapply(Tukey_p_values,significance_codes)
  Tukey_results <-paste0(round(Tukey_p_values, 5), " ",Tukey_significance)
  results_df <- data.frame(
    BPK = BPK_name,
    anovar_p_value= anova_results$`Pr(>F)`[1],
    residuals=paste(round(residuals, 3), collapse = ","),
    shapiro_W= shapirox$statistic,
    shapiro_p_value= shapirox$p.value,
    Bartlett_p_value= bartx$p.value, 
    levene_p_value= levex$`Pr(>F)`[1],
    D4_D2= Tukey_results[1],
    D10_D2= Tukey_results[2],
    PAT2_D2= Tukey_results[3],
    PAT4_D2= Tukey_results[4],
    PAT10_D2= Tukey_results[5],
    D10_D4= Tukey_results[6],
    PAT2_D4= Tukey_results[7],
    PAT4_D4= Tukey_results[8],
    PAT10_D4= Tukey_results[9],
    PAT2_D10= Tukey_results[10]
  )
  results_list[[BPK_name]]<-results_df
}

#exporting the data of the analysis#
final_results<- do.call(rbind, results_list)
write_xlsx(final_results,"02_R_files/statistical_results_Mito_NoTrans_Tukey.xlsx")

##descriptive statistics##
data8<-read.csv("02_R_Files/data_f_Mi.csv")

data8<-data8[data8$count>40,]
data9<-data8%>%
  dplyr::group_by(BPK,stage)%>%
  dplyr::summarise_at(vars(c(per_TrueAlive,median_EGFP,median_MitoPo,)),
                      list(mean=mean,sd=sd,se=std.error))

data9$drug <- ifelse(grepl("Log|Sta",data9$stage),"NoPAT",
                     ifelse(grepl("PAT", data9$stage),"PAT",NA))

data9$stage<-factor(data9$stage, 
                    levels = c("Pro_Log","Pro_Sta_4","Pro_Sta_10",
                               "Pro_PAT 2","Pro_PAT 4","Pro_PAT 10"), 
                    ordered = TRUE)
data9

data9a<-data9[data9$BPK == "BPK026",]
data9b<-data9[data9$BPK == "BPK275",]

##generating the plot BPK026##
plot3a<-ggplot(data = data9a, aes(x = stage, y = median_MitoPo_mean, fill = drug))+
  geom_col(position = "dodge", colour = "black",width = 0.8)+
  geom_errorbar(aes(ymin = median_MitoPo_mean - median_MitoPo_se, 
                    ymax = median_MitoPo_mean + median_MitoPo_se), 
                width = 0.3, position = position_dodge(0.9), color = "black") +
  scale_fill_manual(name = " ", values = c("#525252", "#d9d9d9"),
                    labels = c("No PAT","PAT")) +
  scale_x_discrete(labels = c("Pro_Log" = "NO PAT\nD2", "Pro_Sta_4" = "NO PAT\nD4", 
                              "Pro_Sta_10" = "NO PAT\nD10", "Pro_PAT 2" = "PAT\nD2", 
                              "Pro_PAT 4" = "PAT\nD4", "Pro_PAT 10" = "PAT\nD10")) +
  scale_y_continuous(breaks = seq(0,34000, by =5000),limits = c(0,34000))+
  theme(axis.line = element_line(), panel.background = element_rect(fill = "white"), 
        axis.text.x = element_text(face = "bold", size = 7), 
        axis.text.y = element_text(face = "bold", size = 10),
        axis.title.y = element_text(face="bold",size=10),
        legend.text = element_text(face = "bold", size = 9),
        legend.position = "bottom",
        strip.background = element_rect(fill = "white"),  
        strip.text = element_text(face = "bold", size = 14),
        plot.title = element_text(size = 12, face = "bold", hjust = 0.5)) +
  guides(y = guide_axis(cap = "both", minor.ticks = TRUE)) + 
  labs(title = "BPK026", x = "", y = "Mitochondrial Membrane 
    Potential (MFI)") 

plot3a

plot4a <- plot3a +
  annotate("text", x=2.00, y = 20000, label="*", size=4)+
  annotate("text", x=3.00, y = 17000, label="**", size=4)+
  annotate("text", x=4.00, y = 17000, label="**", size=4)+
  annotate("text", x=5.00, y = 17000, label="***", size=4)+
  annotate("text", x=6.00, y = 3000, label="***", size=4)

plot4a

##generating the plot BPK275##
plot3b<-ggplot(data = data9b, aes(x = stage, y = median_MitoPo_mean, fill = drug))+
  geom_col(position = "dodge", colour = "black",width = 0.8)+
  geom_errorbar(aes(ymin = median_MitoPo_mean - median_MitoPo_se, 
                    ymax = median_MitoPo_mean + median_MitoPo_se), 
                width = 0.3, position = position_dodge(0.9), color = "black") +
  scale_fill_manual(name = " ", values = c("#525252", "#d9d9d9"),
                    labels = c("No PAT","PAT")) +
  scale_x_discrete(labels = c("Pro_Log" = "NO PAT\nD2", "Pro_Sta_4" = "NO PAT\nD4", 
                              "Pro_Sta_10" = "NO PAT\nD10", "Pro_PAT 2" = "PAT\nD2", 
                              "Pro_PAT 4" = "PAT\nD4", "Pro_PAT 10" = "PAT\nD10")) +
  scale_y_continuous(breaks = seq(0,34000, by =5000),limits = c(0,34000))+
  theme(axis.line = element_line(), panel.background = element_rect(fill = "white"), 
        axis.text.x = element_text(face = "bold", size = 7), 
        axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        legend.text = element_text(face = "bold", size = 9),
        legend.position = "bottom",
        strip.background = element_rect(fill = "white"),  
        strip.text = element_text(face = "bold", size = 14),
        plot.title = element_text(size = 12, face = "bold", hjust = 0.5)) +
  guides(y = guide_axis(cap = "both", minor.ticks = TRUE)) + 
  labs(title = "BPK275", x = "", y = "Mitochondrial Membrane 
    Potential (MFI)") 

plot3b

plot4b <- plot3b +
  annotate("segment", x=4.00, xend=6.00 ,y = 15000, yend=15000,size=0.5)+
  annotate("text", x=5.00, y = 16000, label="*", size=4)+
  annotate("segment", x=3.00, xend=4.00 ,y = 16000, yend=16000,size=0.5)+
  annotate("text", x=3.50, y = 16500, label="*", size=4)+
  annotate("segment", x=2.00, xend=4.00 ,y = 18000, yend=18000,size=0.5)+
  annotate("text", x=3.00, y = 19000, label="*", size=4)

plot4b

#combine the plots#
combin_plot1<-ggdraw() +
  draw_plot(plot4a, x = 0, y = 0.2, width = 0.57, height = 0.8) +
  draw_plot(plot4b, x = 0.57, y = 0.2, width = 0.44, height = 0.8)

combin_plot1

#export figures#
ggsave(file="02_R_files/Rplot01_Mito_fig2d.pdf", width = 6, height = 5)
ggsave(file="02_R_files/Rplot01_Mito_fig2.tif", width = 6, height = 5)


####Fig3A_Inh_concen####
##uploading the data##
data0<-read.csv("05_Expe_Inhibitors_B026_B275/03_FCS_Files/2024_BPK_Inh.csv")
head(data0)

##filtering the data and preparing the data##
L9<-list(data1=data0)
L9b<-list()
for (data_name in names(L9)) {
  LX<-L9[[data_name]][-((nrow(L9[[data_name]])-1):nrow(L9[[data_name]])),]
  LX<-LX[,-5]
  colnames(LX)<-c("name","per_TrueAlive","median_EGFP","count")
  LX$BPK<- ifelse(grepl("B026|BPK026",LX$name),"BPK026",
                  ifelse(grepl("B275|BPK275", LX$name),"BPK275",NA))
  LX$repli <- ifelse(grepl("DC1|DD1|DE1|DF1|DG1|DH1|DC4|DD4|DE4|DF4|DG4|DH4", LX$name), "R1",
                     ifelse(grepl("DC7|DD7|DE7|DF7|DG7|DH7", LX$name), "R1",
                            ifelse(grepl("DC10|DD10|DE10|DF10|DG10|DH10", LX$name), "R1", 
                                   ifelse(grepl("DC13|DD13|DE13|DF13|DG13|DH13", LX$name),"R1","UnKnown"))))
  LX$repli <- ifelse(grepl("DC2|DD2|DE2|DF2|DG2|DH2|DC5|DD5|DE5|DF5|DG5|DH5", LX$name), "R2",
                     ifelse(grepl("DC8|DD8|DE8|DF8|DG8|DH8", LX$name), "R2",
                            ifelse(grepl("DC11|DD11|DE11|DF11|DG11|DH11", LX$name), "R2", 
                                   ifelse(grepl("DC14|DD14|DE14|DF14|DG14|DH14", LX$name),"R2",LX$repli))))
  LX$repli <- ifelse(grepl("DC3|DD3|DE3|DF3|DG3|DH3|DC6|DD6|DE6|DF6|DG6|DH6", LX$name), "R3",
                     ifelse(grepl("DC9|DD9|DE9|DF9|DG9|DH9", LX$name), "R3",
                            ifelse(grepl("DC12|DD12|DE12|DF12|DG12|DH12", LX$name), "R3", 
                                   ifelse(grepl("DC15|DD15|DE15|DF15|DG15|DH15", LX$name),"R3",LX$repli))))
  LX$condition<-ifelse(grepl("DC1_|DC2_|DC3_|DD1_|DD2_|DD3_",LX$name),"C01",
                       ifelse(grepl("DC4_|DC5_|DC6_|DD4_|DD5_|DD6_",LX$name), "C02",
                              ifelse(grepl("DC7|DC8|DC9|DD7|DD8|DD9",LX$name), "C03",
                                     ifelse(grepl("DC10|DC11|DC12|DD10|DD11|DD12",LX$name),"C04",
                                            ifelse(grepl("DC13|DC14|DC15|DD13|DD14|DD15",LX$name),"C05","Unknown Condition")))))
  LX$condition<-ifelse(grepl("DE1_|DE2_|DE3_|DF1_|DF2_|DF3_",LX$name),"C01",
                       ifelse(grepl("DE4_|DE5_|DE6_|DF4_|DF5_|DF6_",LX$name), "C06",
                              ifelse(grepl("DE7|DE8|DE9|DF7|DF8|DF9",LX$name), "C07",
                                     ifelse(grepl("DE10|DE11|DE12|DF10|DF11|DF12",LX$name),"C08",
                                            ifelse(grepl("DE13|DE14|DE15|DF13|DF14|DF15",LX$name),"C09",LX$condition)))))
  LX$condition<-ifelse(grepl("DG1_|DG2_|DG3_|DH1_|DH2_|DH3_",LX$name),"C01",
                       ifelse(grepl("DG4|DG5|DG6|DH4|DH5|DH6",LX$name), "C10",
                              ifelse(grepl("DG7|DG8|DG9|DH7|DH8|DH9",LX$name), "C11",
                                     ifelse(grepl("DG10|DG11|DG12|DH10|DH11|DH12",LX$name),"C12",
                                            ifelse(grepl("DG13|DG14|DG15|DH13|DH14|DH15",LX$name),"C13",LX$condition)))))
  LX$expe <- ifelse(grepl("614_|617_", LX$name), "expe1",
                    ifelse(grepl("619_|620_", LX$name), "expe2",
                           ifelse(grepl("621_|625_", LX$name), "expe3", NA)))
  L9b[[data_name]]<-LX
}

data9<-do.call(rbind, L9b)

head(data9)

#Preparing the data and calculating the ratio#
L9c<-split(data9,data9$expe)
L10<-L9c
L10b<-list()
for (expe_name in names(L10)) {
  LX <- L10[[expe_name]]
  LX<-split(LX,LX$condition)
  if("C01" %in% names(LX)) {
    C1_data<- LX[["C01"]]
    for (cond_name in names(LX)) {
      if (cond_name == "C01") {
        LX[[cond_name]]$ratio<-LX[[cond_name]]$per_TrueAlive/C1_data$per_TrueAlive
      } else {
        LX[[cond_name]]$ratio<-LX[[cond_name]]$per_TrueAlive/C1_data$per_TrueAlive
      }
    } 
    L10b[[expe_name]] <- do.call(rbind,LX)
    L10b[[expe_name]]$inh <-ifelse(grepl("C01",L10b[[expe_name]]$condition),"PAT",
                                   ifelse(grepl("C02|C03",L10b[[expe_name]]$condition),"AMPHO",
                                          ifelse(grepl("C04|C05",L10b[[expe_name]]$condition),"MILTE",
                                                 ifelse(grepl("C06|C07",L10b[[expe_name]]$condition),"FCCP",
                                                        ifelse(grepl("C08|C09",L10b[[expe_name]]$condition),"BORTE",
                                                               ifelse(grepl("C10|C11", L10b[[expe_name]]$condition),"PARO", NA))))))
  }
}

data10<-do.call(rbind,L10b)
head(data10)

write.csv(data10,"02_R_Files/data_f_Inh_1.csv")

##Preparing the data for plotting and data analysis##
data10<-read.csv("02_R_Files/data_f_Inh_1.csv")
data10c<-data10[data10$condition == "C01" & data10$expe == "expe1",]
data10d<-data10[!data10$condition == "C01",]
data10e<- rbind(data10c, data10d)

write.csv(data10e,"02_R_Files/data_f_Inh_2.csv")

##Statistical Analysis##
#loading the data#
d2<-read.csv("02_R_Files/data_f_Inh_2.csv")
head(data10)
#be sure that data is numeric and factoring data 
d2$ratio<-as.numeric(d2$ratio)
d2$BPK<-factor(d2$BPK, levels = c("BPK026","BPK275"))
d2$condition<-factor(d2$condition, 
                     levels = c("C01","C02","C03","C04","C05",
                                "C06","C07","C08","C09","C10","C11"))

#statistical analysis
L4<-split(d2,d2$BPK)
results <- list()
for (BPK_name in names(L4)){
  LX<-L4[[BPK_name]]
  fit2<-lm(ratio~condition,LX)
  anova_results <- Anova(fit2)
  residuals <- augment(fit2)$.resid
  normality_test <- ols_test_normality(fit2)
  shapirox<-shapiro.test(fit2$residuals)
  plot_residx<- ols_plot_resid_fit(fit2)
  plot_qqx<-ols_plot_resid_qq(fit2)
  bartx<-bartlett.test(ratio~condition,data=LX)
  levex<-leveneTest(ratio~condition,data=LX)
  Tukeyx<- glht(fit2,linfct = mcp(condition="Tukey"))
  Tukeyy<-summary(Tukeyx)
  results[[BPK_name]] <- list(
    fit=fit2,
    anova= anova_results,
    residuals = residuals,
    normality=normality_test,
    shapiro=shapirox,
    plot_resid=plot_residx,
    plot_qq=plot_qqx,
    barttest=bartx,
    levetest=levex,
    Tukey=Tukeyy)
}
results$BPK275$Tukey

#code for siginificance 
significance_codes <- function(p_value){
  if (p_value < 0.001) {
    return("***")
  } else if (p_value < 0.01) {
    return("**")
  } else if (p_value < 0.05) {
    return("*")
  } else if (p_value < 0.1) {
    return(".")
  } else {
    return("")
  } 
}

#loop for statistical analysis for excell using anova
results_list <- list()
for (BPK_name in names(L4)){
  LX<-L4[[BPK_name]]
  fit2<-lm(ratio~condition,LX)
  anova_results <- Anova(fit2)
  residuals <- augment(fit2)$.resid
  normality_test <- ols_test_normality(fit2)
  shapirox<-shapiro.test(fit2$residuals)
  plot_residx<- ols_plot_resid_fit(fit2)
  plot_qqx<-ols_plot_resid_qq(fit2)
  bartx<-bartlett.test(ratio~condition,data=LX)
  levex<-leveneTest(ratio~condition,data=LX)
  Tukeyx<- glht(fit2,linfct = mcp(condition="Tukey"))
  Tukeyy<-summary(Tukeyx)
  Tukey_p_values <-Tukeyy$test$pvalues
  Tukey_significance <-sapply(Tukey_p_values,significance_codes)
  Tukey_results <-paste0(round(Tukey_p_values, 5), " ",Tukey_significance)
  results_df <- data.frame(
    BPK = BPK_name,
    anovar_p_value= anova_results$`Pr(>F)`[1],
    residuals=paste(round(residuals, 3), collapse = ","),
    shapiro_W= shapirox$statistic,
    shapiro_p_value= shapirox$p.value,
    Bartlett_p_value= bartx$p.value, 
    levene_p_value= levex$`Pr(>F)`[1],
    C02_C01= Tukey_results[1],
    C03_C01= Tukey_results[2],
    C04_C01= Tukey_results[3],
    C05_C01= Tukey_results[4],
    C06_C01= Tukey_results[5],
    C07_C01= Tukey_results[6],
    C08_C01= Tukey_results[7],
    C09_C01= Tukey_results[8],
    C10_C01= Tukey_results[9],
    C11_C01= Tukey_results[10]
  )
  results_list[[BPK_name]]<-results_df
}

#exporting the data of the analysis
final_results<- do.call(rbind, results_list)
write_xlsx(final_results,"statistical_results_Inh_Tukey.xlsx")

##Generating the BOXplot##
data11<-read.csv("02_R_Files/data_f_Inh_2.csv")
head(data11)

#factorizing the data
data11$inh<-factor(data11$inh, levels = c("PAT","AMPHO","MILTE","FCCP","BORTE","PARO"), ordered = TRUE)
data11a<-data11[data11$BPK=="BPK026",]
data11b<-data11[data11$BPK=="BPK275",]

##Generating the BOXplot BPK026##
plot1a<-ggplot(data11a, aes(x=condition, y=ratio))+
  geom_boxplot(aes(colour=inh), fill="white")+
  geom_point(size=2, shape = 21, color="#525252")+
  scale_color_manual(values = c("#525252","#525252","#525252","#525252","#525252","#525252"))+
  scale_x_discrete(labels=c("C01"="C",
                            "C02"="Co",
                            "C03"="S",
                            "C04"="Co",
                            "C05"="S",
                            "C06"="Co",
                            "C07"="S",
                            "C08"="Co",
                            "C09"="S",
                            "C10"="Co",
                            "C11"="S"))+
  scale_y_continuous(limits = c(0,2.5), breaks = seq(-2,2.5, by=0.5))+
  theme(axis.line = element_line(),
        axis.text.x = element_text(face="bold", size=9),
        axis.text.y = element_text(face="bold", size=9),
        strip.background = element_rect(fill="white"),
        strip.text = element_text(face = "bold", size = 14),
        axis.title.y = element_text(face="bold", size = 9),
        panel.background = element_rect(fill="white"),
        plot.title = element_text(size = 14, face = "bold", hjust=0.5))+
  guides(y=guide_axis(cap = "both", minor.ticks = TRUE))+
  labs(title="BPK026", x=" ",y="Ratio of % viable parasite
       (Inhibitor combination / PAT control)")+
  theme(legend.position = "none")+
  geom_vline(xintercept = 1.5, linetype="dashed", color="grey")+
  geom_vline(xintercept = 3.5, linetype="dashed", color="grey")+
  geom_vline(xintercept = 5.5, linetype="dashed", color="grey")+
  geom_vline(xintercept = 7.5, linetype="dashed", color="grey")+
  geom_vline(xintercept = 9.5, linetype="dashed", color="grey")

plot1a

plot2a <- plot1a +
  annotate("text", x=2.00, y = 0.15, label="**", size=5)+
  annotate("text", x=3.00, y = 1.3, label=" ", size=3)+
  annotate("text", x=4.00, y = 0.5, label="*", size=5)+
  annotate("text", x=5.00, y = 1.0, label=" ", size=3)+
  annotate("text", x=6.00, y = 1.4, label=" ", size=3)+
  annotate("text", x=7.00, y = 1.2, label=" ", size=3)+
  annotate("text", x=8.00, y = 0.4, label="**", size=5)+
  annotate("text", x=9.00, y = 1.1, label=" ", size=3)+
  annotate("text", x=10.00, y = 1.9, label=" ", size=3)+
  annotate("text", x=11.00, y = 1.3, label=" ", size=3)

plot2a

##Generating the BOXplot BPK275##
plot1b<-ggplot(data11b, aes(x=condition, y=ratio))+
  geom_boxplot(aes(colour=inh), fill="white")+
  geom_point(size=2, shape = 21, color="#525252")+
  scale_color_manual(values = c("#525252","#525252","#525252","#525252","#525252","#525252"))+
  scale_x_discrete(labels=c("C01"="C",
                            "C02"="Co",
                            "C03"="S",
                            "C04"="Co",
                            "C05"="S",
                            "C06"="Co",
                            "C07"="S",
                            "C08"="Co",
                            "C09"="S",
                            "C10"="Co",
                            "C11"="S"))+
  scale_y_continuous(limits = c(0,2.5), breaks = seq(-2,2.5, by=0.5))+
  theme(axis.line = element_line(),
        axis.text.x = element_text(face="bold", size=9),
        axis.text.y = element_text(face="bold", size=9),
        strip.background = element_rect(fill="white"),
        strip.text = element_text(face = "bold", size = 14),
        axis.title.y = element_blank(),
        panel.background = element_rect(fill="white"),
        plot.title = element_text(size = 14, face = "bold", hjust=0.5))+
  guides(y=guide_axis(cap = "both", minor.ticks = TRUE))+
  labs(title="BPK275", x=" ",y="Ratio of % viable parasite
       (Inhibitor combination / PAT control)")+
  theme(legend.position = "none")+
  geom_vline(xintercept = 1.5, linetype="dashed", color="grey")+
  geom_vline(xintercept = 3.5, linetype="dashed", color="grey")+
  geom_vline(xintercept = 5.5, linetype="dashed", color="grey")+
  geom_vline(xintercept = 7.5, linetype="dashed", color="grey")+
  geom_vline(xintercept = 9.5, linetype="dashed", color="grey")

plot1b

plot2b <- plot1b +
  annotate("text", x=2.00, y = 0.15, label="***", size=5)+
  annotate("text", x=3.00, y = 0.4, label="***", size=5)+
  annotate("text", x=4.00, y = 1.3, label=" ", size=3)+
  annotate("text", x=5.00, y = 1.2, label=" ", size=3)+
  annotate("text", x=6.00, y = 0.3, label="***", size=5)+
  annotate("text", x=7.00, y = 1.0, label=" ", size=5)+
  annotate("text", x=8.00, y = 0.2, label="***", size=5)+
  annotate("text", x=9.00, y = 0.9, label="***", size=5)+
  annotate("text", x=10.00, y = 2.4, label="***", size=5)+
  annotate("text", x=11.00, y = 1.6, label="***", size=5)

plot2b

#combine the plots
combin_plot1<-ggdraw() +
  draw_plot(plot2a, x = 0, y = 0.2, width = 0.50, height = 0.8) +
  draw_plot(plot2b, x = 0.50, y = 0.2, width = 0.47, height = 0.8)

combin_plot1

#export figures
ggsave(file="02_R_files/Rplot01_Inh_figure3.pdf", width = 8, height = 4)
ggsave(file="02_R_files/Rplot01_Inh_BPK_.tif", width = 8, height = 4)

####Fig3B_Bortezomib####
##upploading the data##
data0<-read.csv("06_Expe_Bortezomib_B026_B275/03_FCS_Files/2024_BPK_BOR.csv")
head(data0)

#filtering the data
L9<-list(data1=data0)
L9b<-list()
for (data_name in names(L9)) {
  LX<-L9[[data_name]][-((nrow(L9[[data_name]])-1):nrow(L9[[data_name]])),]
  LX<-LX[,-5]
  colnames(LX)<-c("name","per_TrueAlive","median_EGFP","count")
  LX$BPK<- ifelse(grepl("B026|BPK026",LX$name),"BPK026",
                  ifelse(grepl("B275|BPK275", LX$name),"BPK275",NA))
  LX$condition<-ifelse(grepl("T1_|T2_|T3_|T4_",LX$name),"C1",
                       ifelse(grepl("T5_|T6_|T7_|T8_",LX$name), "C2",
                              ifelse(grepl("T9_|T10_|T11_|T12_",LX$name), "C3",
                                     ifelse(grepl("T13_|T14_|T15_|T16_",LX$name),"C4",
                                            ifelse(grepl("T17_|T18_|T19_|T20_",LX$name),"C5",
                                                   ifelse(grepl("T21_|T22_|T23_|T24_",LX$name),"C6",
                                                          ifelse(grepl("T25_|T26_|T27_|T28_",LX$name),"C7",NA)))))))
  L9b[[data_name]]<-LX
}

data9<-do.call(rbind, L9b)
head(data9)

#Preparing the data and calculating the ratio#
L10<-L9b
L10b<-list()
for (data_name in names(L10)) {
  LX <- L10[[data_name]]
  LX<-split(LX,LX$condition)
  for (cond_name in names(LX)) {
    if (cond_name == "C1") {
      LX[[cond_name]]$ratio<-LX[[cond_name]]$per_TrueAlive/LX[[cond_name]]$per_TrueAlive
    } else {
      LX[[cond_name]]$ratio<-LX[[cond_name]]$per_TrueAlive/LX$C1$per_TrueAlive
    }
  }
  L10b[[data_name]] <- do.call(rbind,LX)
  L10b[[data_name]]$inh <-ifelse(grepl("C1",L10b[[data_name]]$condition),"PAT",
                                 ifelse(grepl("C2|C3",L10b[[data_name]]$condition),"BORT 1",
                                        ifelse(grepl("C4|C5",L10b[[data_name]]$condition),"BORT 2",
                                               ifelse(grepl("C6|C7",L10b[[data_name]]$condition),"BORT 5",NA))))
}


data10<-do.call(rbind,L10b)
head(data10)

#saving the data
write.csv(data10,"02_R_Files/data_f_BOR_1.csv")

##Statistical Analysis##
#loading the data#
d2<-read.csv("02_R_Files/data_f_BOR_1.csv")

#be sure that data is numeric and factoring data 
d2$ratio<-as.numeric(d2$ratio)
d2$BPK<-factor(d2$BPK, levels = c("BPK026","BPK275"))
d2$condition<-factor(d2$condition, 
                     levels = c("C1","C2","C3","C4","C5","C6","C7"))


##statistical analysis PART 2 T-TES##
#a t-test analysis for some samples was done
d2<-read.csv("02_R_Files/data_f_BOR_1.csv")

#be sure that data is numeric and factoring data 
d2$ratio<-as.numeric(d2$ratio)
d2$BPK<-factor(d2$BPK, levels = c("BPK026","BPK275"))
d2$condition<-factor(d2$condition, 
                     levels = c("C1","C2","C3","C4","C5","C6","C7"))


head(d2)
#analysis
#for shapiro wilk
L7<-list(data2=d2[d2$condition == "C2",], 
         data3=d2[d2$condition == "C3",],
         data4=d2[d2$condition == "C4",],
         data5=d2[d2$condition == "C5",],
         data6=d2[d2$condition == "C6",],
         data7=d2[d2$condition == "C7",])

L7b<-list()
for (data_name in names(L7)) {
  LX<-L7[[data_name]]
  LY<-split(LX,LX$BPK)
  L7b[[data_name]]<-lapply(LY, function(x) 
    shapiro.test(x$ratio)$p.value)
}
L7b

#for equal variances#at the end was not variance.
L8<-split(d2, d2$BPK)
L8b<-list()
for (BPK_name in names(L8)) {
  LX<-L8[[BPK_name]]
  L8b[[BPK_name]]<-var.test(LX[LX$condition=="C2","ratio"],
                            LX[LX$condition=="C3","ratio"],
                            LX[LX$condition=="C4","ratio"],
                            LX[LX$condition=="C5","ratio"],
                            LX[LX$condition=="C6","ratio"],
                            LX[LX$condition=="C7","ratio"],)
}
L8b

#for t.test
L9<-list(data1=d2[d2$condition == "C1",], 
         data2=d2[d2$condition == "C2",])

L9<-list(data1=d2)
L9b<-list()
for (data_name in names(L9)) {
  LX<-L9[[data_name]]
  LY<-split(LX,LX$BPK)
  L9b[[data_name]]<-lapply(LY, function(x) t.test(ratio ~ condition, var.equal = FALSE, data=x)$p.value)
}

L9b

d2a<-d2[d2$BPK == "BPK026",]
dataT<-d2a[d2a$condition %in% c("C1","C2"),]
t.test(ratio ~ condition, var.equal = FALSE, data=dataT)


##Generating the BOXplot##
#loading the data
data10<-read.csv("02_R_Files/data_f_BOR_1.csv")
head(data10)
#factorizing the data
data10$inh<-factor(data10$inh, levels = c("PAT","BORT 1","BORT 2","BORT 5"), ordered = TRUE)
data11a<-data10[data10$BPK=="BPK026",]
data11b<-data10[data10$BPK=="BPK275",]

##Generating the BOXplot BPK026##
plot11a<-ggplot(data11a, aes(x=condition, y=ratio))+
  geom_boxplot(aes(colour=inh), fill="white")+
  geom_point(size=2, shape = 21, color="#525252")+
  scale_color_manual(values = c("#525252","#525252","#525252","#525252"))+
  scale_x_discrete(labels=c("C1"="C",
                            "C2"="Co",
                            "C3"="S",
                            "C4"="Co",
                            "C5"="S",
                            "C6"="Co",
                            "C7"="S"))+
  scale_y_continuous(limits = c(0,1.5))+
  theme(axis.line = element_line(),
        axis.text.x = element_text(face="bold", size=9),
        axis.text.y = element_text(face="bold", size=9),
        strip.background = element_rect(fill="white",color = "black"),
        strip.text = element_text(face = "bold", size = 14),
        axis.title.y = element_text(face="bold", size = 10),
        panel.background = element_rect(fill="white"),
        plot.title = element_text(size = 14, face = "bold", hjust=0.5))+
  guides(y=guide_axis(cap = "both", minor.ticks = TRUE))+
  labs(title="BPK026", x=" ",y="Ratio of % viable parasite
       (Inhibitor combination / PAT control)")+
  theme(legend.position = "none")+
  geom_vline(xintercept = 1.5, linetype="dashed", color="grey")+
  geom_vline(xintercept = 3.5, linetype="dashed", color="grey")+
  geom_vline(xintercept = 5.5, linetype="dashed", color="grey")

plot11a

plot12a <- plot11a +
  annotate("text", x=2.00, y = 0.5, label="***", size=5)+
  annotate("text", x=3.00, y = 1.3, label=" ", size=3)+
  annotate("text", x=4.00, y = 0.3, label="***", size=5)+
  annotate("text", x=5.00, y = 1.0, label=" ", size=3)+
  annotate("text", x=6.00, y = 0.25, label="***", size=5)+
  annotate("text", x=7.00, y = 1.2, label=" ", size=3)

plot12a

##Generating the BOXplot BPK275##
plot11b<-ggplot(data11b, aes(x=condition, y=ratio))+
  geom_boxplot(aes(colour=inh), fill="white")+
  geom_point(size=2, shape = 21, color="#525252")+
  scale_color_manual(values = c("#525252","#525252","#525252","#525252"))+
  scale_x_discrete(labels=c("C1"="C",
                            "C2"="Co",
                            "C3"="S",
                            "C4"="Co",
                            "C5"="S",
                            "C6"="Co",
                            "C7"="S"))+
  scale_y_continuous(limits = c(0,1.5))+
  theme(axis.line = element_line(),
        axis.text.x = element_text(face="bold", size=9),
        axis.text.y = element_text(face="bold", size=9),
        strip.background = element_rect(fill="white",color = "black"),
        strip.text = element_text(face = "bold", size = 14),
        axis.title.y = element_blank(),
        panel.background = element_rect(fill="white"),
        plot.title = element_text(size = 14, face = "bold", hjust=0.5))+
  guides(y=guide_axis(cap = "both", minor.ticks = TRUE))+
  labs(title="BPK275", x=" ",y="Ratio of % viable parasite
       (Inhibitor combination / PAT control)")+
  theme(legend.position = "none")+
  geom_vline(xintercept = 1.5, linetype="dashed", color="grey")+
  geom_vline(xintercept = 3.5, linetype="dashed", color="grey")+
  geom_vline(xintercept = 5.5, linetype="dashed", color="grey")

plot11b

plot12b <- plot11b +
  annotate("text", x=2.00, y = 0.5, label="***", size=5)+
  annotate("text", x=3.00, y = 1.3, label=" ", size=3)+
  annotate("text", x=4.00, y = 0.3, label="***", size=5)+
  annotate("text", x=5.00, y = 1.0, label=" ", size=3)+
  annotate("text", x=6.00, y = 0.2, label="***", size=5)+
  annotate("text", x=7.00, y = 1.2, label=" ", size=3)

#combine the plots
combin_plot1<-ggdraw() +
  draw_plot(plot12a, x = 0, y = 0.2, width = 0.50, height = 0.8) +
  draw_plot(plot11b, x = 0.50, y = 0.2, width = 0.45, height = 0.8)

combin_plot1

#export figures
ggsave(file="02_R_files/Rplot01_BOR_figure3b.pdf", width = 7, height = 4)
ggsave(file="02_R_files/Rplot01_BOR_BPK.tif", width = 6, height = 4)



####Table2_IC50calulation####
setwd("~/Project_Managment/01_PQuiescence/02_paper/04_Results/02_IC50_Data_MRPA")

data1<-read.csv("2024_IC50_values_BPKL.csv")
head(data1)

data2<-data1%>%
  dplyr::group_by(BPK,drug)%>%
  dplyr::summarise_at(vars(c(IC)),
                      list(mean=mean,sd=sd,se=std.error))
data2


####Supp_Table4_MRPA####
#MRPA calulation#
setwd("~/Project_Managment/01_PQuiescence/02_paper/04_Results/02_IC50_Data_MRPA")
data1<-read.csv("2024_MRPAvalues.csv")
head(data1)

data2<-data1%>%
  dplyr::group_by(BPK,Inh)%>%
  dplyr::summarise_at(vars(c(MRPA)),
                      list(mean=mean,sd=sd,se=std.error))
data2

####Supp_Fig1_PercentMetabolic####
#make the working directory, upload data#
setwd("~/Project_Managment/01_PQuiescence/02_paper/04_Results/03_ResazurinData")

##preparing the data##
data1 <- read_csv("data_f_0.csv")

data1<- data1[,c(1:4)]
head(data1)

#Simple statistics 
data13<-data1%>%
  dplyr::group_by(BPK,PAT)%>%
  dplyr::summarise_at(vars(c(per)),list(mean=~mean(.x,na.rm=TRUE),
                                        sd=~sd(.x,na.rm=TRUE),
                                        se=~std.error(.x, na.rm=TRUE)))
data13

#factorizing data for plots
data13$ISC<-factor(data13$ISC, levels = c("ISC1","ISC4","Other","ISC6","ISC9","ISC5"), ordered = TRUE)
data13$BPK<-factor(data13$BPK, levels = c("BPK026","BPK031","BPK156","BPK190","BPK087","LdBOB","BPK080","BPK085",
                                          "BPK282","BPK294","BPK275"))

#plot foe the viability for day 0 till day 10
#plot 11
p11<-ggplot(data13,aes(PAT,mean))+
  geom_point(shape=21,size=2.5,colour="black")+
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se),width=0.1)+
  scale_x_continuous(breaks=seq(0,3.5,by=0.5), limits = c(0,3.5))+
  scale_y_continuous(breaks =seq(0,110, by=20),limits = c(-10,110))+  
  geom_line(size=0.4,color="#737373")+
  facet_wrap(~BPK, scales = "free", nrow = 4)+
  theme(axis.line=element_line(), strip.text = element_text(size=12, face ="bold"),
        strip.background = element_blank(),
        axis.text.x = element_text(face="bold", size=10),
        axis.text.y = element_text(face="bold", size=10),
        legend.title = element_text(size=9, face="bold"),
        legend.position = "right",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill="white"),
        axis.title = element_text(size = 11, face = "bold"),
        plot.title = element_text(size = 16, face = "bold"))+
  guides(x=guide_axis(cap = "both",minor.ticks = TRUE), 
         y=guide_axis(cap="both",minor.ticks = TRUE),
         fill=guide_legend(nrow = 11))+
  labs(title= " ",y="% of active cells", x="Log [SbIII] (μM)")

p11

#saving the plot#
ggsave(file="Rplot01_PerceMeta.pdf", width = 8.2, height = 9.7)
ggsave(file="Rplot01_PerceMeta.tif", width = 6.5, height = 5)

####Supp_Table5_ConInh_IC50####
#calculation#
setwd("~/Project_Managment/01_PQuiescence/02_paper/04_Results/02_IC50_Data_MRPA")
data1<-read.csv("2024_ConInh_IC50.csv")
head(data1)

data2<-data1%>%
  dplyr::group_by(BPK,Inh)%>%
  dplyr::summarise_at(vars(c(IC50)),
                      list(mean=mean,sd=sd,se=std.error))
data2

