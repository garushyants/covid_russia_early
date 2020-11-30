library(ggplot2)
library(rstudioapi)
library(Hmisc)
library(grDevices)
library(grid)
library(data.table)
library(ggpubr)


current_path = rstudioapi::getActiveDocumentContext()$path 
setwd(dirname(current_path ))

mycolors<-c("#00ffff","#0066ff","#c300ff","#99cc00",
            "#ffd300","#bcbde8","#8f0000","#ff8c00",
            "#ff0000","#ffffff","#00441b","#ffaaaa")
names(mycolors)<-c("Lineage3","Lineage4","Lineage2","Lineage1",
                   "Lineage8","Lineage9","Lineage7","Lineage5",
                   "Lineage6","Stem","Singleton","SingletonNoStem")
CovMetadata<-read.csv("data/Russian_sequences_metadata_2020_05_26_mod.tsv", header =T, sep ="\t")
#########
#Correlation
SpearmanRho<-cor.test(as.numeric(CovMetadata$Patient.age),
         as.numeric(as.Date(CovMetadata$Collection.date, "%Y-%m-%d")),
         method = "spearman")
RhoS <- format(SpearmanRho$estimate, digits = 2)
SPvalue <- format(SpearmanRho$p.value, digits = 2)
###########
###################################
#Common theme
common = theme_classic()+
  theme(axis.text.x = element_text(size =12), 
        axis.text.y = element_text(size =12),
        axis.title.x = element_text(size =12), 
        axis.title.y = element_text(size =12),
        legend.text = element_text(size =12),
        legend.title = element_blank())

###################################
#Plot correlation
plot<-ggplot(data=CovMetadata, 
       aes(y = as.numeric(Patient.age),
           x = as.Date(Collection.date, "%Y-%m-%d"))) +
  geom_point(aes(fill =ClusterSingletonStem),
             shape = 21,
             size=2,
             #colour = "#000000",
             alpha = 0.8)+
  geom_smooth(method=lm,
              se = T)+
  scale_fill_manual(values = mycolors)+
  ylab("Patient age")+
  scale_x_date(name="Collection date", 
               limits=c(as.Date("2020-03-06", "%Y-%m-%d"),
                        as.Date("2020-04-29", "%Y-%m-%d")),
               date_breaks = "7 days",
               date_labels = "%b %d") +
  common +
  theme(axis.text.x = element_blank(), 
        axis.title.x = element_blank())

##
grob <- grobTree(textGrob(paste0("\U03C1 = ",eval(RhoS),
                                           " (p = ", eval(SPvalue),")"), 
                          x=0.05,y = 0.95,
                          hjust =0,
                          gp=gpar(fontsize=12)))

##################
###Plot lineages
# define the summary function
f <- function(x) {
  r <- quantile(x, probs = c(0.00, 0.25, 0.5, 0.75, 1))
  names(r) <- c("ymin", "lower", "middle", "upper", "ymax")
  r
}
###
#####
#Dates
Dates_additional<-data.frame(
  ClusterSingletonStem = c("Lineage1","Lineage2","Lineage3",
                           "Lineage4","Lineage5","Lineage6",
                           "Lineage7","Lineage8","Lineage9"),
  Date = c("2020-03-20","2020-03-13","2020-03-20",
           "2020-03-24","2020-03-11","2020-04-09",
           "2020-03-20","2020-04-03","2020-03-23"),
  StemDate = c("","","","2020-03-13",
               "","2020-03-13","2020-03-13",
               "2020-03-13","")
)
#########
#Order of y labels
Dates_additional_sorted<-Dates_additional[order(Dates_additional$Date),]
lineages_order<-Dates_additional_sorted$ClusterSingletonStem
lineages_order<-c("Stem","Singleton","SingletonNoStem",lineages_order)
#############
CovMetadataSub<-CovMetadata[CovMetadata$ClusterSingletonStem %like% "Lineage",]
############
##
##Plot intervals and lineages
lines<-ggplot() +
  stat_summary(data=CovMetadata, 
               aes(x = as.Date(Collection.date, "%Y-%m-%d"),
                   y = factor(ClusterSingletonStem,
                              level =lineages_order),
                   fill = ClusterSingletonStem),
               alpha =0.4,fill = "#bdbdbd", width = 0.0,
               fun.data = f, geom="boxplot")+
  geom_point(data=Dates_additional,
             aes(y = factor(ClusterSingletonStem,
                            level =lineages_order), 
                 x=as.Date(StemDate, "%Y-%m-%d")),
             color = "#bdbdbd")+
  geom_segment(data = Dates_additional,
               linetype =2,
               aes(y = factor(ClusterSingletonStem,
                                  level =lineages_order), 
                   x = as.Date(StemDate, "%Y-%m-%d"),
                   yend = ClusterSingletonStem,
                   xend = as.Date(Date, "%Y-%m-%d")),
               color = "#bdbdbd")+
  geom_count(data=CovMetadata,
             aes(x = as.Date(Collection.date, "%Y-%m-%d"),
                 y = factor(ClusterSingletonStem,
                            level =lineages_order),
                 fill = ClusterSingletonStem),
             alpha = 0.8, shape =21)+
  geom_text(data = Dates_additional,
            aes(x = as.Date(Date, "%Y-%m-%d"),
                y = factor(ClusterSingletonStem,
                           level =lineages_order), 
                label=format(as.Date(Date, "%Y-%m-%d"), "%b-%d"),
                color = ClusterSingletonStem),
            hjust=1.1, vjust=-0.4)+
  scale_size(range=c(2,12),
             breaks = c(1, 2, 4, 10))+#max_size = 10)+
  scale_color_manual(values = mycolors)+
  scale_fill_manual(values = mycolors, 
                    breaks = c("Lineage1","Lineage2","Lineage3",
                               "Lineage4","Lineage5","Lineage6",
                               "Lineage7","Lineage8","Lineage9",
                               "SingletonNoStem","Singleton","Stem"),
                    labels = c("Lineage1","Lineage2","Lineage3",
                               "Lineage4","Lineage5","Lineage6",
                               "Lineage7","Lineage8","Lineage9",
                               "Singletons","Stem-derived\nsingletons","Stem clusters"))+
  scale_x_date(name="Collection date", 
               limits=c(as.Date("2020-03-06", "%Y-%m-%d"),
                        as.Date("2020-04-29", "%Y-%m-%d")),
               date_breaks = "7 days",
               date_labels = "%b-%d") +
  common +
  theme(axis.text.y = element_blank(), 
        axis.title.y = element_blank(),
        axis.ticks.y = element_blank()) +
  guides(color = FALSE,
         fill= guide_legend(override.aes = list(size=4)))
#lines
####################################
###Add Russian COVID statistics
RussiaCaseNumber<-read.csv("data/Russia_case_data_by_date.csv", header =T, sep ="\t")
Caseplot<-ggplot(data=RussiaCaseNumber,
                 aes(x = as.Date(Date, "%m/%d/%Y"),
                     y = New_cases))+
  geom_bar(stat="identity", fill = "#08306b",position=position_dodge(),
           width = 1.1) +
  labs(y = "Cases per day") +
  scale_x_date(name="Collection date", 
               limits=c(as.Date("2020-03-06", "%Y-%m-%d"),
                        as.Date("2020-04-29", "%Y-%m-%d")),
               date_breaks = "7 days",
               date_labels = "%b-%d") +
  common +
  theme(axis.text.x = element_blank(), 
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank()) +
  guides(color = FALSE,
         fill= guide_legend(override.aes = list(size=4)))
#Caseplot
  
#############################
#Combine plots in one figure
Figure6<-ggarrange(Caseplot, plot + annotation_custom(grob), lines, 
          labels = c("a", "b", "c"),
          vjust = 0.8,
          hjust = 0,
          ncol = 1, nrow = 3,
          heights = c(0.3,1, 0.8),
          align = "v",
          #common.legend = TRUE, 
          legend="right", 
          legend.grob = get_legend(lines))
#save
ggsave("Fig6abс_AgeVsDate_Introductions_unadjasted.svg",
       path = "Figures",plot= Figure6, 
       width = 18, height = 24, dpi = 350, units = "cm")
