# Loading required packages
library(ggplot2)
library(sf)
library(RColorBrewer)
library(rstudioapi)
library(stringi)
library(dplyr)
library(ggpubr)
library(gridExtra)


#set path to current folder
current_path = rstudioapi::getActiveDocumentContext()$path 
setwd(dirname(current_path ))
#Map of Russian Federation was obtained from GADM:
#https://gadm.org/download_country_v3.html
rus.map.init <- readRDS("data/gadm36_RUS_1_sf.rds")
ukr.map <- readRDS("data/gadm36_UKR_1_sf.rds")
##subset Crimea
Crimea.map <- ukr.map[ukr.map$NAME_1 %in% c("Crimea","Sevastopol'"),]
Crimea.map$NL_NAME_1<-c("Республика Крым","г. Севастополь")
####
rus.map<-rbind(rus.map.init, Crimea.map)
rus.map$NL_NAME_1[rus.map$NL_NAME_1 == "Пермская край"] <- "Пермский край"
rus.map$NL_NAME_1[rus.map$NL_NAME_1 == "Камчатская край"] <- "Камчатский край"
rus.map$NL_NAME_1[rus.map$NL_NAME_1 == "Республика Чечено-Ингушская"] <- "Чеченская республика"
rus.map$NL_NAME_1[rus.map$NL_NAME_1 == "Респу́блика Ингуше́тия"] <- "Республика Ингушетия"
rus.map$NL_NAME_1[rus.map$NL_NAME_1 == "Санкт-Петербург (горсовет)"] <-"г. Санкт-Петербург"
rus.map$NL_NAME_1[is.na(rus.map$NL_NAME_1)] <-"г. Москва"
###
#general graphical values
dimh<-18
dimw<-32
path<-"Figures"
#
lwdp = 0.05
projection =3576 #EPSG code
common = theme_classic()+
  theme(axis.line = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.position="right",
        legend.text = element_text(size =10),
        legend.title = element_text(size =10),
        legend.key.size = unit(0.4, "cm"),
        plot.margin = unit(c(0,0,0,0), "lines"))
####
##################
###Draw number of cases
#Fig1a Cases per region
cases.spotcoron<-read.csv("data/20200526_cases_cum.csv", sep = "\t", header = T, stringsAsFactors = F)
cases.spotcoron_sorted<-cases.spotcoron[match(rus.map$NL_NAME_1, cases.spotcoron$NL_Name),]

labels<-as.data.frame(cbind(rus.map$NL_NAME_1,st_centroid(rus.map$geometry)))
colnames(labels)<-c("label","geom")
labels_sf <- labels %>%
  st_as_sf(crs = 4326)
labels_transformed <- labels_sf %>%
  st_transform(crs = projection)
labels_transformed_with_lat_lon <- cbind(labels_transformed, st_coordinates(labels_transformed))
########
#read data on genome numbers
seqs.by.region<-read.csv("data/20200526_number_seqs.csv",
                         stringsAsFactors = FALSE, sep = "\t", header = T)
##aggregate by regions on the map
seqs.by.region.aggr<-aggregate(seqs.by.region$Seqs, by=list(Category=seqs.by.region$NL_NAME), FUN=sum)
labels_seqs<-subset(labels_transformed_with_lat_lon, labels_transformed_with_lat_lon$label %in% seqs.by.region.aggr$Category)
seqs.by.region.aggr_sort<-seqs.by.region.aggr[match(labels_seqs$label, seqs.by.region.aggr$Category),]
#get number to variable
labels_seqs_num<-cbind(labels_seqs,seqs.by.region.aggr_sort$x)
########Plot
cases_genomes<-ggplot() +
  geom_sf(data = rus.map, aes(fill=log10(cases.spotcoron_sorted$cases)),
          lwd = lwdp) +
  scale_fill_gradient(low="#ffeda0", high="#bd0026",
                      name = "Cases per\nregion (log10)")+
  geom_point(data = labels_seqs_num, 
             aes(x = X, y = Y,
                 size = seqs.by.region.aggr_sort$x),
             color = "#252525", alpha=.3) +
  geom_text(data = labels_seqs_num,
            aes(x = X, y = Y, 
                label=seqs.by.region.aggr_sort.x),
            hjust=-0.4, vjust=-0.4,
            size = 4) +
  scale_size_continuous(limits = c(1,137),
                        range=c(2,12),
                        breaks = c(1, 2, 5, 10, 40,130),
                        labels = c(1,2,5,10,40,130))+
  coord_sf(crs=projection) +
  guides(size =FALSE)+
  common+
  theme(legend.box = "vertical")

###########################
#Figure 1b Draw lineages
###########################
#load cluster data
clusters<-read.csv("data/Lineages_and_locations.csv",
                   stringsAsFactors = FALSE, sep = "\t", header = T)
###aggregate
clusters.aggr<-aggregate(Seq_id ~ LineageID + label, data = clusters, FUN = length)
#merge with coords
clusters_coord<-merge(clusters.aggr,labels)
#transform to sf object
clusters_sf <- clusters_coord %>%
  st_as_sf(crs = 4326)
clusters_transformed <- clusters_sf %>%
  st_transform(crs = projection)
clusters_transformed_with_lat_lon <- cbind(clusters_transformed, 
                                           jitter(st_coordinates(clusters_transformed), factor=22))

##Plot
mycolors<-c("#00ffff","#0066ff","#c300ff","#99cc00",
            "#ffd300","#bcbde8","#8f0000","#ff8c00",
            "#ff0000")
names(mycolors)<-c("Lineage3","Lineage4","Lineage2","Lineage1",
                   "Lineage8","Lineage9","Lineage7","Lineage5",
                   "Lineage6")
levels_order<-c("Lineage8", "Lineage6", "Lineage1",
                "Lineage4", "Lineage3", "Lineage7",
                "Lineage5", "Lineage2", "Lineage9")
####
clusters_sorted<-clusters_transformed_with_lat_lon %>% arrange(factor(LineageID, levels = levels_order))
###
cluster_plot<-ggplot() +
  geom_sf(data = rus.map, lwd = .5, fill ="#ffffff", color ='#bdbdbd') +
  geom_path(data = clusters_transformed_with_lat_lon, 
            aes(x = X, y = Y,
                color = LineageID),
            size = 1.2,
            alpha = 0.7,
            show.legend = F) +
  geom_point(data = clusters_sorted, 
             aes(x = X, y = Y,
                 size = Seq_id,
                 fill = LineageID),
             shape =21,
             alpha =0.7) +
  scale_size_continuous(limits = c(1,137),
                        range = c(2,13),
                        breaks = c(1, 2, 5, 10, 40,130),
                        labels = c(1,2,5,10,40,130))+
  scale_color_manual(values=mycolors)+
  scale_fill_manual(values=mycolors)+
  coord_sf(crs=projection) +
  common +
  labs(size = "Samples",fill = " ")+
  guides(size = guide_legend(order = 1),
         fill= guide_legend(override.aes = list(size=4)))
  
############################
##Combine everything and save
Figure1_pre<-ggarrange(cases_genomes,
                       cluster_plot,
                       labels = c("a", "b"),
                       ncol = 1, nrow = 2,
                       align = "v",
                       legend.grob = get_legend(cluster_plot),
                       legend = "right")
                       #common.legend = TRUE)
Figure1<-Figure1_pre +
  annotation_custom(grob = get_legend(cases_genomes),
                    xmin =0.87,
                    ymin =0.8)
                   
ggsave("Fig1ab_Maps.png",plot= Figure1, 
       path="Figures/", width = 25, height = 25, dpi=350, units = "cm")