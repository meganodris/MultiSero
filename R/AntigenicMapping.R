#----- Antigenic mapping from cross-reactivity estimates -----#
library(Racmacs)
library(ggplot2)
library(ggnewscale)
library(ggrepel)
library(here)


#--- Read in final model estimates 
mus0 <- read.csv(here('Results', 'DENV1+CHIKV+JEV', 'MeanNeg.csv'))
mus1 <- read.csv(here('Results', 'DENV1+CHIKV+JEV', 'MeanPos.csv'))
mus <- rbind(mus0,mus1)


#--- format data for mapping
mus$med <- mus$med + 1 # shift scale up to avoid negative means
mus <- tidyr::spread(mus[,1:3], key='antigen', value='med')
mus <- mus[!mus$pos=='neg',] # remove negative responses
mus <- mus[mus$pos %in% c('DENV1','CHIKV','JEV'),] # keep only monotypic positives
mus <- t(mus) # transform
colnames(mus) <- mus[1,] # set colnames
mus <- mus[!row.names(mus) %in% c('pos','DENV_ELISA'), ] # remove ELISA results


#--- make titer table (mean responses against each antigen)
tt <- acmap(titer_table=mus)
tt$titer_table_flat

#--- run mapping optimizations
map <- optimizeMap(
  map                     = tt,
  number_of_dimensions    = 2,
  number_of_optimizations = 10000,
  minimum_column_basis    = "none"
)

# extract best fitting map
x <- map$optimizations[[1]]
plot_map_table_distance(map, optimization_number = 1,xlim, ylim,line_of_equality = TRUE)


#--- plot map results
ags <- as.data.frame(x$ag_base_coords) # extract antigen coords
srs <- as.data.frame(x$sr_base_coords) # extract sera coords
ags$antigen <- rownames(mus)
srs$sera <- colnames(mus)

# plot
mapplot <- ggplot(srs, aes(V1,V2))+ theme_bw()+ 
  theme(text=element_text(size=18), axis.text=element_blank(), axis.title=element_blank(), axis.ticks=element_blank())+
  geom_point(alpha=0.6,shape=0,size=4,stroke=2)+ 
  geom_point(data=ags, aes(V1,V2),fill='seagreen2',size=4,shape=21,stroke=0.5)+
  geom_text_repel(data=ags, aes(V1,V2,label=antigen), size=5, vjust=-1.2, min.segment.length=Inf)
mapplot


#--- extract antigen-sera distances from map
dists <- mapDistances(map, optimization_number = 1)
dists <- as.data.frame(dists)
colnames(dists) <- colnames(mus)
dists$antigen <- ags$antigen
dists2 <- tidyr::gather(dists, key='sera', value='distance', 1:3)
dists2$agsera <- paste(dists2$antigen, dists2$sera, sep='-')
dists2$distance <- as.numeric(dists2$distance)


#--- read in genetic distances from phylogenetic tree methods
gendists <- read.csv(here('data', 'flavi_genetic_distances.csv'))
gd <- tidyr::gather(gendists, key='antigen2', value='distance', 2:10)
gd$agag <- paste(gd$antigen, gd$antigen2, sep='-')


#--- combine genetic and antigenic dists

# remove CHIKV
dists2 <- dists2[!dists2$antigen=='CHIKV', ]
dists2 <- dists2[!dists2$sera=='CHIKV', ]

# merge antigenic distance results
for(i in 1:nrow(dists2)) dists2$gendist[i] <- gd$distance[gd$agag==dists2$agsera[i]]


# plot comparison
ggplot(dists2, aes(gendist, distance))+ geom_point()+
  ylim(0,NA)+ geom_smooth(method='lm')+ theme_bw()+ ylab('Antigenic distance')+
  xlab('Genetic distance')+ theme(text=element_text(size=12))+
  geom_text_repel(mapping=aes(gendist, distance,label=agsera), hjust=1, vjust=1.3,
                  size=4, box.padding = unit(0, "lines"))

# correlation in genetic and antigenic distances
cor(dists2$distance, dists2$gendist)


