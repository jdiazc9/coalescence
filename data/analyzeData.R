### Juan Diaz-Colunga
### Feb 16 2021

### ----------------------------------------------------------------------
### This script performs all the analyses and produces all the plots
### in the paper.
###    1. Community composition (pie plots & rank plots of ESV frequency)
###    2. Characterization of isolated species from each community
###    3. Matching isolated species to dominants based on sequencing data
###    4. Fractions of isolated species in pairwise competition
###    5. Similarity of coalesced, resident and invasive communities
### ----------------------------------------------------------------------

rm(list=ls())
source('isolates_abundance_from_ESVs.R')
source('myFunctions.R')
library(matlib)
library(testthat)
library(ggplot2)
library(ggh4x)
library(tidyr)
library(RColorBrewer)

# initialize variables to store plots
myplots <- vector(mode='list')
if(!dir.exists('plots')) dir.create('plots')
display_plots <- TRUE
save_plots <- TRUE

# general options
deconvolute_ESVs <- FALSE # should species abundance be estimated from ESV abundance?


### ----------------------------------------------------------------------
### LOAD AND FORMAT DATA
### ----------------------------------------------------------------------

# load OTU table, taxonomy table and metadata
otus <- read.csv('OTU_table.csv',stringsAsFactors=F)
taxa <- read.csv('taxonomy_silva.csv',stringsAsFactors=F)
metadata <- read.csv('metaData.csv',stringsAsFactors=F)

### merge similar taxa?
# merged <- mergeTaxa(otus,taxa)
# otus <- merged[['otus']]
# taxa <- merged[['taxa']]

# check that OTU table and taxonomy table follow the same order
stopifnot(all(otus$X==taxa$X))

# keep only relevant metadata
metadata <- metadata[,c("Sample","Plate","Row","Column","Experiment",
                        "Comm1","StabilizingCarb1","Comm2","StabilizingCarb2",
                        "Carbon","Rep","Transfer")]

# assign IDs to each sequence
otus <- cbind(seq_id = paste('seq_',1:nrow(otus),sep=''),
              otus)

# normalize OTU abundances (so entries in table represent fractions)
otus[,3:ncol(otus)] <- otus[,3:ncol(otus)]/
  matrix(rep(colSums(otus[,3:ncol(otus)]),nrow(otus)),
         nrow=nrow(otus),
         byrow=T)

# change NAs to 'NA' (character) to ease downstream analysis
metadata[is.na(metadata)] <- 'NA'
otus[is.na(otus)] <- 'NA'
taxa[is.na(taxa)] <- 'NA'

# make a 'well' column in metadata to ease visualization
metadata$Well <- paste(metadata$Row,metadata$Column,sep='')

# names and abbreviations of carbon sources, communities and champions
carbon_sources <- c('Glutamine','Citrate')
carbon_sources_abbr <- c('Gln','Cit')
names(carbon_sources_abbr) <- carbon_sources
community_names <- paste('Community',1:8,sep='')
champion_names <- paste('Champion',1:8,sep='')


### ----------------------------------------------------------------------
### COMMUNITY COMPOSITIONS (PIE PLOTS AND RANK PLOTS)
### ----------------------------------------------------------------------

# initialize vector of community compositions
communities <- vector(mode='list',length=2)
names(communities) <- carbon_sources
for (cs in carbon_sources) {
  communities[[cs]] <- vector(mode='list',length=length(community_names))
  names(communities[[cs]]) <- community_names
  for (comm in community_names) {
    communities[[cs]][[comm]] <- data.frame(matrix(nrow=nrow(otus),ncol=0))
    rownames(communities[[cs]][[comm]]) <- otus$seq_id
  }
}

# map original communities (stabilized for 12 transfers) to communities used for coalescence experiments
original_communities <- vector(mode='list',length=2)
names(original_communities) <- carbon_sources
original_communities[['Glutamine']] <- setNames(paste('Community',c(1,3,5,7,8,10,11,12),sep=''),
                                                community_names)
original_communities[['Citrate']] <- setNames(paste('Community',c(1,2,3,7,8,9,10,11),sep=''),
                                              community_names)

# run through all communities and carbon sources, retrieve community composition from sequencing data
for (cs in carbon_sources) {
  for (comm in community_names) {
    
    if (F) { # original communities (before propagation, coalescence etc.)
      wells <- metadata[metadata$Comm1 == original_communities[[cs]][comm] &
                        metadata$Carbon == cs &
                        (grepl('Gln',metadata$Experiment) | grepl('Cit',metadata$Experiment))
                        ,]
      fraction <- otus[,wells$Sample,drop=FALSE]
      colnames(fraction) <- 'original'
      communities[[cs]][[comm]] <- cbind(communities[[cs]][[comm]],
                                         fraction)
    }
    if (T) { # communities propagated for 7 days in control of plates C1 (glutamine) and C4 (citrate)
      wells <- metadata[metadata$Comm1 == comm &
                        metadata$Comm2 == 'None' &
                        metadata$StabilizingCarb1 == cs &
                        metadata$Carbon == cs &
                        (grepl('C1',metadata$Experiment) | grepl('C4',metadata$Experiment))
                        ,]
      fraction <- otus[,wells$Sample,drop=FALSE]
      colnames(fraction) <- 'XvsX-control'
      communities[[cs]][[comm]] <- cbind(communities[[cs]][[comm]],
                                         fraction)
    }
    if (T) { # communities propagated for 7 days in diagonal of plates C1 (glutamine) and C4 (citrate)
      wells <- metadata[metadata$Comm1 == comm &
                        metadata$Comm2 == comm &
                        metadata$StabilizingCarb1 == cs &
                        metadata$StabilizingCarb2 == cs &
                        metadata$Carbon == cs &
                        (grepl('C1',metadata$Experiment) | grepl('C4',metadata$Experiment))
                        ,]
      fraction <- otus[,wells$Sample,drop=FALSE]
      colnames(fraction) <- 'XvsX-diagonal'
      communities[[cs]][[comm]] <- cbind(communities[[cs]][[comm]],
                                         fraction)
    }
    if (T) { # communities propagated for 7 days in control of plates C2 (glutamine) and C5 (citrate)
      wells <- metadata[metadata$Comm1 == comm &
                        metadata$Comm2 == 'None' &
                        metadata$StabilizingCarb1 == cs &
                        metadata$Carbon == cs &
                        (grepl('C2',metadata$Experiment) | grepl('C5',metadata$Experiment))
                        ,]
      fraction <- otus[,wells$Sample,drop=FALSE]
      colnames(fraction) <- 'XvsM-control'
      communities[[cs]][[comm]] <- cbind(communities[[cs]][[comm]],
                                         fraction)
    }
    if (F) { # communities propagated for 7 days in control of plates C7 (glutamine) and C8 (citrate)
      wells <- metadata[metadata$Comm1 == comm &
                        metadata$Comm2 == 'None' &
                        metadata$StabilizingCarb1 == cs &
                        metadata$Carbon ==cs &
                        (grepl('C7',metadata$Experiment) | grepl('C8',metadata$Experiment))
                        ,]
      fraction <- otus[,wells$Sample,drop=FALSE]
      colnames(fraction) <- 'GlnXvsCitX-control'
      communities[[cs]][[comm]] <- cbind(communities[[cs]][[comm]],
                                         fraction)
    }
    
  }
}

# initialize data frame for pie plots and rank plots
plot_this <- data.frame(carbon_source=character(0),
                        community=character(0),
                        sample=character(0),
                        seq_id=character(0),
                        fraction=numeric(0))

# run through all communities in all carbon sources and structure plotting data frame
for (cs in carbon_sources) {
  for (comm in community_names) {
    attach_this <- gather(communities[[cs]][[comm]],
                          sample,fraction,colnames(communities[[cs]][[comm]]))
    attach_this <- data.frame(carbon_source=cs,
                              community=comm,
                              sample=attach_this$sample,
                              seq_id=otus$seq_id,
                              fraction=attach_this$fraction)
    plot_this <- rbind(plot_this,attach_this)
  }
}

# remove sequences with frequency 0 in all samples from the plotting data frame
remove_these <- sapply(otus$seq_id,
                       function(seq_id) {
                         all(plot_this$fraction[plot_this$seq_id == seq_id] == 0)
                       })
remove_these <- otus$seq_id[remove_these]
plot_this <- plot_this[!(plot_this$seq_id %in% remove_these),]

# ranks
plot_this$rank <- NA
for (cs in carbon_sources) {
  for (comm in community_names) {
    for (sample in unique(plot_this$sample)) {
      n <- plot_this$carbon_source==cs & plot_this$community==comm & plot_this$sample==sample
      plot_this$rank[n] <- rank(1 - plot_this$fraction[n],ties.method='random')
    }
  }
}
plot_this$rank[plot_this$fraction < 1e-4] <- NA

# plot: color palette
set.seed(23) # this gives the nicest colors in the pie plots
pl <- brewer.pal(8, "Set1")
pl <- sample(colorRampPalette(pl)(nrow(otus)))
pl_carbon <- c('#d7191c','#2c7bb6')
names(pl_carbon) <- carbon_sources

# show fractions in pies
plot_this$label <- paste(100*round(plot_this$fraction,digits=2),'%',sep='')
plot_this$label[plot_this$fraction < 0.1] <- ' '

# fraction labels positioning
plot_this$label_pos <- NA
for (cs in carbon_sources) {
  for (comm in community_names) {
    for (sample in unique(plot_this$sample)) {
      n <- plot_this$carbon_source==cs & plot_this$community==comm & plot_this$sample==sample
      plot_this$label_pos[n] <- cumsum(plot_this$fraction[n]) - 0.5*plot_this$fraction[n]
    }
  }
}
plot_this$label_pos <- 1 - plot_this$label_pos

# characters as factors
plot_this$carbon_source <- factor(plot_this$carbon_source,levels=c('Glutamine','Citrate'))
plot_this$community <- factor(plot_this$community,levels=community_names)
plot_this$seq_id <- factor(plot_this$seq_id,levels=unique(plot_this$seq_id))
sample_names <- unique(plot_this$sample)
plot_this$sample <- factor(plot_this$sample,levels=sample_names)

# plot pies
myplots[['community-compostion_pieplots-all']] <- 
  ggplot(data=plot_this,
         aes(x="",y=fraction,fill=seq_id)) +
    geom_bar(stat="identity",width=1,color='black') + 
    geom_text(aes(x='',y=label_pos,label=label),color='black',size=3) +
    facet_nested(community ~ carbon_source + sample,
                 switch='y',
                 labeller=labeller(community=setNames(gsub('Community','Community ',community_names),
                                                      community_names))) +
    coord_polar(theta="y",start=pi/2) +
    theme_void() +
    scale_fill_manual(values=pl) +
    theme(legend.position='none',
          strip.background.x=element_rect(fill='white',
                                          color='black'),
          text=element_text(size=15)) +
    ylim(0,1)

# plot ranks
myplots[['community-compostion_rankplots']] <-
  ggplot(data=plot_this[plot_this$fraction>=1e-4 & plot_this$rank<=20,],
         aes(x=rank,y=fraction,
             group=interaction(carbon_source,community,sample),
             color=carbon_source)) +
    geom_line(size=0.25) +
    scale_y_continuous(trans='log10',
                       name='Relative abundance',
                       limits=10^c(-4,0),
                       breaks=10^seq(-4,0,by=1),
                       labels=c(expression(10^-4),
                                expression(10^-3),
                                expression(10^-2),
                                expression(10^-1),
                                expression(10^0))) +
    scale_x_continuous(name='Rank',
                       limits=c(0,20)) +
    scale_color_manual(values=pl_carbon) +
    theme_bw() +
    theme(panel.grid=element_blank(),
          legend.title=element_blank(),
          legend.position=c(0.6,0.9),
          legend.background=element_rect(fill='transparent'),
          text=element_text(size=15),
          axis.text=element_text(size=15),
          axis.line=element_blank(),
          axis.ticks=element_line(size=0.25),
          panel.border=element_rect(size=0.25))

# display and save plots
if (display_plots) {
  print(myplots[['community-compostion_pieplots-all']])
  print(myplots[['community-compostion_rankplots']])
}
if (save_plots) {
  ggsave(file.path('.','plots','community-compostion_pieplots-all.pdf'),
         plot=myplots[['community-compostion_pieplots-all']],
         device='pdf',
         height=297,
         width=210,
         units='mm',dpi=600)
  ggsave(file.path('.','plots','community-compostion_rankplots.pdf'),
         plot=myplots[['community-compostion_rankplots']],
         device='pdf',
         height=70,
         width=70,
         units='mm',dpi=600)
}

### Pie plots show that community compositions are very conserved across
### biological replicates for the most part. There are a couple of
### exceptions:
###    - Community 1 of Glutamine: the community in the diagonal of the
###      XvsX plate displays a different structure than the same community
###      in the control of XvsX and the control of XvsM. The dominant does
###      seem conserved.
###    - Community 8 of Citrate: again, the community in the diagonal of
###      the XvsX plate shows a slightly different composition than the
###      replicates communities in the other plates. The dominant also
###      appears conserved across replicates, but the subdominant species
###      have notably different relative abundances.
### By analyzing the distance between pairs of replicate communities in the
### ESV space, we can more systematically identify the instances where a
### specific replicate displays significant deviations in composition from
### the others.

# initialize vector of distances
communities_dist <- numeric(0)

# run through communities and carbon sources and get distances between replicate pairs
for (cs in carbon_sources) {
  for (comm in community_names) {
    fraction <- communities[[cs]][[comm]]
    fraction_dist <- as.matrix(dist(t(fraction)))
    communities_dist <- c(communities_dist,
                          fraction_dist[upper.tri(diag(ncol(fraction)))])
  }
}

### In a boxplot of the distances between replicates, anything above quantile(x)[4] + 1.5*IQR(x)
### would be considered an outlier. Therefore, if we find that two biological replicates are
### separated by a larger distance, we will exclude (at least) one of them from downstream
### analysis, since we can understand that they are in fact not converging at the same state.

# get threshold
communities_dist_threshold <- quantile(communities_dist)[4] + 1.5*IQR(communities_dist)

# find outlier communities and exclude them
for (cs in carbon_sources) {
  for (comm in community_names) {
    fraction <- communities[[cs]][[comm]]
    fraction_dist <- as.matrix(dist(t(fraction)))
    diag(fraction_dist) <- Inf # values in the diagonal don't matter and this eases downstream analyisis
    these_are_outliers <- apply(fraction_dist,
                                1,
                                function(x) all(x>communities_dist_threshold))
    communities[[cs]][[comm]] <- communities[[cs]][[comm]][,!these_are_outliers]
  }
}

### Indeed we see that this method identified Community 1 from glutamine in the XvsX diagonal and
### Community 8 from citrate in the XvsX diagonal as outliers, just as we had anticipated from the
### pie plots.


### ----------------------------------------------------------------------
### ISOLATED SPECIES
### ----------------------------------------------------------------------

### Every community was plated and one species per community was isolated.
### The isolate, in principle, should be the dominant of the community
### (i.e. the most abundant species). However, in some cases this dominant
### may not have been correctly identified. Additionally, some communities
### may share the same isolated species.
### First, we will characterize the isolated species in terms of their
### sequenced 165 rRNA gene. Note that a single species may yield more
### than one exact sequence variant (ESV) if not all copies of the 16S
### share the exact same sequence.
### Then, we will see how these isolates map to the sequencing data of the
### full communities. This way we can identify the instances where the
### dominant was incorrecly identified.

# initialize vector of isolates
isolates <- vector(mode='list',length=2)
names(isolates) <- carbon_sources
for (cs in carbon_sources) {
  isolates[[cs]] <- vector(mode='list',length=length(community_names))
  names(isolates[[cs]]) <- community_names
  for (comm in community_names) {
    isolates[[cs]][[comm]] <- data.frame(matrix(nrow=nrow(otus),ncol=0))
    rownames(isolates[[cs]][[comm]]) <- otus$seq_id
  }
}

# for every champion (isolate) and carbon source, get sequencing info
for (cs in carbon_sources) {
  for (champ in champion_names) {
    wells <- metadata[metadata$Comm1 == champ &
                      (metadata$Comm2 == champ | metadata$Comm2 == 'None') &
                      metadata$StabilizingCarb1 == cs &
                      metadata$Carbon == cs
                      ,]
    fraction <- otus[,wells$Sample,drop=FALSE]
    isolates[[cs]][[gsub('Champion','Community',champ)]] <- fraction
  }
}

### Similar to what we did before, we want to quantify the distance in the ESV space
### across samples that contain the same single isolate. This way we can quantify the
### variability that comes from sequencing.

# initialize vector of distances
isolates_dist <- numeric(0)

# run through carbon sources and communities, quantify distance across replicates
for (cs in carbon_sources) {
  for (comm in community_names) {
    fraction <- isolates[[cs]][[comm]]
    fraction_dist <- as.matrix(dist(t(fraction)))
    isolates_dist <- c(isolates_dist,
                       fraction_dist[upper.tri(diag(ncol(fraction)))])
  }
}

# estimate technical error from sequencing
isolates_dist_threshold <- quantile(isolates_dist)[4] + 1.5*IQR(isolates_dist)

### Like we did with the communities, we will eliminate any sample that is off all
### other replicates by a distance that is larger than this threshold.

# run through carbon sources and communities and remove samples that are off
for (cs in carbon_sources) {
  for (comm in community_names) {
    fraction <- isolates[[cs]][[comm]]
    fraction_dist <- as.matrix(dist(t(fraction)))
    diag(fraction_dist) <- Inf # values in the diagonal don't matter and this eases downstream analyisis
    these_are_outliers <- apply(fraction_dist,
                                1,
                                function(x) all(x > isolates_dist_threshold))
    #print(sum(these_are_outliers))
    isolates[[cs]][[comm]] <- isolates[[cs]][[comm]][,!these_are_outliers]
  }
}

### This approach yields few eliminated samples, which is good (it indicates that the
### variability across samples is within the sequencing error, which discards the
### possibility of contaminations etc.)
### Now, we can check whether there are any repeated isolates that were obtained from
### more than one community. If the isolates extracted from two communities are closer
### in the ESV space than the threshold we just quantified, we can interpret that they
### are the same isolate.

# initialize table with ESV composition of isolates
isolates_esv <- data.frame(matrix(nrow=nrow(otus),ncol=0))
rownames(isolates_esv) <- otus$seq_id

# build table of isolates ESV composition (averaging across replicates)
for (cs in carbon_sources) {
  for (comm in community_names) {
    isolates_esv <- cbind(isolates_esv,
                          rowMeans(isolates[[cs]][[comm]]))
    colnames(isolates_esv)[ncol(isolates_esv)] <- paste(cs,comm,sep='-')
  }
}

# check which isolates are duplicated (their distance in ESV space is below the threhsold)
isolates_dist <- as.matrix(dist(t(isolates_esv)))
diag(isolates_dist) <- 0 # values in diagonal are irrelevant and this simplifies downstream analysis
isolates_groups <- vector(mode='list',length=nrow(isolates_dist)) # list of unique isolates
for (i in 1:nrow(isolates_dist)) {
  isolates_groups[[i]] <- names(which(isolates_dist[i,] < isolates_dist_threshold))
}

# remove redundant groups
isolates_groups <- unique(isolates_groups)

# name species
names(isolates_groups) <- paste('isolate_',1:length(isolates_groups),sep='')

# check that no isolates have gone missing
stopifnot(all(sort(unlist(isolates_groups)) == sort(colnames(isolates_dist))))

### At this point, the 'isolates_groups' variable is a list containing the information of how many
### unique species were isolated from the plating of the communities, and which communities shared
### the same isolated species.
### Note that out of 16 communities (8 for glutamine and 8 for citrate), there are only 9 unique
### isolates due to repetitions.
### Now we want to identify the ESVs composition of each of those isolates.

# first, average across all the columns that correspond to a same unique isolate
for (iso in names(isolates_groups)) {
  isolates_esv <- cbind(isolates_esv,
                        rowMeans(isolates_esv[,colnames(isolates_esv) %in% isolates_groups[[iso]],drop=FALSE]))
  colnames(isolates_esv)[ncol(isolates_esv)] <- iso
}
isolates_esv <- isolates_esv[,colnames(isolates_esv) %in% names(isolates_groups)]

# round to 2 digits while maintaining normalization
for (i in 1:ncol(isolates_esv)) {
  isolates_esv[,i] <- smart_round(100*isolates_esv[,i])/100
}

### At this point, it is reasonable to apply a sequence frequency threshold. Looking at
### isolates_esv, we can see situations (isolates 1, 2, 5 and 7) in which 99% of the sequenced
### reads match a unique ESV. It is reasonable to round this up to 100%, i.e. assume that all
### copies of the 16S rRNA gene of the isolated species match the same unique sequence.
### Other cases are trickier:
###
###    - Isolate 4 is consistent with a species that has 7 copies of the 16S rRNA gene,
###      out of which 6 match the sequence seq_1 and the other one matches seq_13. Sequences 1
###      and 13 are very similar. They are both 233bp long and 99.57% similar:
###      > nchar(otus$X[1])
###      > nchar(otus$X[13])
###      > seq_sim(otus$X[1],otus$X[13])
###      which indicates a difference of just 1bp.
###      The ratios seen in the isolates_esv table are compatible with a species that has 6
###      copies of seq_1 and an additional one of seq_13 (6/7 = 0.86, 1/7 = 0.14)
###
###    - Isolate 6 is a similar case, this time with sequences 1 and 23 and ratios 2/3
###      (= 0.67) and 1/3 (= 0.33). Sequences 1 and 23 are also just 1bp away.
###      > nchar(otus$X[1])
###      > nchar(otus$X[23])
###      > seq_sim(otus$X[1],otus$X[23])
###
###    - Isolate 3 is the strangest. This isolate was obtained from Community 3 in glutamine
###      and is not repeated in any other community. 88% of its sequenced reads match to 
###      sequence 4, but the remaining 12% match to different sequences. In addition, if we
###      go back to the raw data and we identify the samples that contained this isolate in
###      monoculture:
###      > cs <- 'Glutamine'
###      > champ <- 'Champion3'
###      > wells <- metadata[metadata$Comm1 == champ &
###                          (metadata$Comm2 == champ | metadata$Comm2 == 'None') &
###                          metadata$StabilizingCarb1 == cs &
###                          metadata$Carbon == cs
###                          ,]
###      > fraction <- otus[,wells$Sample,drop=FALSE]
###      > View(fraction)
###      we observe that one of the samples matches seq_4 at a fraction of 99%, while the
###      others show that roughly 87% of the reads match seq_4 but the remaining reads match
###      sequences 1, 5, 10 and 13. Interestingly, even in the sample with 99% matches to
###      seq_4, there are still (few) matches to the same sequences 1, 5, 10 and 13. While
###      sequences 1 and 13 are similar (1bp difference), all other combinations are more
###      different (3-40bp difference).
###      If we look at the sequencing data from the experiment where this isolate_3 was
###      co-cultured with another isolate, for example, isolate_2:
###      > cs <- 'Glutamine'
###      > champ1 <- 'Champion2'
###      > champ2 <- 'Champion3'
###      > wells <- metadata[((metadata$Comm1 == champ1 & metadata$Comm2 == champ2) |
###                          (metadata$Comm1 == champ2 & metadata$Comm2 == champ1)) &
###                          metadata$StabilizingCarb1 == cs &
###                          metadata$StabilizingCarb2 == cs &
###                          metadata$Carbon == cs
###                          ,]
###      > fraction <- otus[,wells$Sample,drop=FALSE]
###      > View(fraction)
###      we again see that some reads match sequences 1 (shared by both isolate_3 and 
###      isolate_2), 4, 5, 10 and 13 (unique to isolate_3).
###      This could be because:
###         a) Isolate 3 has 14 copies of the 16S rRNA gene, out of which 12 match sequence
###            4 (12/14 = 0.86), another one matches sequence 1 (1/14 = 0.07) and another
###            one matches sequence 10 (1/14 = 0.07). This feels unlikely.
###         b) Isolate 3 has all copies of the 16S rRNA gene matching to sequence 4, but
###            some mutations arised in some individuals of certain samples. This is
###            unlikely given that sequences 1, 5, 10 and 13 are not a result of point,
###            single base mutations of sequence 4.
###         c) Whatever other reason (contamination...)
###      Whatever the true reason, the safest option is to just consider that isolate 3 maps
###      directly and unequivocally to sequence 4. This sequence is not shared by any other
###      isolate, and this choice ensures that, if a sample that is supposed to contain
###      isolate_3 shows no reads mapping to sequences 1, 5, 10 or 13 (but does show reads
###      mapping to sequence 4), the deconvolution process will not result in a measurement
###      of zero abundance of isolate_3.

# clean sequence mapping table
isolates_esv[isolates_esv <= 0.1] <- 0
isolates_esv[isolates_esv[,'isolate_3'] < 0.5,'isolate_3'] <- 0
isolates_esv <- isolates_esv/matrix(rep(colSums(isolates_esv),nrow(isolates_esv)),
                                    nrow=nrow(isolates_esv),
                                    byrow=TRUE)

# round to 2 digits while maintaining normalization
for (i in 1:ncol(isolates_esv)) {
  isolates_esv[,i] <- smart_round(100*isolates_esv[,i])/100
}

# split groups by carbon source
isolates_groups <- list(Glutamine = isolates_groups[grepl('Glutamine',isolates_groups)],
                        Citrate = isolates_groups[grepl('Citrate',isolates_groups)])
for (cs in carbon_sources) {
  isolates_groups[[cs]] <- lapply(isolates_groups[[cs]],
                                  function(x) gsub(paste(cs,'-',sep=''),'',x[grepl(cs,x)]))
}

# structure variable for isolates
isolate_names <- colnames(isolates_esv)
isolates <- list(seq=isolates_esv,
                 groups=isolates_groups)


### ----------------------------------------------------------------------
### MATCH ISOLATES AND COMMUNITY DOMINANTS
### ----------------------------------------------------------------------

### The question here is, are we sure that the isolate obtained from each
### community is, indeed, the dominant (most abundant) species of that
### community?
### To answer it, we will cross the sequencing data of the whole
### communities with the information we just retrieved from the isolates.
### We will calculate the relative abundance of the isolate in the
### community, assuming that all the reads that match the sequence of the
### isolate come, indeed, from the isolate (in the most general case, some
### or all of those reads could come from different species that share the
### same sequence in at least some of its copies of the 16S, so this
### relative abundance should be taken as the maximum potential relative
### abundance of the isolate in the community).
### All of this is assuming that the same number of reads per individual
### are produced when sequencing. This need not be the case, as different
### species may carry different copy numbers of their 16S rRNA gene. So
### in practice, this approach only serves us to detect extreme cases
### (e.g., the isolate is present in the community only in a very minor
### abundance).

# initialize plotting table
plot_this <- data.frame(carbon_source = character(0),
                        community = character(0),
                        sample = character(0),
                        isolate = character(0),
                        max_rel_abundance = character(0),
                        is_isolate = logical(0))

# run through all carbon sources, communities and samples
for (cs in carbon_sources) {
  for (comm in community_names) {
    for (sample in colnames(communities[[cs]][[comm]])) {
      
      # community composition array
      cc <- as.matrix(communities[[cs]][[comm]][,sample,drop=FALSE])
      
      # matrix to map species relative abundance to ESV relative abundance
      Q <- as.matrix(isolates[['seq']])
      
      # keep only ESVs with abundance over 0 in either the community or the matrix
      n <- rowSums(Q)>0 | cc>0
      Q <- Q[n,,drop=FALSE]
      cc <- cc[n,,drop=FALSE]
      
      # for every isolate...
      for (iso in isolate_names) {
        
        # maximum potential abundance of the isolate in the community
        max_rel_abundance <- min(cc/Q[,iso],na.rm=TRUE)
        
        # is this the isolate that was obtained from this community?
        is_isolate <- comm %in% isolates[['groups']][[cs]][[iso]]
        
        # add info to plotting table
        plot_this <- rbind(plot_this,
                           data.frame(carbon_source = cs,
                                      community = comm,
                                      sample = sample,
                                      isolate = iso,
                                      max_rel_abundance = max_rel_abundance,
                                      is_isolate = is_isolate))
        
      }
      
    }
  }
}

# factors
plot_this$carbon_source <- factor(plot_this$carbon_source,levels=carbon_sources)
plot_this$community <- factor(plot_this$community,levels=community_names)
plot_this$sample <- factor(plot_this$sample,levels=sample_names)
plot_this$isolate <- factor(plot_this$isolate,levels=isolate_names)
plot_this$is_isolate <- factor(plot_this$is_isolate,levels=c(TRUE,FALSE))

# colors (for is_isolate = TRUE or FALSE)
pl_iso <- c('#ca0020','#404040')

# plot relative abundances
myplots[['isolates-abundance-in-communities']] <-
  ggplot(data=plot_this,
         aes(x=isolate,y=max_rel_abundance,fill=is_isolate)) +
    geom_bar(stat="identity",width=1,color='black') + 
    facet_nested(community ~ carbon_source + sample,
                 labeller=labeller(community=setNames(gsub('Community','Community ',community_names),
                                                      community_names))) +
    scale_x_discrete(name='Isolate #',
                     breaks=isolate_names,
                     labels=gsub('isolate_','',isolate_names),
                     expand=expansion(add=1)) +
    scale_y_continuous(name='Max. relative abundance\nof isolate in community') +
    theme_bw() +
    scale_fill_manual(values=pl_iso,
                      name='Isolate obtained\nfrom community?',
                      labels=c('Yes','No')) +
    theme(text=element_text(size=15))

# display and/or save plot
if (display_plots) {
  print(myplots[['isolates-abundance-in-communities']])
}
if (save_plots) {
  ggsave(file.path('.','plots','isolates-abundance-in-communities.pdf'),
         plot=myplots[['isolates-abundance-in-communities']],
         device='pdf',
         height=270,
         width=300,
         units='mm',dpi=600)
}

### Our criterion will be so that we mantain maximum consistency with the experimental
### protocol. We need to consider that a) the dominants at the moment when the plating
### and isolation was done (transfer 12) may have decreased in abundance after 7
### additional transfers and b) species may carry different copy number of the 16s rRNA
### gene so ESV abundance may not perfectly match species abundance. Therefore, we will
### consider a species to be the dominant of a community if it was isolated from it and is
### still at a significant abundance after the 7 additional transfers (i.e. is visible in
### the plot). If the isolate is not visible, we will exclude the corresponding community
### from downstream analysis. From examining the plot:
### Glutamine:
###    - isolate_1 is compatible with being the dominant of Community 1 and it was
###      isolated from that community.
###    - isolate_2 is compatible with being the dominant of Community 2 and it was
###      isolated from that community.
###    - isolate_3 is compatible with being the dominant of Community 3, but isolate_2
###      is also compatible depending on the sample. Because isolate_3 was isolated
###      directly from Community 3, we will consider it as the dominant.
###    - isolate_4 is compatible with being the dominant of Community 4 and it was
###      isolated from that community.
###    - None of the isolates is really compatible with being dominant in Community 5
###      (all relative abundances below ~0.2 when pie plots show an ESV with frequency
###      ~60% in the community).
###    - isolate_4 is compatible with being the dominant of Community 6. isolate_2
###      could also be compatible but it is consistently in a slightly lower abundance
###      across all samples. isolate_4 was also isolated from Community 4 (and appears
###      to really be the dominant). These two communities share the same dominant.
###    - isolate_5 is compatible with being the dominant of Community 7 and it was
###      isolated from that community.
###    - isolate_5 is compatible with being the dominant of Community 8, but isolate_2
###      is also compatible depending on the sample. Because isolate_5 was isolated
###      directly from Community 8, we will consider it as the dominant. Communities
###      7 and 8 share the same dominant.
### 
### Citrate:
###    - isolate_2 seems to clearly be the dominant of Community 1, even though
###      isolate_6 was the one obtained from this community. We will consider
###      isolate_2 to be the dominant of this community.
###    - isolate_7 is compatible with being the dominant of Community 2, but isolate_2
###      also is. Because isolate_7 was directly obtained from Community 2, we will
###      consider it the dominant of this community.
###    - isolate_2 is compatible with being the dominant of Community 3 and it was
###      isolated from that community. Communities 1 and 3 share the same dominant.
###    - isolate_2 seems to clearly be the dominant of Community 4, even though
###      isolate_8 was the one obtained from this community. isolate_8 is also at a
###      relatively high fraction. (*)
###    - isolate_2 is compatible with being the dominant of Community 5 and it was
###      isolated from that community. Communities 1, 3 and 5 share the same
###      dominant.
###    - isolate_2 is compatible with being the dominant of Community 6. isolate_9 was
###      isolated from this community but does not appear to be at significant
###      abundance. We will consider isolate_2 to be the dominant of this community.
###      Communities 1, 3, 5 and 6 share the same dominant.
###    - isolate_7 is compatible with being the dominant of Community 7 and it was
###      isolated from that community. Communities 2 and 7 share the same dominant.
###    - isolate_2 is compatible with being the dominant of Community 8. isolate_4
###      could also be compatible but it is consistently in a slightly lower abundance
###      across all samples. Neither of these two isolates were isolated from this
###      community. We will consider isolate_2 to be the dominant. Communities 1, 3,
###      5, 6 and 8 share the same dominant.

# dominants
dominants <- vector(mode='list',length=2)
names(dominants) <- carbon_sources
dominants[['Glutamine']] <- list(isolate_1 = 'Community1',
                                 isolate_2 = 'Community2',
                                 isolate_3 = 'Community3',
                                 isolate_4 = c('Community4','Community6'),
                                 isolate_5 = c('Community7','Community8'),
                                 other = 'Community5')
dominants[['Citrate']] <- list(isolate_2 = c('Community3','Community5'),
                               isolate_6 = 'Community1',
                               isolate_7 = c('Community2','Community7'),
                               isolate_8 = 'Community4',
                               other = c('Community6','Community8'))

### Now we have an additional decision to make. Isolates 2, 4 and 6 share a same
### ESV (even though isolates 4 and 6 carry additional copies of the 16S rRNA
### gene mapping to different ESVs). We can use the matrix in
###    > isolates[['seq']]
### to deconvolute ESV abundances into species abundances in those samples
### containing isolates 2, 4 and 6 (see isolates_abundance_from_ESVs.R).
### However, this procedure is not standard and nothing similar seems to have
### been described in the literature. The standard approach would be to just
### consider isolates 2, 4 and 6 as indistinguishable, which would give us:

# should we deconvolute ESV abundance into species abundance? if not...
if (!deconvolute_ESVs) {
  
  # unify isolates 2, 4 and 6 in deconvolution matrix
  isolates[['seq']] <- isolates[['seq']][,
                                         ! colnames(isolates[['seq']]) %in% c('isolate_4','isolate_6')]
  
  # merge isolates groups
  for (cs in carbon_sources) {
    isolates[['groups']][[cs]][['isolate_2']] <- c(isolates[['groups']][[cs]][['isolate_2']],
                                                   isolates[['groups']][[cs]][['isolate_4']],
                                                   isolates[['groups']][[cs]][['isolate_6']])
    isolates[['groups']][[cs]][['isolate_4']] <- NULL
    isolates[['groups']][[cs]][['isolate_6']] <- NULL
  }
  
  # merge dominants
  dominants[['Glutamine']] <- list(isolate_1 = 'Community1',
                                   isolate_2 = c('Community2','Community4','Community6'),
                                   isolate_3 = 'Community3',
                                   isolate_5 = c('Community7','Community8'),
                                   other = 'Community5')
  dominants[['Citrate']] <- list(isolate_2 = c('Community1','Community3','Community5'),
                                 isolate_7 = c('Community2','Community7'),
                                 isolate_8 = 'Community4',
                                 other = c('Community6','Community8'))
  
}


### ----------------------------------------------------------------------
### PAIRWISE COMPETITION VS. COALESCENCE VS. DOMINANT INVADING ALONE
### ----------------------------------------------------------------------

# initialize plotting table
plot_this <- data.frame(carbon_source = character(0),
                        community_1 = character(0),
                        community_2 = character(0),
                        isolate_1 = character(0),
                        isolate_2 = character(0),
                        plate_pairwise = character(0),
                        well_pairwise = character(0),
                        f_pairwise = numeric(0),
                        plate_coalescence = character(0),
                        well_coalescence = character(0),
                        q_bray_curtis = numeric(0),
                        q_jensen_shannon = numeric(0),
                        q_jaccard = numeric(0),
                        q_endemic = numeric(0),
                        q_bray_curtis_cohort = numeric(0),
                        q_jensen_shannon_cohort = numeric(0),
                        q_jaccard_cohort = numeric(0),
                        q_endemic_cohort = numeric(0),
                        f_singleinv = numeric(0))

# run through carbon sources and communities
for (cs in carbon_sources) {
  for (i in 1:(length(community_names)-1)) {
    for (j in (i+1):length(community_names)) {
      
      # communities
      comm_1 <- community_names[i] # by convention, community 1 will be the invasive and 2 the resident
      comm_2 <- community_names[j]
      
      # dominants of the communities
      dom_1 <- names(dominants[[cs]])[grep(comm_1,dominants[[cs]])]
      dom_2 <- names(dominants[[cs]])[grep(comm_2,dominants[[cs]])]
      
      # proceed only if the dominants of both communities are defined and different
      if (dom_1 != dom_2 & dom_1 != 'other' & dom_2 != 'other') {
        
        # find wells of pairwise competition
        champ_1 <- gsub('Community','Champion',comm_1)
        champ_2 <- gsub('Community','Champion',comm_2)
        wells_pairwise <- metadata[((metadata$Comm1 == champ_1 & metadata$Comm2 == champ_2) | (metadata$Comm1 == champ_2 & metadata$Comm2 == champ_1)) &
                                   metadata$Carbon == cs
                                   ,]
        
        # composition of pairwise competition wells
        cc_pairwise <- otus[,wells_pairwise$Sample,drop=F]
        
        # deconvolution matrix (isolates-to-ESVs map of dominants)
        Q = isolates[['seq']][,c(dom_1,dom_2)]
        
        # fraction of dominant 1 (invasive) in pairwise competition
        f_pairwise <- rep(NA,ncol(cc_pairwise))
        for (k in 1:ncol(cc_pairwise)) {
          doms_abundance_pairwise <- species_composition_from_sequencing(Q,cc_pairwise[,k,drop=FALSE])
          doms_abundance_pairwise <- doms_abundance_pairwise[c(dom_1,dom_2),]
          doms_abundance_pairwise <- doms_abundance_pairwise/sum(doms_abundance_pairwise)
          f_pairwise[k] <- doms_abundance_pairwise[dom_1]
        }
        
        # find all instances of community coalescence
        wells_coalescence <- metadata[((metadata$Comm1 == comm_1 & metadata$Comm2 == comm_2) | (metadata$Comm1 == comm_2 & metadata$Comm2 == comm_1)) &
                                      metadata$StabilizingCarb1 == cs &
                                      metadata$StabilizingCarb2 == cs &
                                      metadata$Carbon == cs
                                      ,]
        
        # composition of coalesced community
        cc_coalesced <- otus[,wells_coalescence$Sample,drop=F]
        rownames(cc_coalesced) <- otus$seq_id
        
        # composition of invasive and resident communities (average across replicates)
        cc_invasive <- as.data.frame(rowMeans(communities[[cs]][[comm_1]]))
        cc_resident <- as.data.frame(rowMeans(communities[[cs]][[comm_2]]))
        
        # get species composition from ESV composition for invasive, resident and coalesced communities
        Q <- isolates[['seq']]
        
        cc_invasive <- species_composition_from_sequencing(Q,cc_invasive)
        cc_resident <- species_composition_from_sequencing(Q,cc_resident)
        
        cc_coalesced <- lapply(1:ncol(cc_coalesced),FUN=function(k) species_composition_from_sequencing(Q,cc_coalesced[,k,drop=FALSE]))
        cc_coalesced <- merge(cc_coalesced[[1]],
                              cc_coalesced[[2]],
                              by='row.names',
                              all=T)
        rownames(cc_coalesced) <- cc_coalesced$Row.names
        cc_coalesced <- cc_coalesced[,2:3]
        
        # discard species with abundance 0 in all communities
        all_comms <- merge(cc_invasive,cc_resident,by='row.names',all=T)
        rownames(all_comms) <- all_comms$Row.names
        all_comms <- all_comms[,2:3]
        
        all_comms <- merge(all_comms,cc_coalesced,by='row.names',all=T)
        rownames(all_comms) <- all_comms$Row.names
        all_comms <- all_comms[,2:5]
        all_comms[is.na(all_comms)] <- 0
        
        n <- rowSums(all_comms)>0
        cc_invasive <- all_comms[n,1,drop=FALSE]
        cc_resident <- all_comms[n,2,drop=FALSE]
        cc_coalesced <- all_comms[n,3:4,drop=FALSE]
        
        # check that the row ordering is correct
        stopifnot(all(rownames(cc_invasive)==rownames(cc_resident)))
        stopifnot(all(rownames(cc_invasive)==rownames(cc_coalesced)))
        
        # compositions of community cohorts (dominants removed)
        cc_invasive_cohort <- cc_invasive[!rownames(cc_invasive) %in% c(dom_1,dom_2),,drop=FALSE]
        cc_resident_cohort <- cc_resident[!rownames(cc_resident) %in% c(dom_1,dom_2),,drop=FALSE]
        cc_coalesced_cohort <- cc_coalesced[!rownames(cc_coalesced) %in% c(dom_1,dom_2),,drop=FALSE]
        
        # re-normalize cohort compositions
        cc_invasive_cohort <- cc_invasive_cohort/matrix(rep(colSums(cc_invasive_cohort),nrow(cc_invasive_cohort)),
                                                        nrow=nrow(cc_invasive_cohort),
                                                        byrow=TRUE)
        cc_resident_cohort <- cc_resident_cohort/matrix(rep(colSums(cc_resident_cohort),nrow(cc_resident_cohort)),
                                                        nrow=nrow(cc_resident_cohort),
                                                        byrow=TRUE)
        cc_coalesced_cohort <- cc_coalesced_cohort/matrix(rep(colSums(cc_coalesced_cohort),nrow(cc_coalesced_cohort)),
                                                          nrow=nrow(cc_coalesced_cohort),
                                                          byrow=TRUE)
        
        # relative similarity invasive-to-coalesced (bray_curtis)
        sim_invasive <- cc_similarity(cc_invasive,cc_coalesced,metric='bray_curtis')
        sim_resident <- cc_similarity(cc_resident,cc_coalesced,metric='bray_curtis')
        q_bray_curtis <- sim_invasive/(sim_invasive + sim_resident)
        
        # relative similarity invasive-to-coalesced (jaccard)
        sim_invasive <- cc_similarity(cc_invasive,cc_coalesced,metric='jaccard')
        sim_resident <- cc_similarity(cc_resident,cc_coalesced,metric='jaccard')
        q_jaccard <- sim_invasive/(sim_invasive + sim_resident)
        
        # relative similarity invasive-to-coalesced (jensen_shannon)
        sim_invasive <- cc_similarity(cc_invasive,cc_coalesced,metric='jensen_shannon')
        sim_resident <- cc_similarity(cc_resident,cc_coalesced,metric='jensen_shannon')
        q_jensen_shannon <- sim_invasive/(sim_invasive + sim_resident)
        
        # relative similarity invasive-to-coalesced (endemic species)
        q_endemic <- rep(NA,ncol(cc_coalesced))
        for (k in 1:ncol(cc_coalesced)) {
          q_endemic[k] <- endemic(cc_invasive,cc_resident,cc_coalesced[,k])
        }
        
        # relative similarity invasive-to-coalesced (bray_curtis, cohorts)
        sim_invasive <- cc_similarity(cc_invasive_cohort,cc_coalesced_cohort,metric='bray_curtis')
        sim_resident <- cc_similarity(cc_resident_cohort,cc_coalesced_cohort,metric='bray_curtis')
        q_bray_curtis_cohort <- sim_invasive/(sim_invasive + sim_resident)
        
        # relative similarity invasive-to-coalesced (jaccard, cohorts)
        sim_invasive <- cc_similarity(cc_invasive_cohort,cc_coalesced_cohort,metric='jaccard')
        sim_resident <- cc_similarity(cc_resident_cohort,cc_coalesced_cohort,metric='jaccard')
        q_jaccard_cohort <- sim_invasive/(sim_invasive + sim_resident)
        
        # relative similarity invasive-to-coalesced (jensen_shannon, cohorts)
        sim_invasive <- cc_similarity(cc_invasive_cohort,cc_coalesced_cohort,metric='jensen_shannon')
        sim_resident <- cc_similarity(cc_resident_cohort,cc_coalesced_cohort,metric='jensen_shannon')
        q_jensen_shannon_cohort <- sim_invasive/(sim_invasive + sim_resident)
        
        # relative similarity invasive-to-coalesced (endemic species, cohorts)
        q_endemic_cohort <- rep(NA,ncol(cc_coalesced_cohort))
        for (k in 1:ncol(cc_coalesced_cohort)) {
          q_endemic_cohort[k] <- endemic(cc_invasive_cohort,cc_resident_cohort,cc_coalesced_cohort[,k])
        }
        
        # frequency of invasive dominant invading resident community alone
        wells_singleinv <- metadata[((metadata$Comm1 == champ_1 & metadata$Comm2 == comm_2) | (metadata$Comm1 == comm_2 & metadata$Comm2 == champ_1)) &
                                      metadata$Carbon == cs
                                      ,]
        cc_singleinv <- otus[,wells_singleinv$Sample,drop=F]
        
        # invasive dominant ESV composition
        Q = isolates[['seq']]
        
        # fraction of dominant 1 (invasive) invading resident community alone
        f_singleinv <- rep(NA,ncol(cc_singleinv))
        for (k in 1:ncol(cc_singleinv)) {
          doms_abundance_singleinv <- species_composition_from_sequencing(Q,cc_singleinv[,k,drop=FALSE])
          f_singleinv[k] <- doms_abundance_singleinv[dom_1,]
        }
        
        # fraction of dominant 1 (invasive) invading resident community with cohort
        f_multiinv <- as.numeric(cc_coalesced[dom_1,])
        
        # add to plotting table
        plot_this <- rbind(plot_this,
                           data.frame(
                             carbon_source = cs,
                             community_1 = comm_1,
                             community_2 = comm_2,
                             isolate_1 = dom_1,
                             isolate_2 = dom_2,
                             plate_pairwise = wells_pairwise$Experiment,
                             well_pairwise = wells_pairwise$Well,
                             f_pairwise = f_pairwise,
                             plate_coalescence = wells_coalescence$Experiment,
                             well_coalescence = wells_coalescence$Well,
                             q_bray_curtis = q_bray_curtis,
                             q_jensen_shannon = q_jensen_shannon,
                             q_jaccard = q_jaccard,
                             q_endemic = q_endemic,
                             q_bray_curtis_cohort = q_bray_curtis_cohort,
                             q_jensen_shannon_cohort = q_jensen_shannon_cohort,
                             q_jaccard_cohort = q_jaccard_cohort,
                             q_endemic_cohort = q_endemic_cohort,
                             f_singleinv = f_singleinv,
                             f_multiinv = f_multiinv
                           ))
        
      }
      
    }
  }
}

# check that sample positioning in wells is consistent
stopifnot(all(plot_this$well_pairwise == plot_this$well_coalescence))

# reshape plotting table (this makes it easier to create multipanel plots)
plot_this <- gather(plot_this,metric,value,q_bray_curtis:q_endemic_cohort)

# characters as factors
plot_this$carbon_source <- factor(plot_this$carbon_source,levels=c('Glutamine','Citrate'))
plot_this$community_1 <- factor(plot_this$community_1,levels=community_names)
plot_this$community_2 <- factor(plot_this$community_2,levels=community_names)
plot_this$metric <- factor(plot_this$metric,levels=unique(plot_this$metric))

# parameters for annotations of plots
poly <- list(w=0.10,
             aperture=0.15,
             color=c(rgb(220,245,220,maxColorValue=255),
                     rgb(230,230,230,maxColorValue=255),
                     rgb(250,220,220,maxColorValue=255)))

# make plots
myplots[['q-vs-pairwise_bray-curtis']] <-
  ggplot(data=plot_this[plot_this$metric=='q_bray_curtis',],
         aes(x=f_pairwise,y=value,
             color=carbon_source)) +
  geom_point(size=3,
             shape=1,
             stroke=0.5) +
  geom_smooth(formula = y ~ x,
              method = 'lm',
              se = FALSE,
              size = 0.5) +
  scale_y_continuous(name='Q\nCoalesced - Invasive',
                     limits=c(0,1),
                     breaks=c(0,0.5,1),
                     labels=c('0','0.5','1')) +
  scale_x_continuous(name='Frequency of\ninvasive dominant species\nin pairwise competition',
                     limits=c(0,1),
                     breaks=c(0,0.5,1),
                     labels=c('0','0.5','1')) +
  scale_color_manual(values=pl_carbon) +
  theme_bw() +
  theme(panel.grid=element_blank(),
        legend.title=element_blank(),
        legend.position=c(0.2,0.9),
        legend.background=element_rect(fill='transparent'),
        text=element_text(size=15),
        axis.text=element_text(size=15),
        axis.line=element_blank(),
        axis.ticks=element_line(size=0.25),
        panel.border=element_rect(size=0.25)) +
  coord_fixed() # +
  # ggtitle('Bray-Curtis similarity')

myplots[['q-vs-pairwise_other-metrics']] <-
  ggplot(data=plot_this[plot_this$metric!='q_bray_curtis' & !grepl('cohort',plot_this$metric),],
         aes(x=f_pairwise,y=value,
             color=carbon_source)) +
  geom_point(size=3,
             shape=1,
             stroke=0.5) +
  geom_smooth(formula = y ~ x,
              method = 'lm',
              se = FALSE,
              size = 0.5) +
  facet_grid(~ metric,
             labeller=labeller(metric=setNames(c('Jaccard similarity',
                                                 'Jensen-Shannon\nsimilarity\n(1 - distance)',
                                                 'Endemic species\nsurvival'),
                                               c('q_jaccard',
                                                 'q_jensen_shannon',
                                                 'q_endemic')))) +
  scale_y_continuous(name='Q\nCoalesced - Invasive',
                     limits=c(0,1),
                     breaks=c(0,0.5,1),
                     labels=c('0','0.5','1')) +
  scale_x_continuous(name='Frequency of invasive dominant species\nin pairwise competition',
                     limits=c(0,1),
                     breaks=c(0,0.5,1),
                     labels=c('0','0.5','1')) +
  scale_color_manual(values=pl_carbon) +
  theme_bw() +
  theme(panel.grid=element_blank(),
        legend.title=element_blank(),
        legend.background=element_rect(fill='transparent'),
        text=element_text(size=15),
        axis.line=element_blank(),
        axis.ticks=element_line(size=0.25),
        panel.border=element_rect(size=0.25),
        axis.text=element_text(size=15),
        strip.text=element_text(hjust=-0.01,
                                vjust=-0.01),
        strip.background=element_rect(fill='transparent',
                                      color='transparent')) +
  coord_fixed()

myplots[['q-vs-pairwise_cohorts']] <-
  ggplot(data=plot_this[grepl('cohort',plot_this$metric),],
         aes(x=f_pairwise,y=value,
             color=carbon_source)) +
  geom_point(size=3,
             shape=1,
             stroke=0.5) +
  geom_smooth(formula = y ~ x,
              method = 'lm',
              se = FALSE,
              size = 0.5) +
  facet_wrap(~ metric,
             nrow=2,
             labeller=labeller(metric=setNames(c('Bray-Curtis similarity',
                                                 'Jaccard similarity',
                                                 'Jensen-Shannon similarity\n(1 - distance)',
                                                 'Endemic species survival'),
                                               c('q_bray_curtis_cohort',
                                                 'q_jaccard_cohort',
                                                 'q_jensen_shannon_cohort',
                                                 'q_endemic_cohort')))) +
  scale_y_continuous(name='Q\nCoalesced - Invasive',
                     limits=c(0,1),
                     breaks=c(0,0.5,1),
                     labels=c('0','0.5','1')) +
  scale_x_continuous(name='Frequency of invasive dominant species\nin pairwise competition',
                     limits=c(0,1),
                     breaks=c(0,0.5,1),
                     labels=c('0','0.5','1')) +
  scale_color_manual(values=pl_carbon) +
  theme_bw() +
  theme(panel.grid=element_blank(),
        legend.title=element_blank(),
        legend.background=element_rect(fill='transparent'),
        text=element_text(size=15),
        axis.text=element_text(size=15),
        axis.line=element_blank(),
        axis.ticks=element_line(size=0.25),
        panel.border=element_rect(size=0.25),
        strip.text=element_text(hjust=-0.01,
                                vjust=-0.01),
        strip.background=element_rect(fill='transparent',
                                      color='transparent')) +
  coord_fixed()

myplots[['alone-vs-together']] <-
  ggplot(data=unique(plot_this[,c('carbon_source',
                                  'community_1',
                                  'community_2',
                                  'f_singleinv',
                                  'f_multiinv')]),
         aes(x=f_singleinv,y=f_multiinv,
             color=carbon_source)) +
  annotate('polygon',
           x=c(-0.1,poly$w,poly$w,-0.1,-0.1),
           y=c(-0.1,poly$w,1.1,1.1,-0.1),
           fill=poly$color[1]) +
  annotate('polygon',
           x=c(-0.1,poly$w,1.1,1.1,-0.1),
           y=c(-0.1,poly$w,poly$w,-0.1,-0.1),
           fill=poly$color[3]) +
  annotate('polygon',
           x=c(-0.1,1.1,1.1,1.1-1.2*poly$aperture,-0.1),
           y=c(-0.1,1.1-1.2*poly$aperture,1.1,1.1,-0.1),
           fill=poly$color[2]) +
  geom_abline(intercept=0,
              slope=1,
              color='black', 
              linetype='dashed',
              size=0.25) +
  geom_point(size=3,
             shape=1,
             stroke=0.5) +
  scale_y_continuous(name='Frequency of dominant\nspecies invading with cohort',
                     limits=c(-1,2),
                     breaks=c(0,0.5,1),
                     labels=c('0','0.5','1')) +
  scale_x_continuous(name='Frequency of dominant\nspecies invading alone\n ',
                     limits=c(-1,2),
                     breaks=c(0,0.5,1),
                     labels=c('0','0.5','1')) +
  scale_color_manual(values=pl_carbon) +
  theme_bw() +
  theme(panel.grid=element_blank(),
        legend.title=element_blank(),
        legend.position=c(0.2,0.9),
        legend.background=element_rect(fill='transparent'),
        text=element_text(size=15),
        axis.text=element_text(size=15),
        axis.line=element_blank(),
        axis.ticks=element_line(size=0.25),
        panel.border=element_rect(size=0.25)) +
  coord_fixed(xlim=c(0,1),
              ylim=c(0,1))

# display and save plots
if (display_plots) {
  print(myplots[['q-vs-pairwise_bray-curtis']])
  print(myplots[['q-vs-pairwise_other-metrics']])
  print(myplots[['q-vs-pairwise_cohorts']])
  print(myplots[['alone-vs-together']])
}
if (save_plots) {
  ggsave(file.path('.','plots','q-vs-pairwise_bray-curtis.pdf'),
         plot=myplots[['q-vs-pairwise_bray-curtis']],
         device='pdf',
         height=90,
         width=90,
         units='mm',dpi=600)
  ggsave(file.path('.','plots','q-vs-pairwise_other-metrics.pdf'),
         plot=myplots[['q-vs-pairwise_other-metrics']],
         device='pdf',
         height=90,
         width=180,
         units='mm',dpi=600)
  ggsave(file.path('.','plots','q-vs-pairwise_cohorts.pdf'),
         plot=myplots[['q-vs-pairwise_cohorts']],
         device='pdf',
         height=180,
         width=180,
         units='mm',dpi=600)
  ggsave(file.path('.','plots','alone-vs-together.pdf'),
         plot=myplots[['alone-vs-together']],
         device='pdf',
         height=90,
         width=90,
         units='mm',dpi=600)
}


### ----------------------------------------------------------------------
### SIMULATIONS DATA
### ----------------------------------------------------------------------

# load data from simulations
plot_this <- read.table('simul_data.txt',header=TRUE)

# make plots
myplots[['q-vs-pairwise_bray-curtis_simul']] <-
  ggplot(data=plot_this[!is.na(plot_this$f_pairwise),],
         aes(x=f_pairwise,y=q_bray_curtis)) +
  geom_point(size=3,
             shape=1,
             stroke=0.5,
             color='black') +
  geom_smooth(formula = y ~ x,
              method = 'lm',
              se = FALSE,
              size = 0.5,
              color='black') +
  scale_y_continuous(name='Q\nCoalesced - Invasive',
                     limits=c(0,1),
                     breaks=c(0,0.5,1),
                     labels=c('0','0.5','1')) +
  scale_x_continuous(name='Frequency of\ninvasive dominant species\nin pairwise competition',
                     limits=c(0,1),
                     breaks=c(0,0.5,1),
                     labels=c('0','0.5','1')) +
  theme_bw() +
  theme(panel.grid=element_blank(),
        text=element_text(size=15),
        axis.text=element_text(size=15),
        axis.line=element_blank(),
        axis.ticks=element_line(size=0.25),
        panel.border=element_rect(size=0.25)) +
  coord_fixed() # +
# ggtitle('Bray-Curtis similarity')

myplots[['alone-vs-together_simul']] <-
  ggplot(data=plot_this,
         aes(x=f_singleinv,y=f_coalescence)) +
  annotate('polygon',
           x=c(-0.1,poly$w,poly$w,-0.1,-0.1),
           y=c(-0.1,poly$w,1.1,1.1,-0.1),
           fill=poly$color[1]) +
  annotate('polygon',
           x=c(-0.1,poly$w,1.1,1.1,-0.1),
           y=c(-0.1,poly$w,poly$w,-0.1,-0.1),
           fill=poly$color[3]) +
  annotate('polygon',
           x=c(-0.1,1.1,1.1,1.1-1.2*poly$aperture,-0.1),
           y=c(-0.1,1.1-1.2*poly$aperture,1.1,1.1,-0.1),
           fill=poly$color[2]) +
  geom_abline(intercept=0,
              slope=1,
              color='black', 
              linetype='dashed',
              size=0.25) +
  geom_point(size=3,
             shape=1,
             stroke=0.5,
             color='black') +
  scale_y_continuous(name='Frequency of dominant\nspecies invading with cohort',
                     limits=c(-1,2),
                     breaks=c(0,0.5,1),
                     labels=c('0','0.5','1')) +
  scale_x_continuous(name='Frequency of dominant\nspecies invading alone\n ',
                     limits=c(-1,2),
                     breaks=c(0,0.5,1),
                     labels=c('0','0.5','1')) +
  theme_bw() +
  theme(panel.grid=element_blank(),
        text=element_text(size=15),
        axis.text=element_text(size=15),
        axis.line=element_blank(),
        axis.ticks=element_line(size=0.25),
        panel.border=element_rect(size=0.25)) +
  coord_fixed(xlim=c(0,1),
              ylim=c(0,1))

# display and save plots
if (display_plots) {
  print(myplots[['q-vs-pairwise_bray-curtis_simul']])
  print(myplots[['alone-vs-together_simul']])
}
if (save_plots) {
  ggsave(file.path('.','plots','q-vs-pairwise_bray-curtis_simul.pdf'),
         plot=myplots[['q-vs-pairwise_bray-curtis_simul']],
         device='pdf',
         height=90,
         width=90,
         units='mm',dpi=600)
  ggsave(file.path('.','plots','alone-vs-together_simul.pdf'),
         plot=myplots[['alone-vs-together_simul']],
         device='pdf',
         height=90,
         width=90,
         units='mm',dpi=600)
}













