###### Creation of lists with sample numbers for the following subsets:
## - Cross-sectional
## - Three timepoints
## - Two timepoints 
## - Two timepoints D0-D18
## - Two timepoints D0-D48
## - Two timepoints D18-D48.
## Create lists with sample_ID's, MOMID ID and timepoint 
################
load('C:/Users/arjen/Documents/research/UCLA/METSIM/data/R/FINAL DATA/cleaned datasets/metsim_microbiome_pso.rda')

x <- data.frame(sample_data(ps2)) %>% select(c('Subject', 'metsim_id', 'timepoint'))

str(x)

x$MOMID_timepoint <- paste0(x$metsim_id, '_', x$timepoint)

ps2
o <- otu_table(ps2)
rownames(o)

f<- as.data.frame.matrix(table(x$metsim_id, x$timepoint))
f$sum <- f$D0+f$D18+f$D48

f$METSIM_ID <- rownames(f)


# Select cross-sectional cohort
f$first <- ifelse(f$D0==1, 'D0',
                  ifelse(f$D18==1, "D18",
                  ifelse(f$D48==1, "D48", "missing")))
f$MOMID_timepointfirst <- paste0(f$METSIM_ID,'_', f$first)

sampleID <- subset(x, x$MOMID_timepoint %in% f$MOMID_timepointfirst)
write.csv(sampleID, file = 'sample_list_cross_sectional.csv')


# Select samples for triple cohort
triple <- subset(f, sum == 3)
triplesampleID <- subset(x, x$metsim_id %in% triple$METSIM_ID)

write.csv(triplesampleID, 'sample_list_tripletimepoints.csv')

# Select samples for double cohort
double <- subset(f, sum == 2)
doublesampleID <- subset(x, x$metsim_id %in% double$METSIM_ID)

write.csv(doublesampleID, 'sample_list_doubletimepoints.csv')

D018 <- subset(double, double$D0==1 & double$D18==1)
D018sampleID <- subset(x, x$metsim_id %in% D018$METSIM_ID)



d <- read.csv('C:/Users/arjen/Documents/research/UCLA/METSIM/data/R/FINAL DATA/cleaned datasets/FINAL_MICROBIOME_DATASET.csv')
`%notin%` <- negate(`%in%`)

t <- data.frame(sample_data(ps2))
exclude <- d$Sample_ID[(d$Sample_ID %notin% t$Subject)]
d <- subset(d, Sample_ID %notin% exclude)

write.csv(t, file = "ps2_sample_data.csv")

rownames(d) <- d$Sample_ID
d2 <- sample_data(d)

otu <- otu_table(ps2)
taxtable <- tax_table(ps2)
refsequence <- refseq(ps2)
ps3 <-merge_phyloseq(otu, taxtable, refsequence, d2)

save(ps3, file = 'C:/Users/arjen/Documents/research/UCLA/METSIM/data/R/FINAL DATA/cleaned datasets/metsim_microbiome_pso_updated.rda')
load(file = 'C:/Users/arjen/Documents/research/UCLA/METSIM/data/R/FINAL DATA/cleaned datasets/metsim_microbiome_pso_updated.rda')
l <- read.csv('C:/Users/arjen/Documents/research/UCLA/METSIM/data/R/FINAL DATA/cleaned datasets/sample_list_cross_sectional.csv')

crosssec <- prune_samples(ps3, samples = l$Subject)
