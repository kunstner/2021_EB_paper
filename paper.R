
# Set up output -----------------------------------------------------------

outdir <- 'plots_paper'
dir.create(file.path("./", outdir), showWarnings = FALSE)

# libraries needed --------------------------------------------------------

library(tidyverse)
library(phyloseq)
library(ggplot2)
library(ggpubr)
library(microbiome)
library(pals)
library(gridExtra)

# additional libraries
# library(sjPlot)
# library(performance)
# library(car)

# get data ----------------------------------------------------------------

ps <- readRDS(file = "data/phyloseq.ASV.RDS")
ps

# Estimate Alpha diversity ------------------------------------------------

sample_data(ps)$Shannon <- estimate_richness(ps, measures = "Shannon")$Shannon

# Change RDEB subtype for selected cases ----------------------------------

# EB20, EB35, EB36, EB37 (intermediate -> severe)
ids <- c('EB20', 'EB35', 'EB36', 'EB37')
cova <- sample_data(ps) %>% data.frame()
cova$Severity[cova$MINIONID %in% ids] <- 'severe'
sample_data(ps) <- cova
rm(cova); rm(ids)

# Set age groups ----------------------------------------------------------

sample_data(ps)$AgeGroup <- NA

cova <- sample_data(ps) %>% data.frame()
cova$AgeGroup[ cova$Age < 5 ] <- "01"
cova$AgeGroup[ cova$Age >= 5 & cova$Age < 11 ] <- "02"
cova$AgeGroup[ cova$Age >= 11 ] <- "03"

table(cova$AgeGroup)

cova$Age_Case <- paste(cova$AgeGroup, cova$Location, cova$Case, sep = ":")
cova$Location_Severity <- paste(cova$Location, cova$Severity, sep = ":")

sample_data(ps) <- cova

# Add new parameters ------------------------------------------------------

cova <- gdata::read.xls(xls = "data/MINION-EB_DEB-data.xlsx")
ps_sd <- sample_data(ps) %>% data.frame()
ps_sd$ID <- rownames(ps_sd)
dim(cova)

cova <- cova[, c(1, 38:67)] |>
    dplyr::select(-Stool_taken, -Stool_retrieved, -Stuhlmodus) %>%
    dplyr::mutate( Leukozytes = as.numeric(Leukozytes) )

# Transport stool
cova$Stool_delay %>% summary
cova$Stool_delay %>% sd(., na.rm = TRUE)

cova <- merge(ps_sd, cova, by.x = 'MINIONID', by.y = 'MINIONID', all.x = T, all.y = F)
rownames(cova) <- cova$ID
cova$ID <- NULL
head(cova)

sample_data(ps) <- cova

rm(cova); rm(ps_sd)

# Add species information S. aureus ---------------------------------------

staph <- gdata::read.xls(xls = "data/Staph_blast_asv.xlsx")
staph <- staph %>%
    dplyr::filter(Blast > 99)
tt <- tax_table(ps) %>% data.frame
tt <- merge(x = tt[, -7], y = staph[ c("X", "Species")], by.x = 0, by.y = "X", all.x = TRUE)
rownames(tt) <- tt$Row.names
tt$Row.names <- NULL
tt <- tt[ rownames(data.frame(tax_table(ps))), ]
sum(rownames(tt) == rownames(data.frame(tax_table(ps)))) # check
tax_table(ps) <- as.matrix(tt)

rm(tt)
rm(staph)

# Define global variables and colors --------------------------------------

seedID <- 1253
set.seed(seedID)

cols_alpha <- c("grey25", "grey75")
cols_heatmap <- circlize::colorRamp2(c(-1, 0, 1),  c("#2C7BB6", "white", "#D7191C"))
cols_serverity <- c("steelblue1", "red") #RDEBintermed, RDEBsevere

taxNumberPlot <- 20

cols <- c( pals::glasbey(n = 32), "#666666", "#CCCCCC")
cols[1] <- "#104E8B"
cols[2] <- "#8B1A1A"
cols[3] <- "#008B00"

# pie( rep( 1,length(cols)), col=(cols) )

top_genera <- c("Abiotrophia", "Actinomyces", "Aggregatibacter", "Alistipes", "Anaerococcus", "Bacteroides", "Bifidobacterium", "Blautia A", "Corynebacterium", "Dialister", "Faecalibacterium", "Gemella", "Gemmiger", "Granulicatella", "Haemophilus", "Haemophilus D", "Kocuria", "Lentilactobacillus", "Leptotrichia", "Micrococcus", "Moraxella A", "Neisseria", "Oribacterium", "Phocaeicola", "Prevotella", "Pseudomonas H", "Rothia", "Ruminococcus E", "Staphylococcus", "Streptococcus", "Sutterella", "Unknown", "Veillonella", "Other")

top_phyla <- c("Acidobacteriota", "Actinobacteriota", "Bacteroidota", "Bdellovibrionota", "Campylobacterota", "Deinococcota", "Desulfobacterota", "Firmicutes", "Firmicutes A", "Firmicutes B", "Firmicutes C", "Fusobacteriota", "Planctomycetota", "Proteobacteria", "Spirochaetota", "Verrucomicrobiota", "Other")

cols_top_genera <- cols[1:length(top_genera)]
names(cols_top_genera) <- sort(top_genera)

cols_top_phyla <- cols[1:length(top_phyla)]
names(cols_top_phyla) <- sort(top_phyla)

# Functions ---------------------------------------------------------------

getTaxaAbund <- function(taxRank, phySeq, comp) {

    # get number of contigs for each sample
    contigs <- sample_sums(phySeq)

    # get desired taxonomic rank data
    tg.otu <- microbiome::aggregate_taxa(x = phySeq, level = taxRank) %>% otu_table() %>%  data.frame()

    # any contigs lost due to missing annotations?
    # contigs - apply(tg.otu, 2, sum)
    if ( sum(contigs - apply(tg.otu, 2, sum)) > 0 ) {
        tg.otu <- rbind(tg.otu, (contigs - apply(tg.otu, 2, sum)) )
        rownames(tg.otu)[dim(tg.otu)[1]] <- "NA"
    }

    # transpose data frame
    tg.otu <- data.frame(t(tg.otu))
    # get the numer of taxons
    tax.number <- dim(tg.otu)[2]
    # just sum the number of contigs for each sample
    tg.otu$Contigs <- apply(X = tg.otu, MARGIN = 1, sum)

    rownames(tg.otu) <- rownames( data.frame(sample_data(phySeq)) )


    tg.otu$id <- rownames(tg.otu)

    # get meta data...
    tg.meta <- data.frame(sample_data(phySeq))
    tg.meta$id <- rownames(tg.meta)
    #...and merge
    tg.otu <- merge(tg.otu, tg.meta, by="id")

    compVector <- as.factor(tg.otu[ ,comp])

    ret.df <- data.frame(Tax = as.character())

    # iterate through all taxa
    for ( i in 2:(tax.number+1) ) {
        #  print( colnames(tg.otu)[i] )

        dummy.df <- NULL
        dummy.df <- data.frame(Tax=NULL, Group=NULL,
                               mean=NULL, sd=NULL, median=NULL, mad=NULL)

        for( j in 1 : length(levels(compVector)) ) {
            #  print(levels(compVector)[j])
            t.level  <- levels(compVector)[j]
            t.mean   <- mean( tg.otu[ tg.otu[ , eval(comp)  ] == t.level, i],
                              na.rm = T )
            t.sd    <- sd( tg.otu[ tg.otu[ , eval(comp)  ] == t.level, i],
                           na.rm = T )
            t.median <- median( tg.otu[ tg.otu[ , eval(comp)  ] == t.level, i],
                                na.rm = T )
            t.mad    <- mad( tg.otu[ tg.otu[ , eval(comp)  ] == t.level, i],
                             na.rm = T )

            dummy <- data.frame( colnames(tg.otu)[i], levels(compVector)[j],
                                 t.mean , t.sd, t.median, t.mad)

            colnames(dummy) <- c("Tax", "Group", "mean", "sd", "median", "mad")

            dummy.df <- rbind(dummy.df, dummy )
        }

        dummy.df <- as.data.frame(dummy.df)
        ret.df <- rbind(ret.df, dummy.df )

    }
    ret.df$Tax <- factor(ret.df$Tax)
    return(ret.df)
}

# Figure 1 ----------------------------------------------------------------

ps_ra <- transform_sample_counts(ps, function(x){x / sum(x)})
ps_ra <- subset_samples(physeq = ps_ra, Location == "Forearm" | Location == "Wound")

#
# POOLED SAMPLES
#

taxRank <- "Phylum"
# get abundance data
data <- getTaxaAbund(taxRank, ps_ra, "Location_Case")
# sort taxonomic levels alphabetically
data$Tax <- factor(data$Tax, levels = sort(levels(data$Tax)))
# identify most abundant taxa across all samples
data.aggr <- aggregate( formula = mean ~ Tax, data = data, FUN = mean)
# select n most abundant taxa
best.tax <- data.aggr[ order(-data.aggr$mean), ][1:taxNumberPlot, 1]
# create new data frame
data.plot <- data[ data$Tax %in% best.tax, ]
data.plot <- data.plot[ data.plot$Tax != "NA.",]

levels(data.plot$Tax)

other.df <- data.frame(
    Tax = "Other",
    Group = levels(factor(data.plot$Group)),
    mean = aggregate( formula = mean ~ Group, data = data.plot, FUN = function(x) 1 - sum(x) )[,2],
    sd = NA, median = NA, mad = NA
)
if ( sum(other.df$mean) > 0 ) {
    data.plot <- rbind( data.plot, other.df)
}

data.plot$Tax <- gsub(pattern = "_", replacement = " ", x = data.plot$Tax)
data.plot$Group <- gsub(pattern = "Forearm", replacement = "Unwounded", x = data.plot$Group)
data.plot$Group <- gsub(pattern = ":FALSE", replacement = " CTRL", x = data.plot$Group)
data.plot$Group <- gsub(pattern = ":TRUE",  replacement = " RDEB",   x = data.plot$Group)

fig1_1 <- ggplot(data.plot, aes(x=Group, y=mean, fill=Tax) ) +
    facet_wrap(~., scales="free_x") +
    geom_bar(width=0.75, stat='identity' ) +
    theme_bw() +
    ylim(c(0,1)) +
    theme(axis.line = element_line(colour = "black"),
          legend.text = element_text(face = "italic"),
          legend.justification = "top",
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank(),
          axis.text.x = element_text(angle = 45, hjust = 1, size=12),
          axis.text.y = element_text(size=12),
          strip.text.x = element_text(angle = 0, size = 10, color = "white"),
          axis.ticks.x=element_blank(),
          strip.background = element_rect(colour="white", fill="white", size=1.5, linetype="solid")) +
    ggtitle( paste("", "", sep=" ") ) +
    xlab("") +
    ylab("Proportion mapped contigs") +
    guides(fill = guide_legend(title = taxRank)) +
    scale_fill_manual(values=cols_top_phyla[names(cols_top_phyla) %in% data.plot$Tax])
fig1_1
ggsave(filename = paste0( outdir, "/Fig1_all.", taxRank,".pdf"), width = 4.5, height = 9)

taxRank <- "Genus"
# get abundance data
data <- getTaxaAbund(taxRank, ps_ra, "Location_Case")
# sort taxonomic levels alphabetically
data$Tax <- factor(data$Tax, levels = sort(levels(data$Tax)))
# identify most abundant taxa across all samples
data.aggr <- aggregate( formula = mean ~ Tax, data = data, FUN = mean)
# select n most abundant taxa
best.tax <- data.aggr[ order(-data.aggr$mean), ][1:taxNumberPlot, 1]
# create new data frame
data.plot <- data[ data$Tax %in% best.tax, ]
data.plot <- data.plot[ data.plot$Tax != "NA.",]

levels(data.plot$Tax)

other.df <- data.frame(
    Tax = "Other",
    Group = levels(factor(data.plot$Group)),
    mean = aggregate( formula = mean ~ Group, data = data.plot, FUN = function(x) 1 - sum(x) )[,2],
    sd = NA, median = NA, mad = NA
)
if ( sum(other.df$mean) > 0 ) {
    data.plot <- rbind( data.plot, other.df)
}

data.plot$Tax <- gsub(pattern = "_", replacement = " ", x = data.plot$Tax)
data.plot$Group <- gsub(pattern = "Forearm", replacement = "Unwounded", x = data.plot$Group)
data.plot$Group <- gsub(pattern = ":FALSE", replacement = " CTRL", x = data.plot$Group)
data.plot$Group <- gsub(pattern = ":TRUE",  replacement = " RDEB",   x = data.plot$Group)

fig1_2 <- ggplot(data.plot, aes(factor(Tax), x=Group, y=mean, fill=Tax) ) +
    facet_wrap(~., scales="free_x") +
    geom_bar(width=0.75, stat='identity' ) +
    theme_bw() +
    ylim(c(0,1)) +
    theme(axis.line = element_line(colour = "black"),
          legend.text = element_text(face = "italic"),
          legend.justification = "top",
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank(),
          axis.text.x = element_text(angle = 45, hjust = 1, size=12),
          axis.text.y = element_text(size=12),
          strip.text.x = element_text(angle = 0, size = 10, color = "white"),
          axis.ticks.x=element_blank(),
          strip.background = element_rect(colour="white", fill="white", size=1.5, linetype="solid")) +
    ggtitle( paste("", "", sep=" ") ) +
    xlab("") +
    ylab("Proportion mapped contigs") +
    guides(fill = guide_legend(title = taxRank)) +
    scale_fill_manual(values=cols_top_genera[names(cols_top_genera) %in% data.plot$Tax])
fig1_2
ggsave(filename = paste0( outdir, "/Fig1_all.", taxRank,".pdf"), width = 4.5, height = 9)

#
# Age groups
#

taxRank <- "Phylum"
# get abundance data
data <- getTaxaAbund(taxRank, ps_ra, "Age_Case")
# sort taxonomic levels alphabetically
data$Tax <- factor(data$Tax, levels = sort(levels(data$Tax)))
# identify most abundant taxa across all samples
data.aggr <- aggregate( formula = mean ~ Tax, data = data, FUN = mean)
# select n most abundant taxa
best.tax <- data.aggr[ order(-data.aggr$mean), ][1:taxNumberPlot, 1]
# create new data frame
data.plot <- data[ data$Tax %in% best.tax, ]
data.plot <- data.plot[ data.plot$Tax != "NA.",]

levels(data.plot$Tax)

other.df <- data.frame(
    Tax = "Other",
    Group = levels(factor(data.plot$Group)),
    mean = aggregate( formula = mean ~ Group, data = data.plot, FUN = function(x) 1 - sum(x) )[,2],
    sd = NA, median = NA, mad = NA
)
if ( sum(other.df$mean) > 0 ) {
    data.plot <- rbind( data.plot, other.df)
}

data.plot$Tax <- gsub(pattern = "_", replacement = " ", x = data.plot$Tax)
data.plot$Group <- gsub(pattern = "Forearm", replacement = "Unwounded", x = data.plot$Group)
data.plot$Group <- gsub(pattern = ":FALSE", replacement = " CTRL", x = data.plot$Group)
data.plot$Group <- gsub(pattern = ":TRUE",  replacement = " RDEB",   x = data.plot$Group)

data.plot$AgeGroup <- NA
data.plot$AgeGroup[ grep(pattern = "01:", data.plot$Group) ] <-  "0-4 years"
data.plot$AgeGroup[ grep(pattern = "02:", data.plot$Group) ] <-  "5-10 years"
data.plot$AgeGroup[ grep(pattern = "03:", data.plot$Group) ] <-  "11-22 years"

data.plot$AgeGroup <- factor(data.plot$AgeGroup, levels = c("0-4 years", "5-10 years", "11-22 years"))

data.plot$Group <- gsub(pattern = "01:", replacement = "", x = data.plot$Group)
data.plot$Group <- gsub(pattern = "02:", replacement = "", x = data.plot$Group)
data.plot$Group <- gsub(pattern = "03:", replacement = "", x = data.plot$Group)

fig1_3 <- ggplot(data.plot, aes(x=Group, y=mean, fill=Tax) ) +
    facet_wrap(AgeGroup~., scales="free_x") +
    geom_bar(width=0.75, stat='identity' ) +
    theme_bw() +
    ylim(c(0,1)) +
    theme(axis.line = element_line(colour = "black"),
          legend.text = element_text(face = "italic"),
          legend.justification = "top",
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank(),
          axis.text.x = element_text(angle = 45, hjust = 1, size=12),
          strip.text.x = element_text(angle = 0, size = 10),
          axis.text.y = element_text(size=12),
          axis.ticks.x=element_blank(),
          strip.background = element_rect(colour="white", fill="white", size=1.5, linetype="solid")) +
    ggtitle( paste("", "", sep=" ") ) +
    xlab("") +
    ylab("Proportion mapped contigs") +
    guides(fill = guide_legend(title = taxRank)) +
    scale_fill_manual(values=cols_top_phyla[names(cols_top_phyla) %in% data.plot$Tax])
fig1_3
ggsave(filename = paste0( outdir, "/Fig1_age.", taxRank,".pdf"), width = 9, height = 9)

taxRank <- "Genus"
# get abundance data
data <- getTaxaAbund(taxRank, ps_ra, "Age_Case")
# sort taxonomic levels alphabetically
data$Tax <- factor(data$Tax, levels = sort(levels(data$Tax)))
# identify most abundant taxa across all samples
data.aggr <- aggregate( formula = mean ~ Tax, data = data, FUN = mean)
# select n most abundant taxa
best.tax <- data.aggr[ order(-data.aggr$mean), ][1:taxNumberPlot, 1]
# create new data frame
data.plot <- data[ data$Tax %in% best.tax, ]
data.plot <- data.plot[ data.plot$Tax != "NA.",]

levels(data.plot$Tax)

other.df <- data.frame(
    Tax = "Other",
    Group = levels(factor(data.plot$Group)),
    mean = aggregate( formula = mean ~ Group, data = data.plot, FUN = function(x) 1 - sum(x) )[,2],
    sd = NA, median = NA, mad = NA
)
if ( sum(other.df$mean) > 0 ) {
    data.plot <- rbind( data.plot, other.df)
}

data.plot$Tax <- gsub(pattern = "_", replacement = " ", x = data.plot$Tax)
data.plot$Group <- gsub(pattern = "Forearm", replacement = "Unwounded", x = data.plot$Group)
data.plot$Group <- gsub(pattern = ":FALSE", replacement = " CTRL", x = data.plot$Group)
data.plot$Group <- gsub(pattern = ":TRUE",  replacement = " RDEB",   x = data.plot$Group)

data.plot$AgeGroup <- NA
data.plot$AgeGroup[ grep(pattern = "01:", data.plot$Group) ] <-  "0-4 years"
data.plot$AgeGroup[ grep(pattern = "02:", data.plot$Group) ] <-  "5-10 years"
data.plot$AgeGroup[ grep(pattern = "03:", data.plot$Group) ] <-  "11-22 years"

data.plot$AgeGroup <- factor(data.plot$AgeGroup, levels = c("0-4 years", "5-10 years", "11-22 years"))

data.plot$Group <- gsub(pattern = "01:", replacement = "", x = data.plot$Group)
data.plot$Group <- gsub(pattern = "02:", replacement = "", x = data.plot$Group)
data.plot$Group <- gsub(pattern = "03:", replacement = "", x = data.plot$Group)

fig1_4 <- ggplot(data.plot, aes(factor(Tax), x=Group, y=mean, fill=Tax) ) +
    facet_wrap(AgeGroup~., scales="free_x") +
    geom_bar(width=0.75, stat='identity' ) +
    theme_bw() +
    ylim(c(0,1)) +
    theme(axis.line = element_line(colour = "black"),
          legend.text = element_text(face = "italic"),
          legend.justification = "top",
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank(),
          axis.text.x = element_text(angle = 45, hjust = 1, size=12),
          strip.text.x = element_text(angle = 0, size = 10),
          axis.text.y = element_text(size=12),
          axis.ticks.x=element_blank(),
          strip.background = element_rect(colour="white", fill="white", size=1.5, linetype="solid")) +
    ggtitle( paste("", "", sep=" ") ) +
    xlab("") +
    ylab("Proportion mapped contigs") +
    guides(fill = guide_legend(title = taxRank)) +
    scale_fill_manual(values=cols_top_genera[names(cols_top_genera) %in% data.plot$Tax])
fig1_4
ggsave(filename = paste0( outdir, "/Fig1_age.", taxRank,".pdf"), width = 9, height = 9)

#
# Alpha
#

alpha <- ps %>%
    sample_data() %>%
    data.frame() %>%
    dplyr::select( MINIONID, Location, Case, Location_Case, Age_Case, AgeGroup, Shannon ) %>%
    dplyr::filter( Location == "Forearm" | Location == "Wound") %>%
    dplyr::mutate(Location_Case = as.character(Location_Case))

alpha$Location_Case <- gsub(pattern = "Forearm", replacement = "Unwounded", x = alpha$Location_Case)
alpha$Location_Case <- gsub(pattern = ":FALSE", replacement = " CTRL", x = alpha$Location_Case)
alpha$Location_Case <- gsub(pattern = ":TRUE",  replacement = " RDEB",   x = alpha$Location_Case)
alpha$Location_Case <- factor(alpha$Location_Case, levels = c("Unwounded CTRL", "Unwounded RDEB", "Wound RDEB") )

alpha$AgeGroup <- as.character(alpha$AgeGroup)
alpha$AgeGroup[ grep(pattern = "01", alpha$AgeGroup) ] <-  "0-4 years"
alpha$AgeGroup[ grep(pattern = "02", alpha$AgeGroup) ] <-  "5-10 years"
alpha$AgeGroup[ grep(pattern = "03", alpha$AgeGroup) ] <-  "11-22 years"
alpha$AgeGroup <- factor(alpha$AgeGroup, levels = c("0-4 years", "5-10 years", "11-22 years"))

alpha$Age_Case <- gsub(pattern = "Forearm", replacement = "Unwounded", x = alpha$Age_Case)
alpha$Age_Case <- gsub(pattern = ":FALSE", replacement = " CTRL", x = alpha$Age_Case)
alpha$Age_Case <- gsub(pattern = ":TRUE",  replacement = " RDEB",   x = alpha$Age_Case)
alpha$Age_Case <- gsub(pattern = "01:",  replacement = "",   x = alpha$Age_Case)
alpha$Age_Case <- gsub(pattern = "02:",  replacement = "",   x = alpha$Age_Case)
alpha$Age_Case <- gsub(pattern = "03:",  replacement = "",   x = alpha$Age_Case)
alpha$Age_Case <- factor(alpha$Location_Case, levels = c("Unwounded CTRL", "Unwounded RDEB", "Wound RDEB") )

#
# POOLED
#

mycomp <- list(
    c("Unwounded CTRL", "Unwounded RDEB"),
    c("Unwounded RDEB", "Wound RDEB")
)

fig1_5 <- ggpubr::ggdotplot(data = alpha , x = "Location_Case", y = "Shannon",
                  fill = "Case", palette = cols_alpha,
                  size = 0.75,
                  xlab = "", ylab = "Shannon index of species diversity") +
    facet_wrap(~., scales="free_x") +
    theme(legend.position = "none",
          axis.text.x=element_text(size=12, angle=45, hjust = 1),
          axis.text.y=element_text(size=12, angle=0),
          strip.text.x = element_text(angle = 0, size = 10, color = "white"),
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          axis.line = element_line(colour = "black"),
          strip.background = element_rect(colour="white", fill="white", size=1.5, linetype="solid")) +
    ylim(0,7) +
    #ggpubr::stat_compare_means(method = "kruskal.test", label.y = 8) +
    stat_compare_means(comparisons = mycomp, label.y = c(6, 6))
fig1_5
ggsave(filename = paste0( outdir, "/Fig1_all.Alpha.pdf"), width = 4.5, height = 9)

#
# Age Groups
#

mycomp <- list(
    c("Unwounded CTRL", "Unwounded RDEB"),
    c("Unwounded RDEB", "Wound RDEB")
)

fig1_6 <- ggpubr::ggdotplot(data = alpha , x = "Age_Case", y = "Shannon",
                  fill = "Case", palette = cols_alpha,
                  size = 0.75,
                  xlab = "", ylab = "Shannon index of species diversity") +
    facet_wrap(AgeGroup~., scales="free_x") +
    theme(axis.line = element_line(colour = "black"),
          legend.position = "none",
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank(),
          axis.text.x = element_text(angle = 45, hjust = 1, size=12),
          strip.text.x = element_text(angle = 0, size = 12),
          axis.text.y = element_text(size=12),
          axis.ticks.x=element_blank(),
          strip.background = element_rect(colour="white", fill="white", size=1.5, linetype="solid")) +
    ylim(0,7) +
    #ggpubr::stat_compare_means(method = "kruskal.test", label.y = 8) +
    stat_compare_means(comparisons = mycomp, label.y = c(6, 6))
fig1_6
ggsave(filename = paste0( outdir, "/Fig1_age.Alpha.pdf"), width = 9, height = 9)

wilcox.test(
    alpha$Shannon[alpha$Location_Case == "Wound RDEB" & alpha$AgeGroup == "0-4 years"],
    alpha$Shannon[alpha$Location_Case == "Wound RDEB" & alpha$AgeGroup == "5-10 years"]
    )$p.value
wilcox.test(
    alpha$Shannon[alpha$Location_Case == "Wound RDEB" & alpha$AgeGroup == "5-10 years"],
    alpha$Shannon[alpha$Location_Case == "Wound RDEB" & alpha$AgeGroup == "11-22 years"]
)$p.value
wilcox.test(
    alpha$Shannon[alpha$Location_Case == "Wound RDEB" & alpha$AgeGroup == "0-4 years"],
    alpha$Shannon[alpha$Location_Case == "Wound RDEB" & alpha$AgeGroup == "11-22 years"]
)$p.value

#
# Correlations
#

# Phyla
corr_df <- phylosmith::variable_correlation( ps_ra, treatment = c('Location', 'Case'),
                                             subset = c('Forearm', 'Wound'),
                                             variables = 'Age',
                                             classification = "Phylum", method = 'spearman')
# corr_df <- corr_df %>%
#     dplyr::select(Group = Treatment, Phylum = X, rho) %>%
#     tidyr::spread(data = ., key = Phylum, value = rho) %>%
#     tibble::column_to_rownames(var = "Group") %>%
#     data.frame()
# corr_df[is.na(corr_df)] <- 0

corr_df$Treatment <- gsub(pattern = "Forearm", replacement = " Unwounded", x = corr_df$Treatment)
corr_df$Treatment <- gsub(pattern = " TRUE", replacement = " RDEB", x = corr_df$Treatment)
corr_df$Treatment <- gsub(pattern = " FALSE", replacement = " CTRL", x = corr_df$Treatment)

corr_df$p_clean <- corr_df$p
corr_df$p_clean <- format( round(corr_df$p_clean, 4), scientific = F ) # avoid exponentials
corr_df$p_clean[corr_df$p_clean >= 0.05] <- ""

fig1_7 <- corr_df %>%
    dplyr::mutate(X = gsub(pattern = "_", replacement = " ", X)) %>%
    ggplot(data = ., aes(x=Treatment, y=factor(X, levels = rev(levels(factor(X)))), fill=rho)) +
        geom_tile() +
        geom_text(aes(Treatment, X,
                      label = p_clean),
                  color = "black", size = 4) +
        scale_fill_gradient2(low = "#2C7BB6", high = "#D7191C", mid = "white",
                             midpoint = 0, limit = c(-1,1), space = "Lab",
                             name="Spearman\nCorrelation") +
        theme(axis.line = element_line(colour = "black"),
              legend.position = "right",
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              panel.border = element_blank(),
              panel.background = element_blank(),
              axis.text.x = element_text(angle = 45, hjust = 1, size=12),
              strip.text.x = element_text(angle = 0, size = 12),
              axis.text.y = element_text(size=12, face = "italic"),
              axis.ticks.x=element_blank(),
              strip.background = element_rect(colour="white", fill="white", size=1.5, linetype="solid")) +
        xlab("") + ylab("")
fig1_7
ggsave(filename = paste0( outdir, "/Fig1_all.Correlation_Phylum.pdf"), width = 4.5, height = 9)

# Genera
( ps.sel <- microbiome::aggregate_taxa(x = ps_ra, level = "Genus") %>% phylosmith::taxa_filter(., treatment = c('Location'), frequency = 0.2) )
corr_df <- phylosmith::variable_correlation( ps.sel, treatment = c('Location', 'Case'),
                                             subset = c('Forearm', 'Wound'),
                                             variables = 'Age',
                                             classification = "Genus", method = 'spearman')

corr_df$Treatment <- gsub(pattern = "Forearm", replacement = " Unwounded", x = corr_df$Treatment)
corr_df$Treatment <- gsub(pattern = " TRUE", replacement = " RDEB", x = corr_df$Treatment)
corr_df$Treatment <- gsub(pattern = " FALSE", replacement = " CTRL", x = corr_df$Treatment)

corr_df$p_clean <- corr_df$p
corr_df$p_clean <- format( round(corr_df$p_clean, 4), scientific = F ) # avoid exponentials
corr_df$p_clean[corr_df$p_clean >= 0.05] <- ""
corr_df$p_clean[corr_df$X == "Staphylococcus"] <- format( round(corr_df$p[corr_df$X == "Staphylococcus"], 4), scientific = F )

# keep only genera with p < 0.05
corr_genera <- corr_df$X[corr_df$p < 0.05]

fig1_8 <- corr_df %>%
    dplyr::filter(X %in% corr_genera | X == "Staphylococcus") %>%
    dplyr::mutate(X = gsub(pattern = "_", replacement = " ", X)) %>%
    ggplot(data = ., aes(x=Treatment, y=factor(X, levels = rev(levels(factor(X)))), fill=rho)) +
        geom_tile() +
        geom_text(aes(Treatment, X,
                      label = p_clean),
                  color = "black", size = 4) +
        scale_fill_gradient2(low = "#2C7BB6", high = "#D7191C", mid = "white",
                             midpoint = 0, limit = c(-1,1), space = "Lab",
                             name="Spearman\nCorrelation") +
        theme(axis.line = element_line(colour = "black"),
              legend.position = "right",
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              panel.border = element_blank(),
              panel.background = element_blank(),
              axis.text.x = element_text(angle = 45, hjust = 1, size=12),
              strip.text.x = element_text(angle = 0, size = 12),
              axis.text.y = element_text(size=12, face = "italic"),
              axis.ticks.x=element_blank(),
              strip.background = element_rect(colour="white", fill="white", size=1.5, linetype="solid")) +
        xlab("") + ylab("")
fig1_8
ggsave(filename = paste0( outdir, "/Fig1_all.Correlation_Genus.pdf"), width = 4.5, height = 9)

# Figure 2 ----------------------------------------------------------------

ps_ra <- transform_sample_counts(ps, function(x){x / sum(x)})
ps_ra <- subset_samples(physeq = ps_ra, Location == "Forearm" | Location == "Wound")
ps_ra <- subset_samples(physeq = ps_ra, Severity == "severe" | Severity == "interm")

#
# POOLED SAMPLES
#

taxRank <- "Phylum"
# get abundance data
data <- getTaxaAbund(taxRank, ps_ra, "Location_Severity")
# sort taxonomic levels alphabetically
data$Tax <- factor(data$Tax, levels = sort(levels(data$Tax)))
# identify most abundant taxa across all samples
data.aggr <- aggregate( formula = mean ~ Tax, data = data, FUN = mean)
# select n most abundant taxa
best.tax <- data.aggr[ order(-data.aggr$mean), ][1:taxNumberPlot, 1]
# create new data frame
data.plot <- data[ data$Tax %in% best.tax, ]
data.plot <- data.plot[ data.plot$Tax != "NA.",]

levels(data.plot$Tax)

other.df <- data.frame(
    Tax = "Other",
    Group = levels(factor(data.plot$Group)),
    mean = aggregate( formula = mean ~ Group, data = data.plot, FUN = function(x) 1 - sum(x) )[,2],
    sd = NA, median = NA, mad = NA
)
if ( sum(other.df$mean) > 0 ) {
    data.plot <- rbind( data.plot, other.df)
}

data.plot$Tax <- gsub(pattern = "_", replacement = " ", x = data.plot$Tax)
data.plot$Group <- gsub(pattern = "Forearm", replacement = "Unwounded", x = data.plot$Group)
data.plot$Group <- gsub(pattern = ":interm", replacement = "\nRDEB  intermediate", x = data.plot$Group)
data.plot$Group <- gsub(pattern = ":severe",  replacement = "\nRDEB  severe",   x = data.plot$Group)

fig2_1 <- ggplot(data.plot, aes(x=Group, y=mean, fill=Tax) ) +
    facet_wrap(~., scales="free_x") +
    geom_bar(width=0.75, stat='identity' ) +
    theme_bw() +
    ylim(c(0,1)) +
    theme(axis.line = element_line(colour = "black"),
          legend.text = element_text(face = "italic"),
          legend.justification = "top",
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank(),
          axis.text.x = element_text(angle = 45, hjust = 1, size=12),
          axis.text.y = element_text(size=12),
          strip.text.x = element_text(angle = 0, size = 10, color = "white"),
          axis.ticks.x=element_blank(),
          strip.background = element_rect(colour="white", fill="white", size=1.5, linetype="solid")) +
    ggtitle( paste("", "", sep=" ") ) +
    xlab("") +
    ylab("Proportion mapped contigs") +
    guides(fill = guide_legend(title = taxRank)) +
    scale_fill_manual(values=cols_top_phyla[names(cols_top_phyla) %in% data.plot$Tax])
fig2_1
ggsave(filename = paste0( outdir, "/Fig2_all.", taxRank,".pdf"), width = 4.5, height = 9)

taxRank <- "Genus"
# get abundance data
data <- getTaxaAbund(taxRank, ps_ra, "Location_Severity")
# sort taxonomic levels alphabetically
data$Tax <- factor(data$Tax, levels = sort(levels(data$Tax)))
# identify most abundant taxa across all samples
data.aggr <- aggregate( formula = mean ~ Tax, data = data, FUN = mean)
# select n most abundant taxa
best.tax <- data.aggr[ order(-data.aggr$mean), ][1:taxNumberPlot, 1]
# create new data frame
data.plot <- data[ data$Tax %in% best.tax, ]
data.plot <- data.plot[ data.plot$Tax != "NA.",]

levels(data.plot$Tax)

other.df <- data.frame(
    Tax = "Other",
    Group = levels(factor(data.plot$Group)),
    mean = aggregate( formula = mean ~ Group, data = data.plot, FUN = function(x) 1 - sum(x) )[,2],
    sd = NA, median = NA, mad = NA
)
if ( sum(other.df$mean) > 0 ) {
    data.plot <- rbind( data.plot, other.df)
}

data.plot$Tax <- gsub(pattern = "_", replacement = " ", x = data.plot$Tax)
data.plot$Group <- gsub(pattern = "Forearm", replacement = "Unwounded", x = data.plot$Group)
data.plot$Group <- gsub(pattern = ":interm", replacement = "\nRDEB  intermediate", x = data.plot$Group)
data.plot$Group <- gsub(pattern = ":severe",  replacement = "\nRDEB  severe",   x = data.plot$Group)

fig2_2 <- ggplot(data.plot, aes(factor(Tax), x=Group, y=mean, fill=Tax) ) +
    facet_wrap(~., scales="free_x") +
    geom_bar(width=0.75, stat='identity' ) +
    theme_bw() +
    ylim(c(0,1)) +
    theme(axis.line = element_line(colour = "black"),
          legend.text = element_text(face = "italic"),
          legend.justification = "top",
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank(),
          axis.text.x = element_text(angle = 45, hjust = 1, size=12),
          axis.text.y = element_text(size=12),
          strip.text.x = element_text(angle = 0, size = 10, color = "white"),
          axis.ticks.x=element_blank(),
          strip.background = element_rect(colour="white", fill="white", size=1.5, linetype="solid")) +
    ggtitle( paste("", "", sep=" ") ) +
    xlab("") +
    ylab("Proportion mapped contigs") +
    guides(fill = guide_legend(title = taxRank)) +
    scale_fill_manual(values=cols_top_genera[names(cols_top_genera) %in% data.plot$Tax])
fig2_2
ggsave(filename = paste0( outdir, "/Fig2_all.", taxRank,".pdf"), width = 4.5, height = 9)

#
# Alpha
#

alpha <- ps %>%
    sample_data() %>%
    data.frame() %>%
    dplyr::select( MINIONID, Location, Case, Location_Severity, Severity, Shannon ) %>%
    dplyr::filter( Location == "Forearm" | Location == "Wound" ) %>%
    dplyr::filter( Severity == "severe" | Severity == "interm" ) %>%
    dplyr::mutate( Location_Severity = as.character(Location_Severity) )

alpha$Location_Severity <- gsub(pattern = "Forearm", replacement = "Unwounded", x = alpha$Location_Severity)
alpha$Location_Severity <- gsub(pattern = ":interm", replacement = "\nRDEB intermediate", x = alpha$Location_Severity)
alpha$Location_Severity <- gsub(pattern = ":severe",  replacement = "\nRDEB severe",   x = alpha$Location_Severity)
alpha$Location_Severity <- factor(alpha$Location_Severity,
                                  levels = c("Unwounded\nRDEB intermediate",
                                             "Unwounded\nRDEB severe",
                                             "Wound\nRDEB intermediate",
                                             "Wound\nRDEB severe") )

mycomp <- list(
    c("Unwounded\nRDEB intermediate", "Unwounded\nRDEB severe"),
    c("Wound\nRDEB intermediate", "Wound\nRDEB severe")
)

fig2_3 <- ggpubr::ggdotplot(data = alpha , x = "Location_Severity", y = "Shannon",
                  fill = "Severity", palette = cols_serverity,
                  size = 0.75,
                  xlab = "", ylab = "Shannon index of species diversity") +
    theme(legend.position = "none",
          axis.text.x=element_text(size=12, angle=45, hjust = 1),
          axis.text.y=element_text(size=12, angle=0),
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          axis.line = element_line(colour = "black")) +
    ylim(0,7) +
    #ggpubr::stat_compare_means(method = "kruskal.test", label.y = 8) +
    stat_compare_means(comparisons = mycomp, label.y = c(6, 6))
fig2_3
ggsave(filename = paste0( outdir, "/Fig2_all.Alpha.pdf"), width = 4.5, height = 9)

#
# Correlations
#

# Phyla
corr_df <- phylosmith::variable_correlation( ps_ra %>% subset_samples(., !is.na(EBDASI.total)),
                                             treatment = c('Location'),
                                             subset = c('Forearm', 'Wound'),
                                             variables = 'EBDASI.total',
                                             classification = "Phylum", method = 'spearman')

corr_df$Treatment <- gsub(pattern = "Forearm", replacement = "Unwounded", x = corr_df$Treatment)

corr_df$p_clean <- corr_df$p
corr_df$p_clean <- format( round(corr_df$p_clean, 4), scientific = F ) # avoid exponentials
corr_df$p_clean[corr_df$p_clean >= 0.05] <- ""

fig2_4 <- corr_df %>%
    dplyr::mutate(X = gsub(pattern = "_", replacement = " ", X)) %>%
    ggplot(data = ., aes(x=Treatment, y=factor(X, levels = rev(levels(factor(X)))), fill=rho)) +
    geom_tile() +
    geom_text(aes(Treatment, X,
                  label = p_clean),
              color = "black", size = 4) +
    scale_fill_gradient2(low = "#2C7BB6", high = "#D7191C", mid = "white",
                         midpoint = 0, limit = c(-1,1), space = "Lab",
                         name="Spearman\nCorrelation") +
    theme(axis.line = element_line(colour = "black"),
          legend.position = "right",
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank(),
          axis.text.x = element_text(angle = 45, hjust = 1, size=12),
          strip.text.x = element_text(angle = 0, size = 12),
          axis.text.y = element_text(size=12, face = "italic"),
          axis.ticks.x=element_blank(),
          strip.background = element_rect(colour="white", fill="white", size=1.5, linetype="solid")) +
    xlab("") + ylab("")
fig2_4
ggsave(filename = paste0( outdir, "/Fig2_all.Correlation_EBDASI.total_Phylum.pdf"), width = 4.5, height = 9)

corr_df <- phylosmith::variable_correlation( ps_ra %>% subset_samples(., !is.na(BSA_woundstotal)),
                                             treatment = c('Location'),
                                             subset = c('Forearm', 'Wound'),
                                             variables = 'BSA_woundstotal',
                                             classification = "Phylum", method = 'spearman')

corr_df$Treatment <- gsub(pattern = "Forearm", replacement = "Unwounded", x = corr_df$Treatment)

corr_df$p_clean <- corr_df$p
corr_df$p_clean <- format( round(corr_df$p_clean, 4), scientific = F ) # avoid exponentials
corr_df$p_clean[corr_df$p_clean >= 0.05] <- ""

fig2_5 <- corr_df %>%
    dplyr::mutate(X = gsub(pattern = "_", replacement = " ", X)) %>%
    ggplot(data = ., aes(x=Treatment, y=factor(X, levels = rev(levels(factor(X)))), fill=rho)) +
    geom_tile() +
    geom_text(aes(Treatment, X,
                  label = p_clean),
              color = "black", size = 4) +
    scale_fill_gradient2(low = "#2C7BB6", high = "#D7191C", mid = "white",
                         midpoint = 0, limit = c(-1,1), space = "Lab",
                         name="Spearman\nCorrelation") +
    theme(axis.line = element_line(colour = "black"),
          legend.position = "right",
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank(),
          axis.text.x = element_text(angle = 45, hjust = 1, size=12),
          strip.text.x = element_text(angle = 0, size = 12),
          axis.text.y = element_text(size=12, face = "italic"),
          axis.ticks.x=element_blank(),
          strip.background = element_rect(colour="white", fill="white", size=1.5, linetype="solid")) +
    xlab("") + ylab("")
fig2_5
ggsave(filename = paste0( outdir, "/Fig2_all.Correlation_BSA_woundstotal_Phylum.pdf"), width = 4.5, height = 9)

# Genera
( ps.sel <- microbiome::aggregate_taxa(x = ps_ra, level = "Genus") %>% phylosmith::taxa_filter(., treatment = c('Location'), frequency = 0.2) )
corr_df <- phylosmith::variable_correlation( ps.sel %>% subset_samples(., !is.na(EBDASI.total)),
                                             treatment = c('Location'),
                                             subset = c('Forearm', 'Wound'),
                                             variables = 'EBDASI.total',
                                             classification = "Genus", method = 'spearman')

corr_df$Treatment <- gsub(pattern = "Forearm", replacement = "Unwounded", x = corr_df$Treatment)

corr_df$p_clean <- corr_df$p
corr_df$p_clean <- format( round(corr_df$p_clean, 4), scientific = F ) # avoid exponentials
corr_df$p_clean[corr_df$p_clean >= 0.05] <- ""

# keep only genera with p < 0.05
corr_genera <- corr_df$X[corr_df$p < 0.05]

fig2_6 <- corr_df %>%
    dplyr::filter(X %in% corr_genera) %>%
    dplyr::mutate(X = gsub(pattern = "_", replacement = " ", X)) %>%
    dplyr::filter(X != "Unknown") %>%
    ggplot(data = ., aes(x=Treatment, y=factor(X, levels = rev(levels(factor(X)))), fill=rho)) +
    geom_tile() +
    geom_text(aes(Treatment, X,
                  label = p_clean),
              color = "black", size = 4) +
    scale_fill_gradient2(low = "#2C7BB6", high = "#D7191C", mid = "white",
                         midpoint = 0, limit = c(-1,1), space = "Lab",
                         name="Spearman\nCorrelation") +
    theme(axis.line = element_line(colour = "black"),
          legend.position = "right",
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank(),
          axis.text.x = element_text(angle = 45, hjust = 1, size=12),
          strip.text.x = element_text(angle = 0, size = 12),
          axis.text.y = element_text(size=12, face = "italic"),
          axis.ticks.x=element_blank(),
          strip.background = element_rect(colour="white", fill="white", size=1.5, linetype="solid")) +
    xlab("") + ylab("")
fig2_6
ggsave(filename = paste0( outdir, "/Fig2_all.Correlation_EBDASI.total_Genus.pdf"), width = 4.5, height = 9)

corr_df <- phylosmith::variable_correlation( ps.sel %>% subset_samples(., !is.na(BSA_woundstotal)),
                                             treatment = c('Location'),
                                             subset = c('Forearm', 'Wound'),
                                             variables = 'BSA_woundstotal',
                                             classification = "Genus", method = 'spearman')

corr_df$Treatment <- gsub(pattern = "Forearm", replacement = "Unwounded", x = corr_df$Treatment)

corr_df$p_clean <- corr_df$p
corr_df$p_clean <- format( round(corr_df$p_clean, 4), scientific = F ) # avoid exponentials
corr_df$p_clean[corr_df$p_clean >= 0.05] <- ""

# keep only genera with p < 0.05
corr_genera <- corr_df$X[corr_df$p < 0.05]

fig2_7 <- corr_df %>%
    dplyr::filter(X %in% corr_genera) %>%
    dplyr::mutate(X = gsub(pattern = "_", replacement = " ", X)) %>%
    dplyr::filter(X != "Unknown") %>%
    ggplot(data = ., aes(x=Treatment, y=factor(X, levels = rev(levels(factor(X)))), fill=rho)) +
    geom_tile() +
    geom_text(aes(Treatment, X,
                  label = p_clean),
              color = "black", size = 4) +
    scale_fill_gradient2(low = "#2C7BB6", high = "#D7191C", mid = "white",
                         midpoint = 0, limit = c(-1,1), space = "Lab",
                         name="Spearman\nCorrelation") +
    theme(axis.line = element_line(colour = "black"),
          legend.position = "right",
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank(),
          axis.text.x = element_text(angle = 45, hjust = 1, size=12),
          strip.text.x = element_text(angle = 0, size = 12),
          axis.text.y = element_text(size=12, face = "italic"),
          axis.ticks.x=element_blank(),
          strip.background = element_rect(colour="white", fill="white", size=1.5, linetype="solid")) +
    xlab("") + ylab("")
fig2_7
ggsave(filename = paste0( outdir, "/Fig2_all.Correlation_BSA_woundstotal_Genus.pdf"), width = 4.5, height = 9)


# Figure 3 ----------------------------------------------------------------

ps_ra <- transform_sample_counts(ps, function(x){x / sum(x)})
ps_ra <- subset_samples(physeq = ps_ra, Case == "TRUE")

#
# Correlations
#

#
# CRP, Hb, Leukozytes
#

# Phyla
corr_df_CRP <- phylosmith::variable_correlation( ps_ra %>% subset_samples(., !is.na(CRP)),
                                                 treatment = c('Location'),
                                                 variables = c('CRP'),
                                                 classification = "Phylum", method = 'spearman')
corr_df_Hb <- phylosmith::variable_correlation( ps_ra %>% subset_samples(., !is.na(Hb)),
                                                treatment = c('Location'),
                                                variables = c('Hb'),
                                                classification = "Phylum", method = 'spearman')
corr_df_Leu <- phylosmith::variable_correlation( ps_ra %>% subset_samples(., !is.na(Leukozytes)),
                                                 treatment = c('Location'),
                                                 variables = c('Leukozytes'),
                                                 classification = "Phylum", method = 'spearman')
corr_df <- rbind(corr_df_CRP, corr_df_Hb, corr_df_Leu)

corr_df$Treatment <- gsub(pattern = "Forearm", replacement = "Unwounded", x = corr_df$Treatment)

corr_df$p_clean <- corr_df$p
corr_df$p_clean <- format( round(corr_df$p_clean, 4), scientific = F ) # avoid exponentials
corr_df$p_clean[corr_df$p_clean >= 0.05] <- ""

fig3_1 <- corr_df %>%
    dplyr::mutate(Treatment = factor(x = Treatment,
                                     levels = c("Unwounded", "Wound", "Oral mucosa", "Stool"))) %>%
    dplyr::mutate(X = gsub(pattern = "_", replacement = " ", X)) %>%
    dplyr::mutate(Y = gsub(pattern = "Leukozytes", replacement = "Leucozytes", Y)) %>%
    ggplot(data = ., aes(x=Treatment, y=factor(X, levels = rev(levels(factor(X)))), fill=rho)) +
    geom_tile()  +
    facet_wrap( Y ~ .) +
    geom_text(aes(Treatment, X,
                  label = p_clean),
              color = "black", size = 4) +
    scale_fill_gradient2(low = "#2C7BB6", high = "#D7191C", mid = "white",
                         midpoint = 0, limit = c(-1,1), space = "Lab",
                         name="Spearman\nCorrelation") +
    theme(axis.line = element_line(colour = "black"),
          legend.position = "right",
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank(),
          axis.text.x = element_text(angle = 45, hjust = 1, size=12),
          strip.text.x = element_text(angle = 0, size = 12),
          axis.text.y = element_text(size=12, face = "italic"),
          axis.ticks.x=element_blank(),
          strip.background = element_rect(colour="white", fill="white", size=1.5, linetype="solid")) +
    xlab("") + ylab("")
fig3_1
ggsave(filename = paste0( outdir, "/Fig3_all.Correlation_CRP_Hb_Leu_Phylum.pdf"), width = 12, height = 9)


# Genera
( ps.sel <- microbiome::aggregate_taxa(x = ps_ra, level = "Genus") %>% phylosmith::taxa_filter(., treatment = c('Location'), frequency = 0.2) )

corr_df_CRP <- phylosmith::variable_correlation( ps.sel %>% subset_samples(., !is.na(CRP)),
                                                 treatment = c('Location'),
                                                 variables = c('CRP'),
                                                 classification = "Genus", method = 'spearman')
corr_df_Hb <- phylosmith::variable_correlation( ps.sel %>% subset_samples(., !is.na(Hb)),
                                                treatment = c('Location'),
                                                variables = c('Hb'),
                                                classification = "Genus", method = 'spearman')
corr_df_Leu <- phylosmith::variable_correlation( ps.sel %>% subset_samples(., !is.na(Leukozytes)),
                                                 treatment = c('Location'),
                                                 variables = c('Leukozytes'),
                                                 classification = "Genus", method = 'spearman')
corr_df <- rbind(corr_df_CRP, corr_df_Hb, corr_df_Leu)

corr_df$Treatment <- gsub(pattern = "Forearm", replacement = "Unwounded", x = corr_df$Treatment)

corr_df$p_clean <- corr_df$p
corr_df$p_clean <- format( round(corr_df$p_clean, 4), scientific = F ) # avoid exponentials
corr_df$p_clean[corr_df$p_clean >= 0.05] <- ""

# keep only genera with p < 0.05
corr_genera <- corr_df$X[corr_df$p < 0.05]

fig3_2 <- corr_df %>%
    dplyr::mutate(Treatment = factor(x = Treatment,
                                     levels = c("Unwounded", "Wound", "Oral mucosa", "Stool"))) %>%
    dplyr::filter(X %in% corr_genera) %>%
    dplyr::mutate(X = gsub(pattern = "_", replacement = " ", X)) %>%
    dplyr::mutate(Y = gsub(pattern = "Leukozytes", replacement = "Leucozytes", Y)) %>%
    dplyr::filter(X != "Unknown") %>%
    dplyr::filter(X != "Unclassified") %>%
    ggplot(data = ., aes(x=Treatment, y=factor(X, levels = rev(levels(factor(X)))), fill=rho)) +
    geom_tile() +
    facet_wrap( Y ~ .) +
    geom_text(aes(Treatment, X,
                  label = p_clean),
              color = "black", size = 4) +
    scale_fill_gradient2(low = "#2C7BB6", high = "#D7191C", mid = "white",
                         midpoint = 0, limit = c(-1,1), space = "Lab",
                         name="Spearman\nCorrelation") +
    theme(axis.line = element_line(colour = "black"),
          legend.position = "right",
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank(),
          axis.text.x = element_text(angle = 45, hjust = 1, size=12),
          strip.text.x = element_text(angle = 0, size = 12),
          axis.text.y = element_text(size=12, face = "italic"),
          axis.ticks.x=element_blank(),
          strip.background = element_rect(colour="white", fill="white", size=1.5, linetype="solid")) +
    xlab("") + ylab("")
fig3_2
ggsave(filename = paste0( outdir, "/Fig3_all.Correlation_CRP_Hb_Leu_Genus.pdf"), width = 12, height = 9)


#
# Immunoglobulins
#

# Phyla
corr_df_IgA <- phylosmith::variable_correlation( ps_ra %>% subset_samples(., !is.na(IgA)),
                                                 treatment = c('Location'),
                                                 variables = c('IgA'),
                                                 classification = "Phylum", method = 'spearman')
corr_df_IgE <- phylosmith::variable_correlation( ps_ra %>% subset_samples(., !is.na(IgE)),
                                                 treatment = c('Location'),
                                                 variables = c('IgE'),
                                                 classification = "Phylum", method = 'spearman')
corr_df_IgG <- phylosmith::variable_correlation( ps_ra %>% subset_samples(., !is.na(IgG)),
                                                 treatment = c('Location'),
                                                 variables = c('IgG'),
                                                 classification = "Phylum", method = 'spearman')
corr_df_IgM <- phylosmith::variable_correlation( ps_ra %>% subset_samples(., !is.na(IgG)),
                                                 treatment = c('Location'),
                                                 variables = c('IgM'),
                                                 classification = "Phylum", method = 'spearman')
corr_df <- rbind(corr_df_IgA, corr_df_IgE, corr_df_IgG, corr_df_IgM)

corr_df$Treatment <- gsub(pattern = "Forearm", replacement = "Unwounded", x = corr_df$Treatment)

corr_df$p_clean <- corr_df$p
corr_df$p_clean <- format( round(corr_df$p_clean, 4), scientific = F ) # avoid exponentials
corr_df$p_clean[corr_df$p_clean >= 0.05] <- ""

fig3_3 <- corr_df %>%
    dplyr::mutate(Treatment = factor(x = Treatment,
                                     levels = c("Unwounded", "Wound", "Oral mucosa", "Stool"))) %>%
    dplyr::mutate(X = gsub(pattern = "_", replacement = " ", X)) %>%
    ggplot(data = ., aes(x=Treatment, y=factor(X, levels = rev(levels(factor(X)))), fill=rho)) +
    geom_tile()  +
    facet_wrap( Y ~ ., ncol = 4) +
    geom_text(aes(Treatment, X,
                  label = p_clean),
              color = "black", size = 4) +
    scale_fill_gradient2(low = "#2C7BB6", high = "#D7191C", mid = "white",
                         midpoint = 0, limit = c(-1,1), space = "Lab",
                         name="Spearman\nCorrelation") +
    theme(axis.line = element_line(colour = "black"),
          legend.position = "right",
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank(),
          axis.text.x = element_text(angle = 45, hjust = 1, size=12),
          strip.text.x = element_text(angle = 0, size = 12),
          axis.text.y = element_text(size=12, face = "italic"),
          axis.ticks.x=element_blank(),
          strip.background = element_rect(colour="white", fill="white", size=1.5, linetype="solid")) +
    xlab("") + ylab("")
fig3_3
ggsave(filename = paste0( outdir, "/Fig3_all.Correlation_IgX_Phylum.pdf"), width = 12, height = 9)


# Genera
( ps.sel <- microbiome::aggregate_taxa(x = ps_ra, level = "Genus") %>% phylosmith::taxa_filter(., treatment = c('Location'), frequency = 0.2) )

corr_df_IgA <- phylosmith::variable_correlation( ps.sel %>% subset_samples(., !is.na(IgA)),
                                                 treatment = c('Location'),
                                                 variables = c('IgA'),
                                                 classification = "Genus", method = 'spearman')
corr_df_IgE <- phylosmith::variable_correlation( ps.sel %>% subset_samples(., !is.na(IgE)),
                                                 treatment = c('Location'),
                                                 variables = c('IgE'),
                                                 classification = "Genus", method = 'spearman')
corr_df_IgG <- phylosmith::variable_correlation( ps.sel %>% subset_samples(., !is.na(IgG)),
                                                 treatment = c('Location'),
                                                 variables = c('IgG'),
                                                 classification = "Genus", method = 'spearman')
corr_df_IgM <- phylosmith::variable_correlation( ps.sel %>% subset_samples(., !is.na(IgG)),
                                                 treatment = c('Location'),
                                                 variables = c('IgM'),
                                                 classification = "Genus", method = 'spearman')
corr_df <- rbind(corr_df_IgA, corr_df_IgE, corr_df_IgG, corr_df_IgM)

corr_df$Treatment <- gsub(pattern = "Forearm", replacement = "Unwounded", x = corr_df$Treatment)

corr_df$p_clean <- corr_df$p
corr_df$p_clean <- format( round(corr_df$p_clean, 4), scientific = F ) # avoid exponentials
corr_df$p_clean[corr_df$p_clean >= 0.05] <- ""

# keep only genera with p < 0.05
corr_genera <- corr_df$X[corr_df$p < 0.05]

fig3_4 <- corr_df %>%
    dplyr::mutate(Treatment = factor(x = Treatment,
                                     levels = c("Unwounded", "Wound", "Oral mucosa", "Stool"))) %>%
    dplyr::filter(X %in% corr_genera) %>%
    dplyr::mutate(X = gsub(pattern = "_", replacement = " ", X)) %>%
    dplyr::filter(X != "Unknown") %>%
    dplyr::filter(X != "Unclassified") %>%
    ggplot(data = ., aes(x=Treatment, y=factor(X, levels = rev(levels(factor(X)))), fill=rho)) +
    geom_tile() +
    facet_wrap( Y ~ ., ncol = 4) +
    geom_text(aes(Treatment, X,
                  label = p_clean),
              color = "black", size = 4) +
    scale_fill_gradient2(low = "#2C7BB6", high = "#D7191C", mid = "white",
                         midpoint = 0, limit = c(-1,1), space = "Lab",
                         name="Spearman\nCorrelation") +
    theme(axis.line = element_line(colour = "black"),
          legend.position = "right",
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank(),
          axis.text.x = element_text(angle = 45, hjust = 1, size=12),
          strip.text.x = element_text(angle = 0, size = 12),
          axis.text.y = element_text(size=12, face = "italic"),
          axis.ticks.x=element_blank(),
          strip.background = element_rect(colour="white", fill="white", size=1.5, linetype="solid")) +
    xlab("") + ylab("")
fig3_4
ggsave(filename = paste0( outdir, "/Fig3_all.Correlation_IgX_Genus.pdf"), width = 12, height = 9)


#
# Cytokines
#

# Phyla
corr_df_IFN <- phylosmith::variable_correlation( ps_ra %>% subset_samples(., !is.na(IFN.g)),
                                                 treatment = c('Location'),
                                                 variables = c('IFN.g'),
                                                 classification = "Phylum", method = 'spearman')
corr_df_IL6 <- phylosmith::variable_correlation( ps_ra %>% subset_samples(., !is.na(IL.6)),
                                                 treatment = c('Location'),
                                                 variables = c('IL.6'),
                                                 classification = "Phylum", method = 'spearman')
corr_df_TGF <- phylosmith::variable_correlation( ps_ra %>% subset_samples(., !is.na(TGF.beta)),
                                                 treatment = c('Location'),
                                                 variables = c('TGF.beta'),
                                                 classification = "Phylum", method = 'spearman')
corr_df_TNF <- phylosmith::variable_correlation( ps_ra %>% subset_samples(., !is.na(TNF.a)),
                                                 treatment = c('Location'),
                                                 variables = c('TNF.a'),
                                                 classification = "Phylum", method = 'spearman')
corr_df <- rbind(corr_df_IFN, corr_df_IL6, corr_df_TGF, corr_df_TNF)

corr_df$Treatment <- gsub(pattern = "Forearm", replacement = "Unwounded", x = corr_df$Treatment)

corr_df$p_clean <- corr_df$p
corr_df$p_clean <- format( round(corr_df$p_clean, 4), scientific = F ) # avoid exponentials
corr_df$p_clean[corr_df$p_clean >= 0.05] <- ""

fig3_5 <- corr_df %>%
    dplyr::mutate(Treatment = factor(x = Treatment,
                                     levels = c("Unwounded", "Wound", "Oral mucosa", "Stool"))) %>%
    dplyr::mutate(X = gsub(pattern = "_", replacement = " ", X)) %>%
    dplyr::mutate(Y = gsub(pattern = "IFN.g", replacement = "IFN-gamma", Y)) %>%
    dplyr::mutate(Y = gsub(pattern = "IL.6", replacement = "IL-6", Y)) %>%
    dplyr::mutate(Y = gsub(pattern = "TGF.beta", replacement = "TGF-beta", Y)) %>%
    dplyr::mutate(Y = gsub(pattern = "TNF.a", replacement = "TNF-alpha", Y)) %>%
    ggplot(data = ., aes(x=Treatment, y=factor(X, levels = rev(levels(factor(X)))), fill=rho)) +
    geom_tile()  +
    facet_wrap( Y ~ ., ncol = 4, labeller = label_parsed) +
    geom_text(aes(Treatment, X,
                  label = p_clean),
              color = "black", size = 4) +
    scale_fill_gradient2(low = "#2C7BB6", high = "#D7191C", mid = "white",
                         midpoint = 0, limit = c(-1,1), space = "Lab",
                         name="Spearman\nCorrelation") +
    theme(axis.line = element_line(colour = "black"),
          legend.position = "right",
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank(),
          axis.text.x = element_text(angle = 45, hjust = 1, size=12),
          strip.text.x = element_text(angle = 0, size = 12),
          axis.text.y = element_text(size=12, face = "italic"),
          axis.ticks.x=element_blank(),
          strip.background = element_rect(colour="white", fill="white", size=1.5, linetype="solid")) +
    xlab("") + ylab("")
fig3_5
ggsave(filename = paste0( outdir, "/Fig3_all.Correlation_Cytokines_Phylum.pdf"), width = 12, height = 9)

# Genera
( ps.sel <- microbiome::aggregate_taxa(x = ps_ra, level = "Genus") %>% phylosmith::taxa_filter(., treatment = c('Location'), frequency = 0.2) )

corr_df_IFN <- phylosmith::variable_correlation( ps.sel %>% subset_samples(., !is.na(IFN.g)),
                                                 treatment = c('Location'),
                                                 variables = c('IFN.g'),
                                                 classification = "Genus", method = 'spearman')
corr_df_IL6 <- phylosmith::variable_correlation( ps.sel %>% subset_samples(., !is.na(IL.6)),
                                                 treatment = c('Location'),
                                                 variables = c('IL.6'),
                                                 classification = "Genus", method = 'spearman')
corr_df_TGF <- phylosmith::variable_correlation( ps.sel %>% subset_samples(., !is.na(TGF.beta)),
                                                 treatment = c('Location'),
                                                 variables = c('TGF.beta'),
                                                 classification = "Genus", method = 'spearman')
corr_df_TNF <- phylosmith::variable_correlation( ps.sel %>% subset_samples(., !is.na(TNF.a)),
                                                 treatment = c('Location'),
                                                 variables = c('TNF.a'),
                                                 classification = "Genus", method = 'spearman')
corr_df <- rbind(corr_df_IFN, corr_df_IL6, corr_df_TGF, corr_df_TNF)

corr_df$Treatment <- gsub(pattern = "Forearm", replacement = "Unwounded", x = corr_df$Treatment)

corr_df$p_clean <- corr_df$p
corr_df$p_clean <- format( round(corr_df$p_clean, 4), scientific = F ) # avoid exponentials
corr_df$p_clean[corr_df$p_clean >= 0.05] <- ""

# keep only genera with p < 0.05
corr_genera <- corr_df$X[corr_df$p < 0.05]

fig3_6 <- corr_df %>%
    dplyr::mutate(Treatment = factor(x = Treatment,
                                     levels = c("Unwounded", "Wound", "Oral mucosa", "Stool"))) %>%
    dplyr::filter(X %in% corr_genera) %>%
    dplyr::mutate(X = gsub(pattern = "_", replacement = " ", X)) %>%
    dplyr::mutate(Y = gsub(pattern = "IFN.g", replacement = "IFN-gamma", Y)) %>%
    dplyr::mutate(Y = gsub(pattern = "IL.6", replacement = "IL-6", Y)) %>%
    dplyr::mutate(Y = gsub(pattern = "TGF.beta", replacement = "TGF-beta", Y)) %>%
    dplyr::mutate(Y = gsub(pattern = "TNF.a", replacement = "TNF-alpha", Y)) %>%
    dplyr::filter(X != "Unknown") %>%
    dplyr::filter(X != "Unclassified") %>%
    ggplot(data = ., aes(x=Treatment, y=factor(X, levels = rev(levels(factor(X)))), fill=rho)) +
    geom_tile() +
    facet_wrap( Y ~ ., ncol = 4, labeller = label_parsed) +
    geom_text(aes(Treatment, X,
                  label = p_clean),
              color = "black", size = 4) +
    scale_fill_gradient2(low = "#2C7BB6", high = "#D7191C", mid = "white",
                         midpoint = 0, limit = c(-1,1), space = "Lab",
                         name="Spearman\nCorrelation") +
    theme(axis.line = element_line(colour = "black"),
          legend.position = "right",
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank(),
          axis.text.x = element_text(angle = 45, hjust = 1, size=12),
          strip.text.x = element_text(angle = 0, size = 12),
          axis.text.y = element_text(size=12, face = "italic"),
          axis.ticks.x=element_blank(),
          strip.background = element_rect(colour="white", fill="white", size=1.5, linetype="solid")) +
    xlab("") + ylab("")
fig3_6
ggsave(filename = paste0( outdir, "/Fig3_all.Correlation_Cytokines_Genus.pdf"), width = 12, height = 9)


#
# Model Forearm/Wound vs Staph aureus
#

ps_ra <- transform_sample_counts(ps, function(x){x / sum(x)})

cova <- ps %>%
    sample_data %>%
    data.frame %>%
    dplyr::filter( Case == "TRUE" ) %>%
    dplyr::filter( Location == "Forearm" | Location == "Wound" ) %>%
    dplyr::select(MINIONID,
                  Location, Age, Gender, Weight = Weighjt,
                  BSA_woundstotal, BSA_wounds3weeks,
                  EBDASItotal = EBDASI.total,
                  CRP,
                  Hb, Leukozytes,
                  IFNg = IFN.g, IL6 = IL.6, TGFbeta = TGF.beta, TNFa = TNF.a # Cytokines
    )

ra_sau <- microbiome::aggregate_taxa(x = ps_ra, level = "Species") %>%
    subset_taxa(., Species=="S.aureus") %>%
    otu_table() %>%
    data.frame() %>%
    t()

df_s <- merge(cova, ra_sau, by = 0)

rm(ra_sau); rm(cova); rm(ps_ra)

# Forearm

LOCATION <- "Forearm"
df_s.tmp <- df_s %>%
    dplyr::filter(Location == LOCATION) %>%
    dplyr::rename( S_aureus = S.aureus ) %>%
    dplyr::select( -Row.names, -MINIONID, -Location, -Age, -Gender, -Weight )

lm(formula = S_aureus ~ EBDASItotal, data = df_s.tmp ) %>% summary()

m_sw <- lm(formula = log(S_aureus) ~ ., data = df_s.tmp )
performance::model_performance(m_sw)
summary(m_sw)

# stepwise regression
selectedMod <- step(m_sw)
summary(selectedMod)
car::vif(selectedMod)

sjPlot::tab_model(selectedMod, digits = 4, digits.p = 4)
effectsize::eta_squared(car::Anova(selectedMod, type = 2), partial = FALSE, ci = 0.95)

fig3_7 <- ggplot(df_s.tmp, aes(EBDASItotal, log(S_aureus)) ) +
    geom_point() +
    geom_smooth(method='lm') +
    theme(axis.line = element_line(colour = "black"),
          legend.position = "right",
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank(),
          axis.text.x = element_text(angle = 0, size=12),
          strip.text.x = element_text(angle = 0, size = 12),
          axis.text.y = element_text(size=12),
          strip.background = element_rect(colour="white", fill="white", size=1.5, linetype="solid")) +
    labs(x='EBDASI total', y='log10(Rel. abundance S. aureus)', title='') +
    #ylim(-0.5,1) +
    geom_hline(yintercept=0, linetype="dashed",
               color = "grey45", size = 1)
fig3_7

ggsave( filename = paste0("plots_paper/Fig3.Forearm_Disease_Model.pdf"), width = 6, height = 6 )

# Wound

LOCATION <- "Wound"
df_s.tmp <- df_s %>%
    dplyr::filter(Location == LOCATION) %>%
    dplyr::rename( S_aureus = S.aureus ) %>%
    dplyr::select( -Row.names, -MINIONID, -Location, -Age, -Gender, -Weight )

m_sw <- lm(formula = log10(S_aureus) ~ ., data = df_s.tmp )
performance::model_performance(m_sw)
summary(m_sw)

# stepwise regression
selectedMod <- step(m_sw)
summary(selectedMod)
car::vif(selectedMod)

sjPlot::tab_model(selectedMod, digits = 4, digits.p = 4)
effectsize::eta_squared(car::Anova(selectedMod, type = 2), partial = FALSE, ci = 0.95)

# Figure 4 ----------------------------------------------------------------

ps_ra <- transform_sample_counts(ps, function(x){x / sum(x)})
ps_ra <- subset_samples(physeq = ps_ra, Location == "Oral mucosa" | Location == "Stool")

#
# POOLED SAMPLES
#

taxRank <- "Phylum"
# get abundance data
data <- getTaxaAbund(taxRank, ps_ra, "Location_Case")
# sort taxonomic levels alphabetically
data$Tax <- factor(data$Tax, levels = sort(levels(data$Tax)))
# identify most abundant taxa across all samples
data.aggr <- aggregate( formula = mean ~ Tax, data = data, FUN = mean)
# select n most abundant taxa
best.tax <- data.aggr[ order(-data.aggr$mean), ][1:taxNumberPlot, 1]
# create new data frame
data.plot <- data[ data$Tax %in% best.tax, ]
data.plot <- data.plot[ data.plot$Tax != "NA.",]

levels(data.plot$Tax)

other.df <- data.frame(
    Tax = "Other",
    Group = levels(factor(data.plot$Group)),
    mean = aggregate( formula = mean ~ Group, data = data.plot, FUN = function(x) 1 - sum(x) )[,2],
    sd = NA, median = NA, mad = NA
)
if ( sum(other.df$mean) > 0 ) {
    data.plot <- rbind( data.plot, other.df)
}

data.plot$Tax <- gsub(pattern = "_", replacement = " ", x = data.plot$Tax)
data.plot$Group <- gsub(pattern = ":FALSE", replacement = " CTRL", x = data.plot$Group)
data.plot$Group <- gsub(pattern = ":TRUE",  replacement = " RDEB",   x = data.plot$Group)

fig4_1 <- ggplot(data.plot, aes(x=Group, y=mean, fill=Tax) ) +
    facet_wrap(~., scales="free_x") +
    geom_bar(width=0.75, stat='identity' ) +
    theme_bw() +
    ylim(c(0,1)) +
    theme(axis.line = element_line(colour = "black"),
          legend.text = element_text(face = "italic"),
          legend.justification = "top",
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank(),
          axis.text.x = element_text(angle = 45, hjust = 1, size=11),
          axis.text.y = element_text(size=12),
          strip.text.x = element_text(angle = 0, size = 10, color = "white"),
          axis.ticks.x=element_blank(),
          strip.background = element_rect(colour="white", fill="white", size=1.5, linetype="solid")) +
    ggtitle( paste("", "", sep=" ") ) +
    xlab("") +
    ylab("Proportion mapped contigs") +
    guides(fill = guide_legend(title = taxRank)) +
    scale_fill_manual(values=cols_top_phyla[names(cols_top_phyla) %in% data.plot$Tax])
fig4_1
ggsave(filename = paste0( outdir, "/Fig4_all.", taxRank,".pdf"), width = 4.5, height = 9)

taxRank <- "Genus"
# get abundance data
data <- getTaxaAbund(taxRank, ps_ra, "Location_Case")
# sort taxonomic levels alphabetically
data$Tax <- factor(data$Tax, levels = sort(levels(data$Tax)))
# identify most abundant taxa across all samples
data.aggr <- aggregate( formula = mean ~ Tax, data = data, FUN = mean)
# select n most abundant taxa
best.tax <- data.aggr[ order(-data.aggr$mean), ][1:taxNumberPlot, 1]
# create new data frame
data.plot <- data[ data$Tax %in% best.tax, ]
data.plot <- data.plot[ data.plot$Tax != "NA.",]

levels(data.plot$Tax)

other.df <- data.frame(
    Tax = "Other",
    Group = levels(factor(data.plot$Group)),
    mean = aggregate( formula = mean ~ Group, data = data.plot, FUN = function(x) 1 - sum(x) )[,2],
    sd = NA, median = NA, mad = NA
)
if ( sum(other.df$mean) > 0 ) {
    data.plot <- rbind( data.plot, other.df)
}

data.plot$Tax <- gsub(pattern = "_", replacement = " ", x = data.plot$Tax)
data.plot$Group <- gsub(pattern = ":FALSE", replacement = " CTRL", x = data.plot$Group)
data.plot$Group <- gsub(pattern = ":TRUE",  replacement = " RDEB",   x = data.plot$Group)

fig4_2 <- ggplot(data.plot, aes(factor(Tax), x=Group, y=mean, fill=Tax) ) +
    facet_wrap(~., scales="free_x") +
    geom_bar(width=0.75, stat='identity' ) +
    theme_bw() +
    ylim(c(0,1)) +
    theme(axis.line = element_line(colour = "black"),
          legend.text = element_text(face = "italic"),
          legend.justification = "top",
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank(),
          axis.text.x = element_text(angle = 45, hjust = 1, size=11),
          axis.text.y = element_text(size=12),
          strip.text.x = element_text(angle = 0, size = 10, color = "white"),
          axis.ticks.x=element_blank(),
          strip.background = element_rect(colour="white", fill="white", size=1.5, linetype="solid")) +
    ggtitle( paste("", "", sep=" ") ) +
    xlab("") +
    ylab("Proportion mapped contigs") +
    guides(fill = guide_legend(title = taxRank)) +
    scale_fill_manual(values=cols_top_genera[names(cols_top_genera) %in% data.plot$Tax])
fig4_2
ggsave(filename = paste0( outdir, "/Fig4_all.", taxRank,".pdf"), width = 4.5, height = 9)

#
# Age groups
#

taxRank <- "Phylum"
# get abundance data
data <- getTaxaAbund(taxRank, ps_ra, "Age_Case")
# sort taxonomic levels alphabetically
data$Tax <- factor(data$Tax, levels = sort(levels(data$Tax)))
# identify most abundant taxa across all samples
data.aggr <- aggregate( formula = mean ~ Tax, data = data, FUN = mean)
# select n most abundant taxa
best.tax <- data.aggr[ order(-data.aggr$mean), ][1:taxNumberPlot, 1]
# create new data frame
data.plot <- data[ data$Tax %in% best.tax, ]
data.plot <- data.plot[ data.plot$Tax != "NA.",]

levels(data.plot$Tax)

other.df <- data.frame(
    Tax = "Other",
    Group = levels(factor(data.plot$Group)),
    mean = aggregate( formula = mean ~ Group, data = data.plot, FUN = function(x) 1 - sum(x) )[,2],
    sd = NA, median = NA, mad = NA
)
if ( sum(other.df$mean) > 0 ) {
    data.plot <- rbind( data.plot, other.df)
}

data.plot$Tax <- gsub(pattern = "_", replacement = " ", x = data.plot$Tax)
data.plot$Group <- gsub(pattern = ":FALSE", replacement = " CTRL", x = data.plot$Group)
data.plot$Group <- gsub(pattern = ":TRUE",  replacement = " RDEB",   x = data.plot$Group)

data.plot$AgeGroup <- NA
data.plot$AgeGroup[ grep(pattern = "01:", data.plot$Group) ] <-  "0-4 years"
data.plot$AgeGroup[ grep(pattern = "02:", data.plot$Group) ] <-  "5-10 years"
data.plot$AgeGroup[ grep(pattern = "03:", data.plot$Group) ] <-  "11-22 years"

data.plot$AgeGroup <- factor(data.plot$AgeGroup, levels = c("0-4 years", "5-10 years", "11-22 years"))

data.plot$Group <- gsub(pattern = "01:", replacement = "", x = data.plot$Group)
data.plot$Group <- gsub(pattern = "02:", replacement = "", x = data.plot$Group)
data.plot$Group <- gsub(pattern = "03:", replacement = "", x = data.plot$Group)

fig4_3 <- ggplot(data.plot, aes(x=Group, y=mean, fill=Tax) ) +
    facet_wrap(AgeGroup~., scales="free_x") +
    geom_bar(width=0.75, stat='identity' ) +
    theme_bw() +
    ylim(c(0,1)) +
    theme(axis.line = element_line(colour = "black"),
          legend.text = element_text(face = "italic"),
          legend.justification = "top",
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank(),
          axis.text.x = element_text(angle = 45, hjust = 1, size=11),
          strip.text.x = element_text(angle = 0, size = 12),
          axis.text.y = element_text(size=12),
          axis.ticks.x=element_blank(),
          strip.background = element_rect(colour="white", fill="white", size=1.5, linetype="solid")) +
    ggtitle( paste("", "", sep=" ") ) +
    xlab("") +
    ylab("Proportion mapped contigs") +
    guides(fill = guide_legend(title = taxRank)) +
    scale_fill_manual(values=cols_top_phyla[names(cols_top_phyla) %in% data.plot$Tax])
fig4_3
ggsave(filename = paste0( outdir, "/Fig4_age.", taxRank,".pdf"), width = 9, height = 9)

taxRank <- "Genus"
# get abundance data
data <- getTaxaAbund(taxRank, ps_ra, "Age_Case")
# sort taxonomic levels alphabetically
data$Tax <- factor(data$Tax, levels = sort(levels(data$Tax)))
# identify most abundant taxa across all samples
data.aggr <- aggregate( formula = mean ~ Tax, data = data, FUN = mean)
# select n most abundant taxa
best.tax <- data.aggr[ order(-data.aggr$mean), ][1:taxNumberPlot, 1]
# create new data frame
data.plot <- data[ data$Tax %in% best.tax, ]
data.plot <- data.plot[ data.plot$Tax != "NA.",]

levels(data.plot$Tax)

other.df <- data.frame(
    Tax = "Other",
    Group = levels(factor(data.plot$Group)),
    mean = aggregate( formula = mean ~ Group, data = data.plot, FUN = function(x) 1 - sum(x) )[,2],
    sd = NA, median = NA, mad = NA
)
if ( sum(other.df$mean) > 0 ) {
    data.plot <- rbind( data.plot, other.df)
}

data.plot$Tax <- gsub(pattern = "_", replacement = " ", x = data.plot$Tax)
data.plot$Group <- gsub(pattern = ":FALSE", replacement = " CTRL", x = data.plot$Group)
data.plot$Group <- gsub(pattern = ":TRUE",  replacement = " RDEB",   x = data.plot$Group)

data.plot$AgeGroup <- NA
data.plot$AgeGroup[ grep(pattern = "01:", data.plot$Group) ] <-  "0-4 years"
data.plot$AgeGroup[ grep(pattern = "02:", data.plot$Group) ] <-  "5-10 years"
data.plot$AgeGroup[ grep(pattern = "03:", data.plot$Group) ] <-  "11-22 years"

data.plot$AgeGroup <- factor(data.plot$AgeGroup, levels = c("0-4 years", "5-10 years", "11-22 years"))

data.plot$Group <- gsub(pattern = "01:", replacement = "", x = data.plot$Group)
data.plot$Group <- gsub(pattern = "02:", replacement = "", x = data.plot$Group)
data.plot$Group <- gsub(pattern = "03:", replacement = "", x = data.plot$Group)

fig4_4 <- ggplot(data.plot, aes(factor(Tax), x=Group, y=mean, fill=Tax) ) +
    facet_wrap(AgeGroup~., scales="free_x") +
    geom_bar(width=0.75, stat='identity' ) +
    theme_bw() +
    ylim(c(0,1)) +
    theme(axis.line = element_line(colour = "black"),
          legend.text = element_text(face = "italic" ),
          legend.justification = "top",
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank(),
          axis.text.x = element_text(angle = 45, hjust = 1, size=11),
          strip.text.x = element_text(angle = 0, size = 12),
          axis.text.y = element_text(size=12),
          axis.ticks.x=element_blank(),
          strip.background = element_rect(colour="white", fill="white", size=1.5, linetype="solid")) +
    ggtitle( paste("", "", sep=" ") ) +
    xlab("") +
    ylab("Proportion mapped contigs") +
    guides(fill = guide_legend(title = taxRank, ncol = 1)) +
    scale_fill_manual(values=cols_top_genera[names(cols_top_genera) %in% data.plot$Tax])
fig4_4
ggsave(filename = paste0( outdir, "/Fig4_age.", taxRank,".pdf"), width = 9, height = 9)

#
# Severity
#

ps_ra <- transform_sample_counts(ps, function(x){x / sum(x)})
ps_ra <- subset_samples(physeq = ps_ra, Location == "Oral mucosa" | Location == "Stool")
ps_ra <- subset_samples(physeq = ps_ra, Severity == "severe" | Severity == "interm")

taxRank <- "Phylum"
# get abundance data
data <- getTaxaAbund(taxRank, ps_ra, "Location_Severity")
# sort taxonomic levels alphabetically
data$Tax <- factor(data$Tax, levels = sort(levels(data$Tax)))
# identify most abundant taxa across all samples
data.aggr <- aggregate( formula = mean ~ Tax, data = data, FUN = mean)
# select n most abundant taxa
best.tax <- data.aggr[ order(-data.aggr$mean), ][1:taxNumberPlot, 1]
# create new data frame
data.plot <- data[ data$Tax %in% best.tax, ]
data.plot <- data.plot[ data.plot$Tax != "NA.",]

levels(data.plot$Tax)

other.df <- data.frame(
    Tax = "Other",
    Group = levels(factor(data.plot$Group)),
    mean = aggregate( formula = mean ~ Group, data = data.plot, FUN = function(x) 1 - sum(x) )[,2],
    sd = NA, median = NA, mad = NA
)
if ( sum(other.df$mean) > 0 ) {
    data.plot <- rbind( data.plot, other.df)
}

data.plot$Tax <- gsub(pattern = "_", replacement = " ", x = data.plot$Tax)
data.plot$Group <- gsub(pattern = ":interm", replacement = "\nRDEB intermediate", x = data.plot$Group)
data.plot$Group <- gsub(pattern = ":severe",  replacement = "\nRDEB severe",   x = data.plot$Group)

fig4_5 <- ggplot(data.plot, aes(x=Group, y=mean, fill=Tax) ) +
    facet_wrap(~., scales="free_x") +
    geom_bar(width=0.75, stat='identity' ) +
    theme_bw() +
    ylim(c(0,1)) +
    theme(axis.line = element_line(colour = "black"),
          legend.text = element_text(face = "italic"),
          legend.justification = "top",
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank(),
          axis.text.x = element_text(angle = 45, hjust = 1, size=12),
          axis.text.y = element_text(size=12),
          strip.text.x = element_text(angle = 0, size = 10, color = "white"),
          axis.ticks.x=element_blank(),
          strip.background = element_rect(colour="white", fill="white", size=1.5, linetype="solid")) +
    ggtitle( paste("", "", sep=" ") ) +
    xlab("") +
    ylab("Proportion mapped contigs") +
    guides(fill = guide_legend(title = taxRank)) +
    scale_fill_manual(values=cols_top_phyla[names(cols_top_phyla) %in% data.plot$Tax])
fig4_5
ggsave(filename = paste0( outdir, "/Fig4_all.", taxRank,".severity.pdf"), width = 4.5, height = 9)

taxRank <- "Genus"
# get abundance data
data <- getTaxaAbund(taxRank, ps_ra, "Location_Severity")
# sort taxonomic levels alphabetically
data$Tax <- factor(data$Tax, levels = sort(levels(data$Tax)))
# identify most abundant taxa across all samples
data.aggr <- aggregate( formula = mean ~ Tax, data = data, FUN = mean)
# select n most abundant taxa
best.tax <- data.aggr[ order(-data.aggr$mean), ][1:taxNumberPlot, 1]
# create new data frame
data.plot <- data[ data$Tax %in% best.tax, ]
data.plot <- data.plot[ data.plot$Tax != "NA.",]

levels(data.plot$Tax)

other.df <- data.frame(
    Tax = "Other",
    Group = levels(factor(data.plot$Group)),
    mean = aggregate( formula = mean ~ Group, data = data.plot, FUN = function(x) 1 - sum(x) )[,2],
    sd = NA, median = NA, mad = NA
)
if ( sum(other.df$mean) > 0 ) {
    data.plot <- rbind( data.plot, other.df)
}

data.plot$Tax <- gsub(pattern = "_", replacement = " ", x = data.plot$Tax)
data.plot$Group <- gsub(pattern = ":interm", replacement = "\nRDEB intermediate", x = data.plot$Group)
data.plot$Group <- gsub(pattern = ":severe",  replacement = "\nRDEB severe",   x = data.plot$Group)

fig4_6 <- ggplot(data.plot, aes(factor(Tax), x=Group, y=mean, fill=Tax) ) +
    facet_wrap(~., scales="free_x") +
    geom_bar(width=0.75, stat='identity' ) +
    theme_bw() +
    ylim(c(0,1)) +
    theme(axis.line = element_line(colour = "black"),
          legend.text = element_text(face = "italic"),
          legend.justification = "top",
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank(),
          axis.text.x = element_text(angle = 45, hjust = 1, size=12),
          axis.text.y = element_text(size=12),
          strip.text.x = element_text(angle = 0, size = 10, color = "white"),
          axis.ticks.x=element_blank(),
          strip.background = element_rect(colour="white", fill="white", size=1.5, linetype="solid")) +
    ggtitle( paste("", "", sep=" ") ) +
    xlab("") +
    ylab("Proportion mapped contigs") +
    guides(fill = guide_legend(title = taxRank)) +
    scale_fill_manual(values=cols_top_genera[names(cols_top_genera) %in% data.plot$Tax])
fig4_6
ggsave(filename = paste0( outdir, "/Fig4_all.", taxRank,".severity.pdf"), width = 4.5, height = 9)

#
# Alpha
#

alpha <- ps %>%
    sample_data() %>%
    data.frame() %>%
    dplyr::select( MINIONID, Location, Case, Location_Severity, Severity, Shannon ) %>%
    dplyr::filter( Location == "Oral mucosa" | Location == "Stool" ) %>%
    dplyr::filter( Severity == "severe" | Severity == "interm" ) %>%
    dplyr::mutate( Location_Severity = as.character(Location_Severity) )

alpha$Location_Severity <- gsub(pattern = ":interm", replacement = "\nRDEB intermediate", x = alpha$Location_Severity)
alpha$Location_Severity <- gsub(pattern = ":severe",  replacement = "\nRDEB severe",   x = alpha$Location_Severity)
alpha$Location_Severity <- factor(alpha$Location_Severity,
                                  levels = c("Oral mucosa\nRDEB intermediate",
                                             "Oral mucosa\nRDEB severe",
                                             "Stool\nRDEB intermediate",
                                             "Stool\nRDEB severe") )

mycomp <- list(
    c("Oral mucosa\nRDEB intermediate", "Oral mucosa\nRDEB severe"),
    c("Stool\nRDEB intermediate", "Stool\nRDEB severe")
)

fig4_7 <- ggpubr::ggdotplot(data = alpha , x = "Location_Severity", y = "Shannon",
                  fill = "Severity", palette = cols_serverity,
                  size = 0.75,
                  xlab = "", ylab = "Shannon index of species diversity") +
    theme(legend.position = "none",
          axis.text.x=element_text(size=12, angle=45, hjust = 1),
          axis.text.y=element_text(size=12, angle=0),
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          axis.line = element_line(colour = "black")) +
    ylim(0,7) +
    #ggpubr::stat_compare_means(method = "kruskal.test", label.y = 8) +
    stat_compare_means(comparisons = mycomp, label.y = c(6, 6))
fig4_7
ggsave(filename = paste0( outdir, "/Fig4_all.Alpha.pdf"), width = 4.5, height = 9)


alpha <- ps %>%
    sample_data() %>%
    data.frame() %>%
    dplyr::select( MINIONID, Location, Case, Location_Case, Severity, Shannon ) %>%
    dplyr::filter( Location == "Oral mucosa" | Location == "Stool" ) %>%
    dplyr::mutate( Location_Case = as.character(Location_Case) )

alpha$Location_Case <- gsub(pattern = ":FALSE", replacement = " CTRL", x = alpha$Location_Case)
alpha$Location_Case <- gsub(pattern = ":TRUE",  replacement = " RDEB",   x = alpha$Location_Case)
alpha$Location_Case <- factor(alpha$Location_Case,
                                  levels = c("Oral mucosa CTRL",
                                             "Oral mucosa RDEB",
                                             "Stool CTRL",
                                             "Stool RDEB") )

mycomp <- list(
    c("Oral mucosa CTRL", "Oral mucosa RDEB"),
    c("Stool CTRL", "Stool RDEB")
)

fig4_7a <- ggpubr::ggdotplot(data = alpha , x = "Location_Case", y = "Shannon",
                            fill = "Case", palette = cols_alpha,
                            size = 0.75,
                            xlab = "", ylab = "Shannon index of species diversity") +
    theme(legend.position = "none",
          axis.text.x=element_text(size=12, angle=45, hjust = 1),
          axis.text.y=element_text(size=12, angle=0),
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          axis.line = element_line(colour = "black")) +
    ylim(0,7) +
    #ggpubr::stat_compare_means(method = "kruskal.test", label.y = 8) +
    stat_compare_means(comparisons = mycomp, label.y = c(6, 6))
fig4_7a
ggsave(filename = paste0( outdir, "/Fig4_all2.Alpha.pdf"), width = 4.5, height = 9)

# just to test this as well
alpha <- ps %>%
    sample_data() %>%
    data.frame() %>%
    dplyr::select( MINIONID, Location, Case, Location_Case, Severity, Shannon, AgeGroup ) %>%
    dplyr::filter( Location == "Oral mucosa" | Location == "Stool" ) %>%
    dplyr::mutate( Location_Case = as.character(Location_Case) ) %>%
    dplyr::mutate( Location_Case = paste(Location_Case, AgeGroup, sep = ":") )

alpha$Location_Case <- gsub(pattern = ":FALSE", replacement = " CTRL", x = alpha$Location_Case)
alpha$Location_Case <- gsub(pattern = ":TRUE",  replacement = " RDEB",   x = alpha$Location_Case)
alpha$Location_Case <- gsub(pattern = ":01",  replacement = " 0-4 years",   x = alpha$Location_Case)
alpha$Location_Case <- gsub(pattern = ":02",  replacement = " 5-10 years",   x = alpha$Location_Case)
alpha$Location_Case <- gsub(pattern = ":03",  replacement = " 11-22 years",   x = alpha$Location_Case)
alpha$Location_Case <- factor(alpha$Location_Case,
                              levels = c("Oral mucosa CTRL 0-4 years",
                                         "Oral mucosa RDEB 0-4 years",
                                         "Stool CTRL 0-4 years",
                                         "Stool RDEB 0-4 years",
                                         "Oral mucosa CTRL 5-10 years",
                                         "Oral mucosa RDEB 5-10 years",
                                         "Stool CTRL 5-10 years",
                                         "Stool RDEB 5-10 years",
                                         "Oral mucosa CTRL 11-22 years",
                                         "Oral mucosa RDEB 11-22 years",
                                         "Stool CTRL 11-22 years",
                                         "Stool RDEB 11-22 years")
)

mycomp <- list(
    c("Oral mucosa CTRL 0-4 years", "Oral mucosa CTRL 5-10 years"),
    c("Oral mucosa CTRL 5-10 years", "Oral mucosa CTRL 11-22 years"),
    c("Oral mucosa CTRL 0-4 years", "Oral mucosa CTRL 11-22 years"),

    c("Oral mucosa RDEB 0-4 years", "Oral mucosa RDEB 5-10 years"),
    c("Oral mucosa RDEB 5-10 years", "Oral mucosa RDEB 11-22 years"),
    c("Oral mucosa RDEB 0-4 years", "Oral mucosa RDEB 11-22 years"),

    c("Stool CTRL 0-4 years", "Stool CTRL 5-10 years"),
    c("Stool CTRL 5-10 years", "Stool CTRL 11-22 years"),
    c("Stool CTRL 0-4 years", "Stool CTRL 11-22 years"),

    c("Stool RDEB 0-4 years", "Stool RDEB 5-10 years"),
    c("Stool RDEB 5-10 years", "Stool RDEB 11-22 years"),
    c("Stool RDEB 0-4 years", "Stool RDEB 11-22 years")


)

fig4_7_supp <- ggpubr::ggdotplot(data = alpha , x = "Location_Case", y = "Shannon",
                                 fill = "Case", palette = cols_alpha,
                                 size = 0.75,
                                 xlab = "", ylab = "Shannon index of species diversity") +
    theme(legend.position = "none",
          axis.text.x=element_text(size=12, angle=45, hjust = 1),
          axis.text.y=element_text(size=12, angle=0),
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          axis.line = element_line(colour = "black")) +
    ylim(0,15) +
    #ggpubr::stat_compare_means(method = "kruskal.test", label.y = 8) +
    stat_compare_means(comparisons = mycomp, label.y = c(6, 6, 7, 8, 8, 9, 10, 10, 11, 12, 12, 13))
fig4_7_supp
ggsave(filename = paste0( outdir, "/Fig4_all3.Alpha.pdf"), width = 15, height = 9, plot = fig4_7_supp)


#
# Correlations Age Oral Mucosa
#

ps_ra <- transform_sample_counts(ps, function(x){x / sum(x)})
ps_ra <- subset_samples(physeq = ps_ra, Location == "Oral mucosa")

# Phyla
corr_df <- phylosmith::variable_correlation( ps_ra, treatment = c('Location', 'Case'),
                                             subset = c('Oral mucosa'),
                                             variables = 'Age',
                                             classification = "Phylum", method = 'spearman')

corr_df$Treatment <- gsub(pattern = " TRUE", replacement = " RDEB", x = corr_df$Treatment)
corr_df$Treatment <- gsub(pattern = " FALSE", replacement = " CTRL", x = corr_df$Treatment)

corr_df$p_clean <- corr_df$p
corr_df$p_clean <- format( round(corr_df$p_clean, 4), scientific = F ) # avoid exponentials
corr_df$p_clean[corr_df$p_clean >= 0.05] <- ""

fig4_8 <- corr_df %>%
    dplyr::mutate(X = gsub(pattern = "_", replacement = " ", X)) %>%
    ggplot(data = ., aes(x=Treatment, y=factor(X, levels = rev(levels(factor(X)))), fill=rho)) +
    geom_tile() +
    geom_text(aes(Treatment, X,
                  label = p_clean),
              color = "black", size = 4) +
    scale_fill_gradient2(low = "#2C7BB6", high = "#D7191C", mid = "white",
                         midpoint = 0, limit = c(-1,1), space = "Lab",
                         name="Spearman\nCorrelation") +
    theme(axis.line = element_line(colour = "black"),
          legend.position = "right",
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank(),
          axis.text.x = element_text(angle = 45, hjust = 1, size=12),
          strip.text.x = element_text(angle = 0, size = 12),
          axis.text.y = element_text(size=12, face = "italic"),
          axis.ticks.x=element_blank(),
          strip.background = element_rect(colour="white", fill="white", size=1.5, linetype="solid")) +
    xlab("") + ylab("")
fig4_8
ggsave(filename = paste0( outdir, "/Fig4_all.Correlation_Oral_Mucosa_Phylum.pdf"), width = 4.5, height = 9)

# Genera
( ps.sel <- microbiome::aggregate_taxa(x = ps_ra, level = "Genus") %>% phylosmith::taxa_filter(., treatment = c('Location'), frequency = 0.2) )
corr_df <- phylosmith::variable_correlation( ps.sel, treatment = c('Location', 'Case'),
                                             subset = c('Oral mucosa'),
                                             variables = 'Age',
                                             classification = "Genus", method = 'spearman')

corr_df$Treatment <- gsub(pattern = " TRUE", replacement = " RDEB", x = corr_df$Treatment)
corr_df$Treatment <- gsub(pattern = " FALSE", replacement = " CTRL", x = corr_df$Treatment)

corr_df$p_clean <- corr_df$p
corr_df$p_clean <- format( round(corr_df$p_clean, 4), scientific = F ) # avoid exponentials
corr_df$p_clean[corr_df$p_clean >= 0.05] <- ""

# keep only genera with p < 0.05
corr_genera <- corr_df$X[corr_df$p < 0.05]

fig4_9 <- corr_df %>%
    dplyr::filter(X %in% corr_genera) %>%
    dplyr::mutate(X = gsub(pattern = "_", replacement = " ", X)) %>%
    ggplot(data = ., aes(x=Treatment, y=factor(X, levels = rev(levels(factor(X)))), fill=rho)) +
    geom_tile() +
    geom_text(aes(Treatment, X,
                  label = p_clean),
              color = "black", size = 4) +
    scale_fill_gradient2(low = "#2C7BB6", high = "#D7191C", mid = "white",
                         midpoint = 0, limit = c(-1,1), space = "Lab",
                         name="Spearman\nCorrelation") +
    theme(axis.line = element_line(colour = "black"),
          legend.position = "right",
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank(),
          axis.text.x = element_text(angle = 45, hjust = 1, size=12),
          strip.text.x = element_text(angle = 0, size = 12),
          axis.text.y = element_text(size=12, face = "italic"),
          axis.ticks.x=element_blank(),
          strip.background = element_rect(colour="white", fill="white", size=1.5, linetype="solid")) +
    xlab("") + ylab("")
fig4_9
ggsave(filename = paste0( outdir, "/Fig4_all.Correlation_Oral_Mucosa_Genus.pdf"), width = 4.5, height = 9)

#
# Correlations Age Oral Mucosa
#

ps_ra <- transform_sample_counts(ps, function(x){x / sum(x)})
ps_ra <- subset_samples(physeq = ps_ra, Location == "Stool")

# Phyla
corr_df <- phylosmith::variable_correlation( ps_ra, treatment = c('Location', 'Case'),
                                             subset = c('Stool'),
                                             variables = 'Age',
                                             classification = "Phylum", method = 'spearman')

corr_df$Treatment <- gsub(pattern = " TRUE", replacement = " RDEB", x = corr_df$Treatment)
corr_df$Treatment <- gsub(pattern = " FALSE", replacement = " CTRL", x = corr_df$Treatment)

corr_df$p_clean <- corr_df$p
corr_df$p_clean <- format( round(corr_df$p_clean, 4), scientific = F ) # avoid exponentials
corr_df$p_clean[corr_df$p_clean >= 0.05] <- ""

fig4_10 <- corr_df %>%
    dplyr::mutate(X = gsub(pattern = "_", replacement = " ", X)) %>%
    ggplot(data = ., aes(x=Treatment, y=factor(X, levels = rev(levels(factor(X)))), fill=rho)) +
    geom_tile() +
    geom_text(aes(Treatment, X,
                  label = p_clean),
              color = "black", size = 4) +
    scale_fill_gradient2(low = "#2C7BB6", high = "#D7191C", mid = "white",
                         midpoint = 0, limit = c(-1,1), space = "Lab",
                         name="Spearman\nCorrelation") +
    theme(axis.line = element_line(colour = "black"),
          legend.position = "right",
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank(),
          axis.text.x = element_text(angle = 45, hjust = 1, size=12),
          strip.text.x = element_text(angle = 0, size = 12),
          axis.text.y = element_text(size=12, face = "italic"),
          axis.ticks.x=element_blank(),
          strip.background = element_rect(colour="white", fill="white", size=1.5, linetype="solid")) +
    xlab("") + ylab("")
fig4_10
ggsave(filename = paste0( outdir, "/Fig4_all.Correlation_Stool_Phylum.pdf"), width = 4.5, height = 9)

# Genera
( ps.sel <- microbiome::aggregate_taxa(x = ps_ra, level = "Genus") %>% phylosmith::taxa_filter(., treatment = c('Location'), frequency = 0.2) )
corr_df <- phylosmith::variable_correlation( ps.sel, treatment = c('Location', 'Case'),
                                             subset = c('Stool'),
                                             variables = 'Age',
                                             classification = "Genus", method = 'spearman')

corr_df$Treatment <- gsub(pattern = " TRUE", replacement = " RDEB", x = corr_df$Treatment)
corr_df$Treatment <- gsub(pattern = " FALSE", replacement = " CTRL", x = corr_df$Treatment)

corr_df$p_clean <- corr_df$p
corr_df$p_clean <- format( round(corr_df$p_clean, 4), scientific = F ) # avoid exponentials
corr_df$p_clean[corr_df$p_clean >= 0.05] <- ""

# keep only genera with p < 0.05
corr_genera <- corr_df$X[corr_df$p < 0.05]

fig4_11 <- corr_df %>%
    dplyr::filter(X %in% corr_genera) %>%
    dplyr::mutate(X = gsub(pattern = "_", replacement = " ", X)) %>%
    ggplot(data = ., aes(x=Treatment, y=factor(X, levels = rev(levels(factor(X)))), fill=rho)) +
    geom_tile() +
    geom_text(aes(Treatment, X,
                  label = p_clean),
              color = "black", size = 4) +
    scale_fill_gradient2(low = "#2C7BB6", high = "#D7191C", mid = "white",
                         midpoint = 0, limit = c(-1,1), space = "Lab",
                         name="Spearman\nCorrelation") +
    theme(axis.line = element_line(colour = "black"),
          legend.position = "right",
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank(),
          axis.text.x = element_text(angle = 45, hjust = 1, size=12),
          strip.text.x = element_text(angle = 0, size = 12),
          axis.text.y = element_text(size=12, face = "italic"),
          axis.ticks.x=element_blank(),
          strip.background = element_rect(colour="white", fill="white", size=1.5, linetype="solid")) +
    xlab("") + ylab("")
fig4_11
ggsave(filename = paste0( outdir, "/Fig4_all.Correlation_Stool_Genus.pdf"), width = 4.5, height = 9)

#
# Differential abundance Stool healthy vs RDEB
#

siglevel <- 0.05
location <- "Stool"

ps_sub <- subset_samples(ps, Location == location)

TAXLEVEL <- "Genus"
da_corncob <-
    corncob::differentialTest(formula = ~ Case + Age,
                              phi.formula = ~ Case + Age, # model to be fitted to the dispersion
                              formula_null = ~ Age, # Formula for mean under null, without response
                              phi.formula_null = ~ Case + Age, # Formula for overdispersion under null, without response
                              test = "Wald", boot = FALSE, B = 0,
                              data = ps_sub %>%  tax_glom(TAXLEVEL),
                              fdr_cutoff = siglevel)
da_corncob$significant_taxa
corncob::otu_to_taxonomy(OTU = da_corncob$significant_taxa, data = ps, level = TAXLEVEL)
da_plot <- plot(da_corncob, level = c("Phylum", "Genus")) +
    #ggplot2::xlim(c(-10,10)) +
    geom_point(size = 3)
da_plot
da_data <- data.frame(
    Taxa = corncob::otu_to_taxonomy(OTU = da_corncob$significant_taxa, data = ps_sub, level = TAXLEVEL),
    Taxa_veri = da_plot$data$taxa,
    p_fdr = da_corncob$p_fdr[da_corncob$p_fdr < siglevel & !is.na(da_corncob$p_fdr) ],
    Effect = da_plot$data$x,
    Error_min = da_plot$data$xmin,
    Error_max = da_plot$data$xmax
)

fig4_12 <- da_data %>%
    dplyr::mutate(Taxa = gsub(pattern = "_", replacement = " ", Taxa)) %>%
    ggplot(data = ., aes(x = Effect, y = factor(Taxa, levels = rev(levels(factor(Taxa)))))) +
    geom_point(size = 3) +
    geom_errorbar(aes(xmin=Error_min, xmax=Error_max), width=.3,
                  position=position_dodge(.9)) +
    geom_vline(xintercept=0, linetype="dashed",
               color = "grey45", size = 1) +
    theme(axis.line = element_line(colour = "black"),
          legend.position = "none",
          panel.grid.major = element_line(size = 0.25, linetype = 'solid',
                                          colour = "grey85"),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank(),
          axis.text.x = element_text(size=12),
          strip.text.x = element_text(angle = 0, size = 12),
          axis.text.y = element_text(size=12, face = "italic"),
          #axis.ticks.x = element_blank(),
          text = element_text(size = 12),
          strip.background = element_rect(colour="white", fill="white", size=1.5, linetype="solid")
    ) +
    xlim(-11, 11) +
    xlab("Effect size") + ylab("")
fig4_12
ggsave(filename = paste0( outdir, "/Fig4_all.Diff_Abundance_Stool_Genus.pdf"), width = 4.5, height = 4.5)

# Save workspace ----------------------------------------------------------

save.image(file='data/99_plots_paper.RData')

#
# Some pretty panel plotting
# Fig1
# Fig2
# Fig3
# Fig4
#

# Arrange plots -----------------------------------------------------------

load(file = "data/99_plots_paper.RData")

# Fig 1 (1..8) ------------------------------------------------------------

lay <- rbind(c(1,1,2,2,3,3,3,4,4,4),
             c(1,1,2,2,3,3,3,4,4,4),
             c(1,1,2,2,3,3,3,4,4,4),
             c(5,5,5,5,6,6,6,6,6,6),
             c(5,5,5,5,6,6,6,6,6,6),
             c(7,7,7,8,8,8,NA,NA,NA,NA),
             c(7,7,7,8,8,8,NA,NA,NA,NA))

# use arrangeGrob instead of grid.arrange to not draw the figure
fig1 <- gridExtra::arrangeGrob(fig1_1, fig1_2, fig1_3, fig1_4,
                                fig1_5, fig1_6,
                                fig1_7, fig1_8,
                                layout_matrix = lay) %>%
    as_ggplot(.) +
    cowplot::draw_plot_label(label = c("a", "b", "c", "d", "e", "f", "g", "h"), size = 15,
                             x = c(0, 0.2, 0.4, 0.7, 0, 0.4, 0, 0.3),
                             y = c(1, 1, 1, 1, 0.58, 0.58, 0.3, 0.3)
                    )
ggsave(filename = "plots_paper/panelplot_Fig1.pdf", width = 50, height = 40, units = 'cm', plot = fig1)

# Fig 2 (1..7) ------------------------------------------------------------

lay <- rbind(c(1,1,1,2,2,2),
             c(1,1,1,2,2,2),
             c(1,1,1,2,2,2),
             c(3,3,4,4,5,5),
             c(3,3,4,4,5,5),
             c(3,3,6,6,7,7),
             c(3,3,6,6,7,7))

fig2 <- gridExtra::arrangeGrob(fig2_1, fig2_2, fig2_3, fig2_4, fig2_6, fig2_5, fig2_7,
                                layout_matrix = lay) %>%
    as_ggplot(.) +
    cowplot::draw_plot_label(label = c("a", "b", "c", "d", "e", "f", "g"), size = 15,
                             x = c(0, 0.5, 0, 0.35, 0.68, 0.35, 0.68),
                             y = c(1, 1, 0.58, 0.58, 0.58, 0.3, 0.3))
ggsave(filename = "plots_paper/panelplot_Fig2.pdf", width = 40, height = 40, units = 'cm', plot = fig2)

# Fig 3 (1..7) ------------------------------------------------------------

lay <- rbind(c(1,2), c(3,4))

fig3 <- gridExtra::arrangeGrob(fig3_2 + theme(legend.position = "none"),
                                fig3_4 + theme(legend.position = "none"),
                                fig3_6 + theme(legend.position = "none"),
                                fig3_7,
                                layout_matrix = lay) %>%
    as_ggplot(.) +
    cowplot::draw_plot_label(label = c("a", "b", "c", "d"),
                             x = c(0, 0.5, 0, 0.5),
                             y = c(1, 1, 0.5, 0.5))
ggsave(filename = "plots_paper/panelplot_Fig3.pdf", width = 50, height = 30, units = 'cm', plot = fig3)

lay <- rbind(c(1), c(2), c(3))
fig3_supp <- gridExtra::arrangeGrob(fig3_1 + theme(legend.position = "none"),
                                fig3_3 + theme(legend.position = "none"),
                                fig3_5 + theme(legend.position = "none"),
                                layout_matrix = lay) %>%
    as_ggplot(.) +
    cowplot::draw_plot_label(label = c("a", "b", "c"),
                             x = c(0, 0, 0),
                             y = c(1, 0.68, 0.34))
ggsave(filename = "plots_paper/panelplot_Fig3_suppl.pdf", width = 30, height = 60, units = 'cm', plot = fig3_supp)

# Fig 4 (1..12) ------------------------------------------------------------

lay <- rbind(c(1,1,2,2,3,3),
             c(1,1,2,2,3,3),
             c(1,1,2,2,3,3),
             c(4,4,5,5,6,6),
             c(4,4,5,5,6,6),
             c(7,7,8,8,9,9),
             c(7,7,8,8,9,9))

fig4 <- gridExtra::arrangeGrob(fig4_2, fig4_6, fig4_4,
                                fig4_7a, fig4_8, fig4_9,
                                fig4_12, fig4_10, fig4_11,
                                layout_matrix = lay) %>%
    as_ggplot(.) +
    cowplot::draw_plot_label(label = c("a", "b", "c", "d", "e", "f", "g", "h", "i"),
                                            x = c(0, 0.33, 0.66, 0, 0.33, 0.66, 0, 0.33, 0.66),
                                            y = c(1, 1, 1, 0.58, 0.58, 0.58, 0.3, 0.3, 0.3))
ggsave(filename = "plots_paper/panelplot_Fig4.pdf", width = 40, height = 40, units = 'cm', plot = fig4)

lay <- rbind(c(1,2,3))
fig4_supp <- gridExtra::arrangeGrob(fig4_1, fig4_5, fig4_3,
                                layout_matrix = lay) %>%
    as_ggplot(.) +
    cowplot::draw_plot_label(label = c("a", "b", "c"),
                             x = c(0, 0.33, 0.66),
                             y = c(1,1,1))
ggsave(filename = "plots_paper/panelplot_Fig4_supp.pdf", width = 40, height = 20, units = 'cm', plot = fig4_supp)
