# load('~/Documents/Github/test.dbpmm/Real Data/[Used] Koerber_et_al/data/H043-W99H5K.RData')
#
 require(tidyverse)
#
# snvs = dataset$primary_SNV %>%
#   rename(
#     chr = CHROM,
#     from = POS,
#     ref = REF,
#     alt = ALT
#   ) %>%
#   mutate(to = as.numeric(from) + nchar(alt))
#
# cna = dataset$primary_CNA %>%
#   separate(col = genotype, into = c('Major', 'minor'))

# example_dataset_CNAqc = list(snvs = snvs, cna = cna, purity = 0.89)
# usethis::use_data(example_dataset_CNAqc)

load('/Users/gcaravagna/Documents/Github/CNAqc/data/example_dataset_CNAqc.rda')

x = CNAqc::init(example_dataset_CNAqc$snvs, example_dataset_CNAqc$cna,example_dataset_CNAqc$purity)
x = CNAqc::analyze_peaks(x)

CNAqc::plot_segments(x)

CNAqc::plot_karyotypes(x)

CNAqc::plot_peaks_analysis(x)

cowplot::plot_grid(
  CNAqc::plot_counts(x),
  CNAqc::plot_vaf(x, N = 2000),
  CNAqc::plot_depth(x, N = 2000),
  CNAqc::plot_segments(x),
  align = 'v', nrow = 4, rel_heights = c(.15, .15, .15, .8))

install.packages('BioCircos')
library(BioCircos)

BioCircos()

tracklist = BioCircosTextTrack('myTextTrack', 'Some text', size = "2em", opacity = 0.5,
                               x = -0.67, y = -0.5)

BioCircos(tracklist, genomeFillColor = "PuOr",
          chrPad = 0, displayGenomeBorder = FALSE,
          genomeTicksLen = 2, genomeTicksTextSize = 0, genomeTicksScale = 1e+8,
          genomeLabelTextSize = "9pt", genomeLabelDy = 0)


snvs

library(BioCircos)

# Chromosomes on which the points should be displayed
points_chromosomes = c(1:5)
# Coordinates on which the points should be displayed
points_coordinates = c(102621, 140253678, 98567307, 28937403, 20484611)
# Values associated with each point, used as radial coordinate
#   on a scale going to minRadius for the lowest value to maxRadius for the highest value
points_values = 0:4

tracklist = BioCircosSNPTrack('mySNPTrack', points_chromosomes, points_coordinates,
                              points_values, colors = c("tomato2", "darkblue"), minRadius = 0.5, maxRadius = 0.9)

# Background are always placed below other tracks
tracklist = tracklist + BioCircosBackgroundTrack(
  "myBackgroundTrack",
  minRadius = 0.5, maxRadius = 0.9,
  borderColors = "#AAAAAA", borderSize = 0.6, fillColors = "#B3E6FF")

BioCircos(tracklist, genomeFillColor = "PuOr",
          chrPad = 0.05, displayGenomeBorder = FALSE, yChr =  FALSE,
          genomeTicksDisplay = FALSE,  genomeLabelTextSize = 18, genomeLabelDy = 0)


snvs_1 = snvs %>% filter(chr == 'chr1')

# Chromosomes on which the points should be displayed
points_chromosomes = rep('1', nrow(snvs_1))
# Coordinates on which the points should be displayed
points_coordinates = as.numeric(snvs_1$from)
# Values associated with each point, used as radial coordinate
#   on a scale going to minRadius for the lowest value to maxRadius for the highest value
points_values = snvs_1$VAF

tracklist = BioCircosSNPTrack(
  'mySNPTrack',
  points_chromosomes, points_coordinates,
  points_values, colors = c("tomato2", "darkblue"), minRadius = 0.5, maxRadius = 0.9)

# Background are always placed below other tracks
tracklist = tracklist + BioCircosBackgroundTrack(
  "myBackgroundTrack",
  minRadius = 0.5, maxRadius = 0.9,
  borderColors = "#AAAAAA", borderSize = 0.6, fillColors = "#B3E6FF")

BioCircos(tracklist, genomeFillColor = "PuOr",
          chrPad = 0.05, displayGenomeBorder = FALSE, yChr =  FALSE,
          genomeTicksDisplay = FALSE,  genomeLabelTextSize = 18, genomeLabelDy = 0)

x = rnorm(160) * 1e4
factors = paste0('chr', sample(1:16, 160, replace = TRUE))
circos.initialize(factors = factors, x = x)
circos.trackHist(factors = factors, x = x, col = "#999999",
                 border = "#999999")
circos.trackHist(factors = factors, x = x, bin.size = 0.1,
                 col = "#999999", border = "#999999")
circos.trackHist(factors = factors, x = x, draw.density = TRUE,
                 col = "#999999", border = "#999999")


set.seed(999)
bed = generateRandomBed()
head(bed)

circos.initializeWithIdeogram()
circos.trackHist(factors = snvs$chr, x = as.numeric(snvs$from), col = "#999999",
                 bin.size = 1e6,
                 border = "#999999")

# circos.initializeWithIdeogram(chromosome.index = paste0("chr", c(3,5,17,8)))
circos.initializeWithIdeogram(chromosome.index = paste0("chr", 1))

snvs_1

circos.trackHist(factors = snvs_1$chr, x = as.numeric(snvs_1$from), col = "#999999",
                 bin.size = 1e6,
                 border = "#999999")


cna_1 = cna %>% filter(chr == 'chr1')

circos.genomicTrackPlotRegion(cna_1, ylim = c(0, 1),
                              panel.fun = function(region, value, ...) {
                                circos.genomicLines(region, value, ...)
                              })



