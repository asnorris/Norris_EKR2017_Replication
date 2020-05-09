# stylized_seg.R: script to produce Figure 1 in the paper.

# Load packages
library(gridExtra)
library(grid)
library(RColorBrewer)
library(rgdal)
library(rgeos)
library(maptools)
library(seg)
library(splancs)
library(tidyverse)

## This creates the dimensions of the hypothetical electoral districts
# Set dimensions of grid, cell size, and number of localities
grid.dim <- c(4, 4)
cell.size <- 1
no.loc <- grid.dim[1] * grid.dim[2]

# Generate grid
grd <- GridTopology(c(0, 0), c(cell.size, cell.size), grid.dim)
grd.sp <- as.SpatialPolygons.GridTopology(grd)
polyg <- fortify(grd.sp) #for ggplot

# Now we’ll generate a population for each of two electoral districts,
# assuming we know the demographic makeup of localities within each
# district.

# Set number of members of each group
g1 <- 1000
g2 <- 1000

# Number of individuals in each locality
nloc <- g1 / (no.loc/2)

# Data for segregated electoral district
dta.seg <- data.frame(
  id = 1:no.loc,
  group1 = c(rep(nloc, 3), 0, rep(nloc, 3), 0, rep(nloc, 2), rep(0, 2), rep(0, 4)),
  j = "(a)"
)

dta.seg$group2 <- ifelse(dta.seg$group1 == 0, nloc, 0)

# Data for integrated electoral district
dta.cb <- data.frame(
  id = 1:no.loc,
  group1 = rep(c(g1 / (no.loc/2), 0), no.loc/2),
  group2 = rep(c(0, g2 / (no.loc/2)), no.loc/2),
  j = "(b)"
)

dta.cb$group1[c(5:8, 13:16)] <- dta.cb$group1[c(2:5, 2:5)]
dta.cb$group2[c(5:8, 13:16)] <- dta.cb$group2[c(2:5, 2:5)]

# Plot ----------------------

# Function for graphing circle
circleFn <- function(center = c(1,2), diameter = 1, npoints = 100){
  r = diameter / 2
  tt <- seq(0,2*pi,length.out = npoints)
  xx <- center[1] + r * cos(tt)
  yy <- center[2] + r * sin(tt)
  return(data.frame(x = xx, y = yy))
}

# Graphing theme
theme <- theme_bw() + theme(
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  panel.border = element_blank(),
  axis.ticks = element_blank(),
  axis.text = element_blank(),
  strip.background = element_blank(),
  strip.text.x = element_text(size = 10),
  plot.title = element_text(size = 11)
)

# Plot
dta.list <- list(dta.seg, dta.cb)

plots <- lapply(dta.list, function(dta) {

  g1 <- dotsInPolys(grd.sp, as.integer(dta$group1))
  g2 <- dotsInPolys(grd.sp, as.integer(dta$group2))

  points <- rbind(data.frame(coordinates(g1), group = "group 1"),
                  data.frame(coordinates(g2), group = "group 2"))

  pg <- data.frame(x = 1, y = 2)
  catchment <- circleFn(diameter = 2.3)

  ggplot(data = polyg, aes(x = long, y = lat)) +
    geom_polygon(aes(group = id), color = "black", fill = "white", size = 1.1) +
    geom_point(data = points, aes(x = x, y = y, color = group, shape = group), size = 1.4) +
    geom_point(data = pg, aes(x, y), size = 5, shape = 18) +
    geom_polygon(data = catchment, aes(x, y), alpha = 0.15, fill = "red") +
    scale_color_manual(values = c("gray75", "black"), guide = F) +
    scale_shape_manual(values = c(0, 16), guide = F) +
    coord_equal() +
    labs(x = NULL, y = NULL) +
    theme
})

pdf("Replication/figures/stylized_seg.pdf", width = 5.7, height = 3.2)
do.call(grid.arrange, c(plots, nrow = 1))
dev.off()
