###############################################################################
### Figure 2B - Critical curvature of competition interfaces ##################
###############################################################################

# Nick Lowery and Tristan Ursell
# 2018
# 
# This script generates Figure 2B from our preprint:
# Structured environments fundamentally alter dynamics and stability of ecological communities
# https://www.biorxiv.org/content/early/2018/07/10/366559
# 
# The script reads in metadata files generated from the LV_competition_two_species_asym_curvature.m
# function. The functions to generate the cleaned data from the raw .mat data 
# files are included but commented out.

###############################################################################

#library(R.matlab)
library(tidyverse)
library(cowplot)
library(colorspace)

# # Read in data
# curve.data.raw <- list.files("E:/ecological_modelling/2sp_asym_curve", pattern = ".mat", full.names = T) %>%
#   lapply(., readMat)
# 
# # organize into data frame
# curve.data <- curve.data.raw %>% {
#   tibble(
#     H = map_dbl(., "H"),
#     P = map_dbl(., "P"),
#     delta = map_dbl(., "delta"),
#     R = map_dbl(., "R"),
#     K = map_dbl(., "K"),
#     extinct = map_dbl(., "extinct")
#   )
# } %>%
#   mutate(dx = H - 2*R,
#          gap = H - 4*R)
# 
# write_csv(curve.data, "E:/ecological_modelling/2sp_critical_curvature.csv")
curve.data <- read_csv("E:/ecological_modelling/2sp_critical_curvature.csv")

# define color palette (generated from colorspace GUI)
pal <- function (n, h = c(300, 75), c. = c(35, 95), l = c(15, 90), power = c(0.8, 
                 1.2), fixup = TRUE, gamma = NULL, alpha = 1, ...) 
  {
    if (!is.null(gamma)) 
        warning("'gamma' is deprecated and has no effect")
    if (n < 1L) 
        return(character(0L))
    h <- rep(h, length.out = 2L)
    c <- rep(c., length.out = 2L)
    l <- rep(l, length.out = 2L)
    power <- rep(power, length.out = 2L)
    rval <- seq(1, 0, length = n)
    rval <- hex(polarLUV(L = l[2L] - diff(l) * rval^power[2L], 
        C = c[2L] - diff(c) * rval^power[1L], H = h[2L] - diff(h) * 
            rval), fixup = fixup, ...)
    if (!missing(alpha)) {
        alpha <- pmax(pmin(alpha, 1), 0)
        alpha <- format(as.hexmode(round(alpha * 255 + 1e-04)), 
            width = 2L, upper.case = TRUE)
        rval <- paste(rval, alpha, sep = "")
    }
    return(rval)
  }

# create plot
curve.data %>%
  # filter out data for which curve estimation is unreliable (due to small finite number of pixels)
  filter(gap >= 30) %>%
  ggplot(aes(x = delta, y = K, color = gap, group = gap)) +
  geom_point(size = 3) + 
  scale_color_gradientn(colors = pal(60)) + 
  labs(x = expression(paste("competition asymmetry (", delta, ")")), 
       y = "critical curvature (K)", 
       color = paste("edge-edge","pillar spacing", sep = "\n")) +
  theme(legend.position = c(0.8, 0.8))
