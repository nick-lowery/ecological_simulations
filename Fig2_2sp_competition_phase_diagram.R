#
###############################################################################
### Figure 2 - Phase diagram of competition between two species ###############
###############################################################################

# Nick Lowery and Tristan Ursell
# 2018
# 
# This script generates Figure 2A from our preprint:
# Structured environments fundamentally alter dynamics and stability of ecological communities
# https://www.biorxiv.org/content/early/2018/07/10/366559
# 
# The script reads in correlation and metadata files generated from the LV_competition_two_species_asym.m
# (or LV_competition_two_species_asym_regLV.m, for the supplemental figure showing results of the 
# classic competitive LV model) script(s). The functions to generate the data from the
# raw .mat data files are included but commented out.  Note the code below is for the 'regLV' 
# version of the model, but will generate Figure 2A when fed in the corresponding data.

###############################################################################

#######################################################
### read in and clean raw data ########################
#######################################################

library(tidyverse)
library(cowplot)
library(colorspace)
# library(R.matlab)
# 
# # Read in data
# asym.lattice.raw <- list.files("2sp_asym_regLV", pattern = ".mat", full.names = T) %>%
#   lapply(., readMat)
# 
# # organize into a data frame
# asym.lattice <- asym.lattice.raw %>% {
#   tibble(
#     L = map_dbl(., "L"),
#     D = map_dbl(., "D"),
#     N = map_dbl(., "N"),
#     dt = map_dbl(., "dt"),
#     alpha = map_dbl(., "alpha"),
#     pillars = map_dbl(., "pillarq"),
#     R = map_dbl(., "R"), 
#     dx = map_dbl(., "dx"),
#     rep = map_dbl(., "rep"),
#     t.stop = map_dbl(., "t.stop"), 
#     pillar.weight = map_dbl(., "filt.all.weight"),
#     data = map(., `[`, c("meanA.out", "meanB.out")) %>% map(data.frame)
#   )
#   } %>%
#   mutate(data = map(data, ~mutate(.x, t.step = seq(0.4, by = 0.4, length.out = nrow(.x)))),
#          t.stop = t.stop * dt) %>%
#   unnest() %>%
#   rename(A = meanA.out, B = meanB.out) %>%
#   mutate(A.extinct = ifelse(A <= ((2*R)^2 - pi*R^2)/4/pillar.weight, 1, 0),  
#          B.extinct = ifelse(B <= ((2*R)^2 - pi*R^2)/4/pillar.weight, 1, 0))

# write_csv(asym.lattice, "2sp_asym_regLV_cleaned_data.csv")
asym.lattice <- read_csv("E:/ecological_modelling/2sp_asym_regLV_cleaned_data.csv")

#######################################################
### calculate extinction stats and generate figure ####
#######################################################

# summarize extinction occurences
asym.lattice.summary <- asym.lattice %>%
  group_by(D, alpha, pillars, R, dx, rep, t.stop) %>%
  summarise(A.extinction = max(A.extinct),
            B.extinction = max(B.extinct),
            t.A.extinct = ifelse(max(A.extinct), t.step[which.max(A.extinct)], Inf),
            t.B.extinct = ifelse(max(B.extinct), t.step[which.max(B.extinct)], Inf))

# create color palette (function generated from colorspace GUI)
pal2 <- function (n, h = c(300, 75), c. = c(35, 95), l = c(15, 90), 
                 power = c(0.8, 1.2), fixup = TRUE, gamma = NULL, alpha = 1, 
                 ...) 
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

# calculate extinction frequencies, and generate phase diagram
asym.lattice.phase <- asym.lattice.summary %>%
  group_by(D, alpha, R) %>%
  summarise(prop.extinct = (sum(A.extinction)+sum(B.extinction))/n()) %>%
  mutate(R.nl = round(R/sqrt(D)/1.29),
         R.nl = fct_relevel(factor(R.nl), "0", after = Inf)) %>%
  ggplot(., aes(x = alpha, y = R.nl, fill = prop.extinct)) +
  geom_tile() +
  geom_hline(yintercept = 9.5, color = "grey50", size = 2) +
  scale_fill_gradientn(colors = pal2(7)) +
  labs(x = expression(paste("competition coefficient (", alpha, ")")),
       y = expression(paste("pillar radius (", Delta, "x - 2R/ 1.29", lambda, ")")),
       fill = paste("extinction","frequency","", sep = "\n")) +
  theme(panel.background = element_rect(fill = "grey80"))

# display the figure
asym.lattice.phase

