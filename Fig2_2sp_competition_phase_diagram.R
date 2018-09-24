# 2sp asymmetric competition, classic competition model (as opposed to just spacing)



#######################################################
### read in and clean raw data ########################
#######################################################

library(R.matlab)

# Read in data
asym.lattice.raw <- list.files("E:/ecological_modelling/2sp_asym_regLV2", pattern = ".mat", full.names = T) %>%
  lapply(., readMat)

# get things organized
asym.lattice <- asym.lattice.raw %>% {
  tibble(
    L = map_dbl(., "L"),
    D = map_dbl(., "D"),
    N = map_dbl(., "N"),
    P = map_dbl(., "P"),
    dt = map_dbl(., "dt"),
    delta = map_dbl(., "delta"),
    pillars = map_dbl(., "pillarq"),
    R = map_dbl(., "R"),
    dx = map_dbl(., "dx"),
    rep = map_dbl(., "rep"),
    t.stop = map_dbl(., "t.stop"),
    pillar.weight = map_dbl(., "filt.all.weight"),
    data = map(., `[`, c("meanA.out", "meanB.out")) %>% map(data.frame)
  )
} %>%
  mutate(data = map(data, ~mutate(.x, t.step = seq(nrow(.x)))),
         t.stop = t.stop * dt) %>%
  unnest() %>%
  rename(A = meanA.out, B = meanB.out) %>%
  mutate(A.extinct = ifelse(A <= ((2*R)^2 - pi*R^2)/4/pillar.weight, 1, 0),
         B.extinct = ifelse(B <= ((2*R)^2 - pi*R^2)/4/pillar.weight, 1, 0))

write_csv(asym.lattice, "E:/ecological_modelling/2sp_asym_regLV_cleaned_data2.csv")

library(tidyverse)
library(cowplot)
library(colorspace)

asym.lattice <- read_csv("E:/ecological_modelling/2sp_asym_regLV_cleaned_data2.csv")

asym.lattice.summary <- asym.lattice %>%
  group_by(D, delta, pillars, R, dx, rep, t.stop) %>%
  summarise(A.extinction = max(A.extinct),
            B.extinction = max(B.extinct),
            t.A.extinct = ifelse(max(A.extinct), t.step[which.max(A.extinct)], Inf),
            t.B.extinct = ifelse(max(B.extinct), t.step[which.max(B.extinct)], Inf))

# phase plot
pal2 <- choose_palette()

asym.lattice.phase <- asym.lattice.summary %>%
  group_by(D, delta, R) %>%
  summarise(prop.extinct = (sum(A.extinction)+sum(B.extinction))/n()) %>%
  mutate(R.nl = round(R/sqrt(D)/1.29),
         R.nl = fct_relevel(factor(R.nl), "0", after = Inf)) %>%
  ggplot(., aes(x = delta, y = R.nl, fill = prop.extinct)) +
  geom_tile() +
  geom_hline(yintercept = 9.5, color = "grey50", size = 2) +
  scale_fill_gradientn(colors = pal2(7)) +
  scale_x_continuous(breaks = seq(0, 0.02, by = 0.002), 
                     labels = c("1","","1.1","","1.2","","1.3","","1.4","","1.5")) +
  labs(x = expression(paste("competitive asymmetry (", P[A] / P[B], ")")),
       y = expression(paste("edge-edge pillar spacing (", Delta, "x - 2R/ 1.29", lambda, ")")),
       fill = paste("extinction","frequency","", sep = "\n")) +
  theme(panel.background = element_rect(fill = "grey80"))
asym.lattice.phase

