# 2sp asymmetric competition - calculation of critical curvature

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

# visualise
pal <- choose_palette()

curve.data %>%
  filter(gap >= 30) %>%
  ggplot(aes(x = delta, y = K, color = gap, group = gap)) +
  geom_point(size = 3) + 
  scale_color_gradientn(colors = pal(60)) + 
  labs(x = expression(paste("competition asymmetry (", delta, ")")), 
       y = "critical curvature (K)", 
       color = paste("edge-edge","pillar spacing", sep = "\n")) +
  theme(legend.position = c(0.8, 0.8))
