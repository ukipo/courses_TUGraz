#load library
library(tidyr)
library(dplyr)
library(ggplot2)

#import data
mystery <- read_tsv("/export/home/up0289/courses/scripting_for_biotechnologists/report/R_hw/mystery.csv")

#summary statistics x
mystery %>%
  summarize(min_x = min(mystery$x), max_x = max(mystery$x), mean_x = mean(mystery$x), sd_x = sd(mystery$x))

#summary statistics y
mystery %>%
  summarize(min_y = min(mystery$y), max_y = max(mystery$y), mean_y = mean(mystery$y), sd_y = sd(mystery$y))

#correlation between variables
cor(mystery$x, mystery$y)

#reorder datasets
mystery_new <- mystery
mystery_new$dataset <- factor(mystery_new$dataset,
                              levels = c("d1", "d2", "d3", "d4", "d5", "d6", "d7", "d8", "d9", "d10", "d11", "d12", "d13"))

#variables
x <- mystery_new$x
y <- mystery_new$y
dataset <- mystery_new$dataset

#group by dataset and summarize statistics again
mystery_new %>%
  group_by(dataset) %>%
  summarize(min_x = min(x), max_x = max(x), mean_x = mean(x), sd_x = sd(x))


mystery_new %>%
  group_by(dataset) %>%
  summarize(min_y = min(y), max_y = max(y), mean_y = mean(y), sd_y = sd(y))


#plot
ggplot(data = mystery_new, mapping = aes(x = mystery_new$x, y = mystery_new$y, color = mystery_new$dataset)) +
  geom_point()+
  facet_wrap(facets = vars(mystery_new$dataset)) +
  labs(title = "Plots of the various datasets",
       x = "X coordinate",
       y = "Y coordinate",
       color = "Dataset") +
  theme_minimal()

#What is the problem when only looking at summary statistics?
