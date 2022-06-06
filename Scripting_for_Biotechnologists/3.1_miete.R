#load tidyverse library
library(tidyverse)

#import data, skip 1st row
miete <- read_tsv("/export/home/up0289/courses/scripting_for_biotechnologists/report/R_hw/miete.tsv", skip = 1)



#apartment 1: less than 50 sqm, 1 room, built after 1966
apt_1 <- miete  %>% 
  filter(miete$wfl < 50, rooms == 1, miete$bj > 1966)

#apartment 2: 3 rooms, more than 80 sqm, built after 1966
apt_2 <- miete %>% 
  filter(miete$wfl > 80, rooms >= 3, miete$bj > 1966)

#my (sadly fictional) apartment: between 48 and 55 sqm, 2 rooms, built after 1980
apt_my <- miete %>% 
  filter(miete$wfl > 50, miete$wfl < 55, rooms == 2, miete$bj > 1980)



#Price range of these apartments
range(apt_1$nm)
range(apt_2$nm)
range(apt_my$nm)

#Price range per square meter of these apartments, rounded to 2 digits
round(digits = 2, range(apt_1$nmqm))
round(digits = 2, range(apt_2$nmqm))
round(digits = 2, range(apt_my$nmqm))
