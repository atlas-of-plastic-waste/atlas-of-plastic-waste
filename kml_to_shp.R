library(sf)
library(tidyverse)

apw.shp <- st_read("atlas_of_plastic.kml")

apw.shp <- separate(apw.shp, col = "Description", into = c("Country", "site_type", "waste_type", "year", "observer"), sep = ", ")

apw.shp <- st_zm(apw.shp, drop = TRUE, what = "ZM")

st_write(apw.shp, "apw/atlas_of_plastic.shp", append=FALSE)
