library(mapdeck)
library(sf)
key <- 'pk.eyJ1IjoiY2hyaXN0aW5hYnVlbG93IiwiYSI6ImNsZ3k5aXVzaDA0azUzcm83czBvNG5scTcifQ.KWljR6tNZGlHhSPZlntzFw'

typ <- st_read('data/typologies/typologies_220km2_plus.gpkg') %>% st_simplify(1000, preserveTopology = T)
mapdeck(token = key, style = mapdeck_style("light")) %>% add_polygon(data = typ, fill_colour = 'area_ha', palette = 'inferno')
