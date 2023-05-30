# plot histogram of typology size distribution (log-scale)

library(sf)
library(tidyverse)

typ <- st_read('data/typologies/Mangrove_Typology_v3_Composite_valid.gpkg')

sf_use_s2(TRUE)

typ_area <- typ %>% 
  mutate(area_ha = as.numeric(st_area(.))/10000) # calculate area with s2 on, so approximating a sphere

df <- st_drop_geometry(typ_area)
st_write(select(df, Type, area_ha), 'outputs/processed-data/typology-area.csv', row.names = F)

ggplot(df) +
  geom_histogram(aes(log(area_ha))) +
  theme_classic()

ggsave('outputs/typology-size-histogram.png')

# filter large sites

typ2 <- typ_area %>% 
  filter(area_ha > exp(10))

st_write(typ2, 'data/typologies/typologies_220km2_plus.gpkg')
