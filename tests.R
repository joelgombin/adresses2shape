# test API
library(dplyr)
library(multidplyr)
library(rgeos)
data(paris2012)

paris <- rgdal::readOGR("/media/Data/Dropbox/Cartographie/arrondissements paris", "ARRONDISSEMENT")

# exécution en parallèle

cl <- get_default_cluster()
cl %>%
  cluster_library("dplyr") %>%
  cluster_library("adresses2shape") %>%
  cluster_library("maptools") %>%
  cluster_library("rgeos") %>%
  cluster_assign_value("paris", value = paris)

test <- paris2012 %>%
      mutate(adresse = paste(numero, voie, nom)) %>%
      mutate(insee = paste0("751", arrondissement)) %>%
      partition(insee, cluster = cl) %>%
      group_by(insee) %>%
      do(arrdts = {
        ID <- .[["insee"]]
        geocode(., "adresse", "insee") %>%
           create_shapes("ID", threshold = 0.85) %>%
           crop_voronois(paris[as.character(paris@data$C_ARINSEE) %in% ID,])
      }
         )

test <- test %>% collect()

test2 <- do.call(rbind, test$arrdts)
