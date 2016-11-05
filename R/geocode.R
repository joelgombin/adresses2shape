#' géocode des adresses grâce à l'API BANO
#'
#' @param data un \code{dataframe} contenant, a minima, les adresses à géocoder
#' @param adresses nom de la colonne contenant les adresses à géocoder
#' @param codeINSEE nom d'une colonne contenant le code INSEE de la commune. Si \code{NULL} (par défaut), ce paramètre n'est pas utilisé.
#'
#' @return un dataframe avec des colonnes supplémentaires avec la géolocalisation
#' @import magrittr
#' @export
#'
geocode <- function(data, adresses, codeINSEE = NULL) {
  # faire une sorte de test sur les adresses et éventuellement sur les IDs

  # TODO : vérifier que le CSV est inférieur à 8 MO (limite de l'API)

  # sauver adresses dans un CSV pour les envoyer à l'API BANO
  tmp <- paste0(tempfile(), ".csv")
  vars <- adresses
  if (!is.null(codeINSEE))
    vars <- c(vars, codeINSEE)

  write.csv(data %>% dplyr::select_(.dots = vars), file = tmp, row.names = FALSE) # be nice with the server/bandwidth, just upload what's necessary!

  # envoyer à BANO
  baseURL <- "http://api-adresse.data.gouv.fr/search/csv/"
  message("Géocodage en cours...")
  queryResults <- httr::POST(baseURL,
                       body = list(data = httr::upload_file(tmp),
                       columns = adresses,
                       citycode = codeINSEE)
  )
  # TODO : gérer les retours en erreur de l'API
  if (stringr::str_detect(queryResults$headers$`content-type`, "text/csv")) {
    locations <- readr::read_csv(queryResults$content, col_types = paste0("cc", paste0(rep("?", 13), collapse = "")))
  } else {
    stop("Erreur sur le résultat, code de l'erreur : ", httr::status_code(queryResults))
  }

  return(dplyr::left_join(data, locations))

}



#' Create Voronoi polygons
#'
#' @param x un objet de type \code{Spatial*}, ou bien un dataframe dont les deux première colonnes comportent les coordonnées x et y des points à partir desquels générer les polygones de Voronoi.
#' @import sp
#'
voronoipolygons <- function(x) {
  # source : http://stackoverflow.com/questions/9403660/how-to-create-thiessen-polygons-from-points-using-r-packages
  # à paralléliser !
  if (.hasSlot(x, 'coords')) {
    crds <- x@coords
  } else crds <- x
  z <- deldir::deldir(crds[[1]], crds[[2]])
  w <- deldir::tile.list(z)
  #   polys <- vector(mode='list', length=length(w))
  polys <- plyr::llply(seq(along=w), .fun = function(i) {
    pcrds <- cbind(w[[i]]$x, w[[i]]$y)
    pcrds <- rbind(pcrds, pcrds[1,])
    Polygons(list(Polygon(pcrds)), ID=as.character(i))
    #  }, .progress = "text", .parallel = TRUE, .paropts=list(.packages="sp", .export=c("w")))
  }) # TODO : réintégrer parallalélisme + barre de progression
  #   for (i in seq(along=polys)) {
  #     pcrds <- cbind(w[[i]]$x, w[[i]]$y)
  #     pcrds <- rbind(pcrds, pcrds[1,])
  #     polys[[i]] <- Polygons(list(Polygon(pcrds)), ID=as.character(i))
  #   }
  SP <- SpatialPolygons(polys)
  voronoi <- SpatialPolygonsDataFrame(SP, data=data.frame(x=crds[,1],
                                                          y=crds[,2], row.names=sapply(slot(SP, 'polygons'),
                                                                                       function(x) slot(x, 'ID'))))
}

#' Crée des polygones de Voronoi à partir de coordonnées géographiques lon/lat
#'
#' @param df un dataframe avec les coordonnées géographiques des points ainsi qu'un identifiant par polygone à constituer
#' @param ID nom de la colonne contenant les identifiants
#' @param x nom de la colonne contenant les longitudes
#' @param y nom de la colonne contenant les latitudes
#' @param threshold seuil de confiance dans les données géolocalisées à respecter. Si \code{NA} (défaut), pas de seuil.
#' @param threshold_score nom de la colonne contenant le degré de confiance dans les données géolocalisées. Par défaut, \code{"result_score"}, nom donné par l'API adresses
#' @param clip_by_chull valeur logique. Si \code{TRUE}, les polygones de Voronoi sont réduits à leur enveloppe complexe (le pavage ne couvre donc pas l'intégralité de l'espace mais se limite aux espaces sur lesquels il y a des adresses). Cela est utile par exemple lorsque d'importants espaces non habités existent. \code{FALSE} par défaut.
#'
#' @return un objet \code{SpatialPolygonsDataFrame} contenant les polygones de Voronoi recherchés
#' @export
#' @import maptools rgeos
#'
create_shapes <- function(df, ID, x = "longitude", y = "latitude", threshold = NA, threshold_score = "result_score", clip_by_chull = FALSE) {
  if (!is.na(threshold)) {
    filter_criteria <- lazyeval::interp(~ threshold_score > threshold, threshold_score = as.name(threshold_score))
    df <- df %>% filter_(filter_criteria)
  }
#  df <- df %>% dplyr::group_by_(IDs)
  df <- df %>% dplyr::distinct_(x, y, .keep_all = TRUE)
  res <- voronoipolygons(df %>% ungroup %>% select_(x, y) %>% as.data.frame)
  res2 <- unionSpatialPolygons(res, df[[ID]])
  res2 <- sp::SpatialPolygonsDataFrame(res2, data.frame(ID = names(res2), row.names = names(res2), stringsAsFactors = FALSE))

  if (clip_by_chull) {
    ch <- chull(df %>% ungroup %>% select_(x, y) %>% as.data.frame)
    ch <- df %>% ungroup %>% select_(x, y) %>% slice(c(ch, ch[1])) %>% as.data.frame()
    ch <- SpatialPolygons(list(Polygons(list(Polygon(ch)), ID=1)))
    proj4string(ch) <- proj4string(res2)
    res2 <- gIntersection(res2, ch, byid = TRUE)
  }
  res2 <- sp::spChFIDs(res2, as.character(unique(df[[ID]])))
  res2 <- sp::SpatialPolygonsDataFrame(res2, data.frame(ID = names(res2), row.names = names(res2), stringsAsFactors = FALSE))
  return(res2)
}

#' Découper les polygones de Voronoi pour respecter les limites d'un polygone
#'
#' @param voronoi_shp polygones de Voronoi (objet \code{SpatialPolygons*})
#' @param outline contour à respecter (objet \code{SpatialPolygons*})
#'
#' @return un objet \code{SpatialPolygons}
#' @export
#'
crop_voronois <- function(voronoi_shp, outline) {
  if (is.na(proj4string(voronoi_shp))) {
    proj4string(voronoi_shp) <- CRS("+init=epsg:4326")
    warning("Attention, le shape voronoi_shp n'a pas de projection explicite. On suppose qu'il s'agit de coordonnées géographiques")
  }
  if (is.na(proj4string(outline))) {
    warning("Attention, le shape outline n'a pas de projection explicite. On suppose qu'il s'agit de coordonnées géographiques")
    proj4string(outline) <- CRS("+init=epsg:4326")
  } else if (!stringr::str_detect(proj4string(outline), "longlat")) {
    outline <- spTransform(outline, CRS("+init=epsg:4326"))
    message("Transformation du shape en coordonnées géographiques")
  }
  ids <- names(voronoi_shp)
  res <- rgeos::gIntersection(voronoi_shp, outline, byid = TRUE)
  res <- spChFIDs(res, stringr::str_split_fixed(names(res), " ", n = 2)[,1])
  res <- sp::SpatialPolygonsDataFrame(res, data.frame(ID = stringr::str_split_fixed(names(res), " ", n = 2)[,1], row.names = stringr::str_split_fixed(names(res), " ", n = 2)[,1], stringsAsFactors = FALSE))
  res <- sp::spChFIDs(res, ids)
  return(res)
}
