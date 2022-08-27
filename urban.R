install.packages("sf")
install.packages("rlang") #Para los paquetes que fueron desarrollados en versiones distintas
library(rlang)
library(mapedit)
library(raster)
library(mapview)
library(dplyr)
library(rgee)
library(sp)
library(sf)

#Dibujar nuestro ambito de estudio
drawing <- drawFeatures()
mapview(drawing)
#Inicializar Earth Engine
ee_Initialize("Junior")
draw_ee <- drawing %>% sf_as_ee()
Map$centerObject(draw_ee)
Map$addLayer(draw_ee)

#Funcion de factor de escala para imagenes Surface Reflectance (SR) Landsat 5
applyScaleFactorsl5 <- function(image) {
  opticalBands <- image$select('SR_B.')$multiply(0.0000275)$add(-0.2)
  thermalBand <- image$select('ST_B6')$multiply(0.00341802)$add(149.0)
  return(image$addBands(opticalBands, NULL, TRUE)$
           addBands(thermalBand, NULL, TRUE))
}

#Funcion de factor de escala para imagenes Surface Reflectance (SR) Landsat 8
applyScaleFactorsl89 <- function(image) {
  opticalBands <- image$select('SR_B.')$multiply(0.0000275)$add(-0.2)
  thermalBand <- image$select('ST_B.*')$multiply(0.00341802)$add(149.0)
  return(image$addBands(opticalBands, NULL, TRUE)$
           addBands(thermalBand, NULL, TRUE))
}


#Landsat 5 TM Collection 2 reflectancia superficial corregida atmosféricamente.

# Año 1990 
img1991 <- ee$ImageCollection('LANDSAT/LT05/C02/T1_L2')$
  filterDate('1991-01-01', '1992-01-01')$
  filterBounds(draw_ee)$
  filterMetadata('CLOUD_COVER', 'less_than', 10)$
  map(applyScaleFactorsl5)$
  median()$
  clip(draw_ee)

# Año 2000 
img1998 <- ee$ImageCollection('LANDSAT/LT05/C02/T1_L2')$
  filterDate('1998-01-01', '1999-02-01')$
  filterBounds(draw_ee)$
  filterMetadata('CLOUD_COVER', 'less_than', 10)$
  map(applyScaleFactorsl5)$
  median()$
  clip(draw_ee)

# Parametros de visualizacion
visparal5 <- list(
  bands = c('SR_B4', 'SR_B5', 'SR_B3'),
  min = 0.1,
  max = 0.3
)
Map$addLayer(img1991, visparal5, 'img1991')|
  Map$addLayer(img1998, visparal5, 'img1998')
#Landsat 8 OLI/TIRS Collection 2 reflectancia superficial corregida atmosféricamente

# Año 2015 
img2016 <- ee$ImageCollection('LANDSAT/LC08/C02/T1_L2')$
  filterDate('2016-01-01', '2017-01-01')$
  filterBounds(draw_ee)$
  filterMetadata('CLOUD_COVER', 'less_than', 10)$
  map(applyScaleFactorsl89)$
  median()$
  clip(draw_ee)

# Parametros de visualizacion
visparal89 <- list(
  bands = c('SR_B5', 'SR_B6', 'SR_B4'),
  min = 0.1,
  max = 0.3
)
Map$addLayer(img2016, visparal89, 'img2016')

# Landsat 9 OLI-2/TIRS-2 Colección 2 Reflectancia superficial corregida atmosféricamente.

# Año 2021 
img2021 <- ee$ImageCollection('LANDSAT/LC08/C02/T1_L2')$
  filterDate('2021-01-01', '2022-01-01')$
  filterBounds(draw_ee)$
  filterMetadata('CLOUD_COVER', 'less_than', 10)$
  map(applyScaleFactorsl89)$
  median()$
  clip(draw_ee)
Map$addLayer(img2016, visparal89, 'img2016')|
Map$addLayer(img2021, visparal89, 'img2021')

# function for calculating built-up area
funbuil <- function(img, swir2, swir1, nir, red, green, blue) {
  # Calculating NDVI
  ndvi <- img$normalizedDifference(c(nir, red))
  ndvi_gt <- ndvi$gt(0.15)
  ndvi_mask <- ndvi_gt$updateMask(ndvi_gt)
  
  # Calculating NDBI
  ndbi <- img$normalizedDifference(c(swir2, nir))
  ndbi_gt <- ndbi$gt(0)
  ndbi_mask <- ndbi_gt$updateMask(ndbi_gt)
  
  # Calculating MNDWI
  mndwi <- img$normalizedDifference(c(green, nir))
  mndwi_gt <- mndwi$gt(0)
  mndwi_mask <- mndwi_gt$updateMask(mndwi_gt) 
  
  # Calculating BI
  bi <- img$expression(
    "(SWIR + RED) − (NIR + BLUE) / (SWIR + RED) + (NIR + BLUE)", list(
      "SWIR" = img$select(swir1),
      "RED" = img$select(red),
      "NIR" = img$select(nir),
      "BLUE" = img$select(blue)
    )
  )
  bi_gt <- bi$gt(0)
  bi_mask <- bi_gt$updateMask(bi_gt)
  
  # Blend MNDWI and BI
  ble <- mndwi_mask$blend(bi_mask)$blend(ndvi_mask)$selfMask()
  
  # Remove blend areas of MNDWI and BI of NDVI without vegetation
  built <- ndbi_mask$updateMask(ble$unmask()$Not())$selfMask()
  
  return(built)
}

# Calculating Built-up 

buil_1991 <- funbuil(
  img = img1991,
  swir2 = "SR_B7",
  swir1 = "SR_B5", 
  nir = "SR_B4",
  red = "SR_B3",
  green = "SR_B2",
  blue = "SR_B1"
)$rename("built")
buil_1998 <- funbuil(
  img = img1998,
  swir2 = "SR_B7",
  swir1 = "SR_B5", 
  nir = "SR_B4",
  red = "SR_B3",
  green = "SR_B2",
  blue = "SR_B1"
)$rename("built")
buil_2016 <- funbuil(
  img = img2016,
  swir2 = "SR_B7",
  swir1 = "SR_B6", 
  nir = "SR_B5",
  red = "SR_B4",
  green = "SR_B3",
  blue = "SR_B2"
)$rename("built")
buil_2021 <- funbuil(
  img = img2021,
  swir2 = "SR_B7",
  swir1 = "SR_B6",
  nir = "SR_B5",
  red = "SR_B4",
  green = "SR_B3",
  blue = "SR_B2"
)$rename("built")

Map$addLayer(buil_1991, list(palette="blue"), "Built-up 1991")|
  Map$addLayer(buil_1998, list(palette="red"), "Built-up 1998")

Map$addLayer(buil_2016, list(palette="blue"), "Built-up 2016")|
Map$addLayer(buil_2021, list(palette="red"), "Built-up 2021")

#functions to Calculating area km2
funarea <- function(imgbuil) {
  #Multiplicacion matricial a la clasificacion
  areaimage <- imgbuil$multiply(ee$Image$pixelArea())
  
  #sumamos todos los valores de la región 
  area <- areaimage$reduceRegion(
    reducer = ee$Reducer$sum(),
    geometry = draw_ee$geometry(),
    scale = 30,
    maxPixels = 1e9
  )
  
  #Obtencion de las area de cambio en Km2
  clasAreaKm2 <- ee$Number(area$get('built'))$divide(1e6)
  
  return(ee$Number$getInfo(clasAreaKm2))
}
funarea(buil_2021)
areabuilt <- c(6.16, 7.27, 13.90, 14.83)

#Reducer raster data to vector and  Earth engine to Local

funlocal <- function(ras, dic, nom) {
  # Reduce to vector from raster date
  vect <- ras$reduceToVectors(
    reducer = ee$Reducer$countEvery(),
    geometry = draw_ee, 
    scale = 30,
    maxPixels = 1e12
  )
  
  # Vector data from Earth Engine to local
  loc <- ee_as_sf(x = vect, dsn = dic, via = "drive", container = nom)
  return(loc)
}
 
loc_1991 <- funlocal(buil_1991, "Resultados/built_1991.shp", "built_1991")
loc_1998 <- funlocal(buil_1998, "Resultados/built_1998.shp", "built_1998")
loc_2016 <- funlocal(buil_2016, "Resultados/built_2016.shp", "built_2016")
loc_2021 <- funlocal(buil_2021, "Resultados/built_2021.shp", "built_2021")

# Visualizing result data con ggplot2 
library(ggplot2)

# Function to mapping with ggplot2 
funggplot <- function(loc, are, titl ) {
  ggpl <- ggplot()+
    geom_sf(data = loc,  aes(fill = are), 
            pch = 21,    
            col = "red", # Color del borde
            cex = 0)+
    theme_bw()+
    labs(x='Logitud', y='Latitud', title = titl, subtitle = "Huancayo", fontface = "bold")+
    scale_fill_discrete(name = "Area")+
    theme(plot.title = element_text(size=15))
  
  return(ggpl)
}

# Plotting built up of all year 
funggplot(loc = loc_1991, are = "6.16 Km2", titl = "Areas Construidas 1991")
funggplot(loc = loc_1998, are = "7.27 Km2", titl = "Areas Construidas 1998")
funggplot(loc = loc_2016, are = "13.90 Km2", titl = "Areas Construidas 2016")
funggplot(loc = loc_2021, are = "14.83 Km2", titl = "Areas Construidas 2021")
