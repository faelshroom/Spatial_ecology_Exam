- library("rgbif") #Downloading raw GBIF data.

- library("sf") #To transform data frames into spatial objects, and other useful functions.

- library("spatstat") #To calculate KDE and ppp. \Lcross??

- library("viridis") #To use color-blind friendly color palettes in the plots.

- library("rnaturalearth") #provides Australia's borders used as the observation window.

- library("ggplot2") #To visualize the occurrence and density in Australia. 

-library ("taxize") #to use the backbone function

# 1. Spatial Setup

#We load the Australia map and transform it into a metric system, better for calculus.
australia <- ne_countries (
country = "Australia",
 scale = "medium",      
returnclass = "sf"
) |>
  st_transform(3577)

#We transform the map into the observation window used later for the ppp

australia_poly <- as.owin(australia)

#We visualize the Australia map without occurrence data.

ggplot() + geom_sf(
data = australia,
fill = "#f8f9fa",          
color = "grey80", 
linewidth = 0.2
) + theme_minimal() + theme(panel.grid = element_blank())


# 2. Data Loading Function

#We retireve the occurrence data from GBIF. 


load_species_sf <- function(taxonKey, exclude_pattern = NULL) {   

#We put he downloaded data into a table with $data, and filter it for the  coordinates in Australia.

data <- occ_search(
  taxonKey = taxonKey,
  country = "AU",
  hasCoordinate = TRUE, 
  limit = 10000

)$data

#We keeep only the columns useful for the analysis. 

data <- data[, c(
  "decimalLongitude", 
  "decimalLatitude", 
  "scientificName"
)]

#Filter to cancel records with no coordinates.
 
data <- data[
 !is.na(data$decimalLongitude) &
 !is.na(data$decimalLatitude), 
]

#We deliberately choose not to exclude domestic animals as they also are fundamental to prove our point
 
#We remove duplicated coordinates to avoid repeated records to produce artifical density hotspots

data <- data[
 !duplicated(
  data[, c("decimalLongitude","decimalLatitude")]
 ),
]

  #We convert now the table into an sf spatial object, and again transform it in metric units. 

sf_points <- st_as_sf(
  data,
  coords = c("decimalLongitude", "decimalLatitude"), 
  crs = 4326
) |>
  st_transform(3577)

#Finally, we only keep the points inside the Australian border with a vector. 
#Sparse needs to be false for it to return a logical vector.

  return(
    sf_points[
     st_intersects(
       sf_points,
       australia, 
       sparse = FALSE
     ),
    ]
  )
}   

# 3. Data Acquisition

#We apply the function and download both species data, first we find the numeric keys and don't assume them manually, as they can change on GBIF

cat_taxon <- name_backbone(
  name= "Felis catus"
)

quokka_taxon <- name_backbone (
  name = "Setonix barchyurus"
)

cat_taxon$usageKey
quokka_taxon$usagekey

feral_cat_sf <- load_species_sf(
  cat_taxon$usageKey
)

quokka_sf <- load_species_sf(
  quokka_taxon$usageKey
)

# 4. Density Calculation (Sigma: 50km)

#We create our ppp objects necessary for the KDE. We extract X and Y coordinates in meters from the spatial object, 1 being the longitude and 2 the latitude.
#Then, we also define the observation window created before, our australia polygon. We do this process for cat and quokka.

cat_ppp <- ppp(
  st_coordinates(feral_cat_sf)[,1],
  st_coordinates(feral_cat_sf)[,2],
  window = australia_poly
)

quokka_ppp <- ppp(
  st_coordinates(quokka_sf)[,1], 
  st_coordinates(quokka_sf)[,2], 
  window = australia_poly
)

#With the ppp objects ready, we can now compute the density calculus. We use a common sigma of 50km for both species, with a 512x512 grid. 

cat_dens <- density(
  cat_ppp, 
  sigma = 50000, 
  dimyx = 512
)

quokka_dens <- density(
  quokka_ppp, 
  sigma = 50000, 
  dimyx = 512
)

# 5. Log-Normalization

#We create the Log-Normalization function for our densities, to visualize and compare them better in the plots.

apply_log_norm <- function(dens_obj) {

  #We create and add a small offset to add to every pixel, to avoid log(0) problem, so every pixel has a value.
  #na.rm ignores any NA values. 

  offset <- max(
    dens_obj$v, 
    na.rm = TRUE
  ) / 1000

  # Log transformation.
  dens_obj$v <- log(
    dens_obj$v + offset
  )

  #We apply a Min-Max scaling for the normalization so that every value is contained between 1.0 and 0.0 for both density scales. 

  dens_obj$v <- (
    dens_obj$v - 
      min(dens_obj$v, na.rm=TRUE) #true o t?
  ) / 
    (
      max(dens_obj$v, na.rm=TRUE) - 
        min(dens_obj$v, na.rm=TRUE)
    )
  return(dens_obj)

}

#Apply the Log-Normalization to our densities. 

cat_dens_log <- apply_log_norm(cat_dens)
quokka_dens_log <- apply_log_norm(quokka_dens)

# 6. Plotting Functions

#The first plotting function is used to plot the occurrence points in Australia.

plot_occ <- function(
  sf_points, 
  species_label, 
  color_p
) {

  ggplot() +

    #We first draw the Australian map to put the points into.
    geom_sf(
      data = australia, 
      fill = "#f8f9fa", 
      color = "grey80", 
      linewidth = 0.2
    ) +

    #We now draw the points, taking the data from our sf object created before.

    geom_sf(
      data = sf_points, 
      color = color_p, 
      size = 0.3, 
      alpha = 0.4
    ) + labs(
      title = paste(
        "Occurrences:", 
        species_label
      )
     ) + 

    #We remove the default grey background and grid lines of R, making the map look cleaner. 

    theme_minimal() + theme(
      panel.grid = element_blank()
    )
}

#The second plotting function is used to plot the density estimations in Aunstralia with a "heatmap". 

plot_dens <- function(
  dens_obj, 
  species_label, 
  palette
) {

  #The density function creates an image, but we need a data table for ggplot2. We use the as.data.frame function to obtain it. 

  df <- as.data.frame(dens_obj)

  #We put names onto the columns. 

  colnames(df) <- c(
    "x", 
    "y",  
    "value"
  )
  
  ggplot() +

    #We first draw the shape of Australia, filled with a light grey, acting as a background, so areas with zero density still look like part of the country.

    geom_sf(
      data = australia, 
      fill = "#eeeeee", 
      color = NA
    ) +

    #We then draw our grid of pixels, the aes function tells R that the color of the pixel should be determined by the density number (0 to 1). Alpha is at 85% so you can still see the map.

    geom_raster(
      data = df, 
      aes(
        x = x, 
        y = y, 
        fill = value
      ), 
      alpha = 0.85
    ) +

    #We set our viridis color palette

    scale_fill_viridis(
      option = palette, 
      name = "Log Density"
    ) +

    #We draw the australian border again, but this time with fill NA, just putting the outline on top of the heatmap so the borders look sharp. 

    geom_sf(
      data = australia, 
      fill = NA, 
      color = "white", 
      linewidth = 0.1
    ) + labs(
      title = paste(
        "KDE:", 
        species_label
      ), 
      subtitle = "Sigma: 50km",
      x = "Longitude"
      y = "Latitude"
    ) +

    #We remove the default grey backgorund and grid lines of R, to make the map look cleaner. 

    theme_minimal() + 
    theme(
      panel.grid = element_blank()
    )
    aspect.ratio = 1.3
}

# 7. Final Layout

#We create all 4 plots

p1 <- plot_occ(
  feral_cat_sf, 
  "Feral Cat", 
  "#440154"
)

p2 <- plot_dens(
  cat_dens_log, 
  "Feral Cat", 
  "viridis"
) 

p3 <- plot_occ(
  quokka_sf, 
  "quokka", 
  "#35b779"
)

p4 <- plot_dens(
  quokka_dens_log, 
  "quokka", 
  "magma")

#Finally, we print the maps

p1
p2
p3
p4


# 8. Statistical Analysis
#We compute our Spearman Analysis. We convert our density matrix into a single column vector, so each grid cell becomes a single observation.
 
#We use complete.obs to ignore NA values. 

#We calculate the correlation coefficient between the two vectors with the Spearman method. 

spearman_rho <- cor(
 as.vector(cat_dens_log$v),
 as.vector(quokka_dens_log$v), 
 method = "spearman", 
 use = "complete.obs"
)

#We round our result to 2 decimal places, and then see the result. 
print(
  paste(
    "Spearman Correlation:", 
    round(spearman_rho, 2)
  )
)

Spearman correlation: 0.37

# 9. Density Difference Map

#We create the Density Difference Map, showing areas where the quokka or the feral cat are mainly present. 

#We again turn the feral cat density matrix into a data table, giving it specific names. 

diff_df <- as.data.frame(
  cat_dens_log
)

colnames(diff_df) <- c(
  "x",
  "y",
  "cat_val"
)

#Since both the feral cat and quokka density maps were created using the exact same grid, we can simply "paste" the quokka density values as a new column in our table. They line up perfectly pixel-for-pixel.

diff_df$quokka_val <- 
 as.data.frame(
 quokka_dens_log
)$value

#Since both variables are normalized, here is where we calculate the density. With this order, positive values will mean higher feral cat density, and negative will mean higher quokka density.

diff_df$diff <- 
  diff_df$cat_val -
  diff_df$quokka_val

#We use the same logic used for the normal density maps, just changing the data used. 

ggplot() + geom_sf(
    data = australia,
    fill = "grey95", 
    color = NA
  ) + geom_raster(
    data = diff_df, 
    aes(
     x = x,
     y = y,
     fill = diff
    ), 
    alpha = 0.9
  ) + scale_fill_viridis_c(
    option = "mako", 
    name = "Difference"
  ) + geom_sf(
    data = australia, 
    fill = NA, 
    color = "white",
    linewidth = 0.1
  ) + labs(
    title = "Relative Spatial Intensity: Feral Cat vs Quokka", 
    subtitle = 
      "Brighter = Higher cat intensity | Negative = Higher quokka intensity"
  ) + theme_minimal() + 
  theme(panel.grid = element_blank()
)


