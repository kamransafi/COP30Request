#' Download and annotate Copernicus DSM (30 m) for a move2 object
#'
#' This function:
#' \itemize{
#'   \item Expands the spatial extent (bounding box) of a \code{move2} object.
#'   \item Requests a Copernicus DEM GLO-30 (COPERNICUS_30) DSM subset from the
#'         Copernicus Data Space Ecosystem via the Sentinel Hub Process API.
#'   \item Extracts elevation values at each animal location (nearest neighbour
#'         or bilinear interpolation).
#'   \item Optionally plots a quick overview map (DEM + locations coloured by track ID).
#' }
#'
#' Authentication is done via the \pkg{CDSE} package using a Sentinel Hub
#' OAuth client ID and secret stored securely in the system keyring via
#' the \pkg{keyring} package (entries named \code{"sentinelhub_client_id"}
#' and \code{"sentinelhub_client_secret"} under the chosen \code{key_service}).
#'
#' @param x A \code{move2} object containing animal trajectory data.
#' @param extent_factor Numeric > 0. Factor by which the bounding box of \code{x}
#'   should be expanded before requesting the DSM (1 = original extent, 1.1 = +10\%).
#' @param method Character; interpolation method used by \code{terra::extract()}.
#'   One of \code{"nearest"} (nearest neighbour; internally \code{"simple"}) or
#'   \code{"bilinear"} (bilinear interpolation).
#' @param plot_overview Logical; if \code{TRUE} (default), produce a base R
#'   overview plot with the DSM as background and \code{move2} points coloured
#'   by track ID.
#' @param dem_file Optional character file path for the output GeoTIFF.
#'   If \code{NULL} (default), a temporary file is created.
#' @param key_service Character; name of the keyring service under which the
#'   Sentinel Hub client ID and secret are stored. Default is \code{"SentinelHUB"}.
#'
#' @return An (invisible) list with three elements:
#'   \itemize{
#'     \item \code{move2}: the input \code{move2} object with a new numeric column
#'           \code{dem_copernicus30} containing DSM elevations (in metres) at
#'           each location (NA where geometry is missing or outside DEM).
#'     \item \code{dem}: a \code{terra::SpatRaster} containing the downloaded DSM.
#'     \item \code{demfile}: character path to the GeoTIFF file on disk.
#'   }
#'
#' @details
#' The function assumes that the DEM is requested and returned in WGS84
#' longitude/latitude (EPSG:4326) and therefore does not reproject the points
#' when extracting. If you change the DEM CRS in the Process API request,
#' you must adapt the extraction step accordingly.
#'
#' This function is intended as a convenience wrapper for exploratory analysis
#' and not as a production-grade client. Use with appropriate API quotas
#' and follow Copernicus Data Space terms of use.
#'
#' @examples
#' \dontrun{
#' library(move2)
#'
#' # Example with bundled move2 data
#' fishers <- mt_read(mt_example())
#'
#' res <- get_dsm_for_move2(
#'   x             = fishers,
#'   extent_factor = 1.2,
#'   method        = "bilinear",
#'   plot_overview = TRUE,
#'   key_service   = "SentinelHUB"
#' )
#'
#' fishers_annotated <- res$move2
#' head(fishers_annotated$dem_copernicus30)
#' }
#'
#' @export

library(keyring)
library(CDSE)
library(httr2)
library(terra)
library(sf)
library(move2)
library(dplyr)

# -------------------------------------------------------------------
# Package setup: install missing packages, then load
# -------------------------------------------------------------------
required_pkgs <- c(
  "keyring",
  "CDSE",
  "httr2",
  "terra",
  "sf",
  "move2",
  "dplyr"
)

missing_pkgs <- required_pkgs[!required_pkgs %in% installed.packages()[, "Package"]]

if (length(missing_pkgs) > 0) {
  message("Installing missing packages: ",
          paste(missing_pkgs, collapse = ", "))
  install.packages(missing_pkgs)
}

# load all packages (stop with a clear error if something fails)
invisible(lapply(required_pkgs, function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop("Package '", pkg, "' is required but could not be loaded.")
  }
  library(pkg, character.only = TRUE)
}))


.dem_evalscript <- "
//VERSION=3
function setup() {
  return {
    input: [\"DEM\"],
    output: {
      id: \"default\",
      bands: 1,
      sampleType: SampleType.FLOAT32
    }
  };
}
function evaluatePixel(sample) {
  return [sample.DEM];
}
"

get_dsm_for_move2 <- function(x,
                              extent_factor = 1.1,
                              method = c("nearest", "bilinear"),
                              plot_overview = TRUE,
                              dem_file = NULL,
                              key_service = "SentinelHUB") {
  method <- match.arg(method)
  
  if (!inherits(x, "move2")) {
    stop("x must be a move2 object.")
  }
  
  # 0. Track IDs from move2 (independent of column naming)
  track_ids_event <- mt_track_id(x)
  
  # 1. Reproject to WGS84 and build expanded bbox
  x_wgs <- st_transform(x, 4326)
  
  bb <- st_bbox(x_wgs)
  cx <- (bb["xmin"] + bb["xmax"]) / 2
  cy <- (bb["ymin"] + bb["ymax"]) / 2
  dx <- (bb["xmax"] - bb["xmin"]) * extent_factor / 2
  dy <- (bb["ymax"] - bb["ymin"]) * extent_factor / 2
  
  bbox_vec <- c(cx - dx, cy - dy, cx + dx, cy + dy)
  
  # 2. Get Sentinel Hub token via CDSE using keyring credentials
  sh_client_id  <- key_get("sentinelhub_client_id",     service = key_service)
  sh_client_sec <- key_get("sentinelhub_client_secret", service = key_service)
  
  if (is.null(sh_client_id) || is.null(sh_client_sec)) {
    stop("Sentinel Hub credentials not found in keyring under service '",
         key_service, "'.")
  }
  
  token <- GetOAuthToken(id = sh_client_id, secret = sh_client_sec)
  
  # 3. Process API body for COPERNICUS_30 DEM
  process_body <- list(
    input = list(
      bounds = list(
        properties = list(
          crs = "http://www.opengis.net/def/crs/OGC/1.3/CRS84"
        ),
        bbox = as.numeric(bbox_vec)
      ),
      data = list(list(
        type = "dem",
        dataFilter = list(
          demInstance = "COPERNICUS_30"
        ),
        processing = list(
          upsampling   = "BILINEAR",
          downsampling = "BILINEAR"
        )
      ))
    ),
    output = list(
      resx = 0.0003,
      resy = 0.0003,
      responses = list(list(
        identifier = "default",
        format = list(
          type = "image/tiff"
        )
      ))
    ),
    evalscript = .dem_evalscript
  )
  
  process_url <- "https://sh.dataspace.copernicus.eu/api/v1/process"
  
  # 4. Call Process API, save DEM, read raster
  resp <- request(process_url) |>
    req_headers(
      Authorization  = paste("Bearer", token),
      `Content-Type` = "application/json"
    ) |>
    req_body_json(process_body, auto_unbox = TRUE) |>
    req_perform()
  
  resp_check_status(resp)
  
  if (is.null(dem_file)) {
    dem_file <- tempfile(fileext = ".tif")
  }
  writeBin(resp_body_raw(resp), dem_file)
  
  dem_rast <- rast(dem_file)
  
  # 5. Extract DEM values at locations (robust to missing geometry)
  sf_pts <- st_as_sf(x_wgs)
  
  # identify rows with valid point geometry
  coords <- st_coordinates(sf_pts)
  valid  <- !is.na(coords[, "X"]) & !is.na(coords[, "Y"])
  
  # initialise with NA_real_
  dem_vals <- rep(NA_real_, nrow(sf_pts))
  
  if (any(valid)) {
    pts_valid <- vect(sf_pts[valid, ])
    
    method_terra <- if (method == "nearest") "simple" else "bilinear"
    dem_valid    <- extract(dem_rast, pts_valid, method = method_terra)[, 2]
    
    dem_vals[valid] <- dem_valid
  }
  
  x$dem_copernicus30 <- dem_vals
  
  # 6. Optional overview plot: DEM + move2 points coloured by track id
  if (plot_overview) {
    op <- par(no.readonly = TRUE)
    on.exit(par(op), add = TRUE)
    
    plot(dem_rast,
         main = "Copernicus DSM (GLO-30) with tracks",
         col  = terrain.colors(64))
    
    track_ids <- as.factor(track_ids_event)
    cols <- rainbow(length(levels(track_ids)))[track_ids]
    
    plot(st_geometry(x_wgs),
         add = TRUE,
         col = cols,
         pch = 16,
         cex = 0.6)
    
    legend("topright",
           legend = levels(track_ids),
           col    = rainbow(length(levels(track_ids))),
           pch    = 16,
           bg     = "white",
           cex    = 0.8)
  }
  
  invisible(list(
    move2   = x,
    dem     = dem_rast,
    demfile = dem_file
  ))
}
