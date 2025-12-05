# Copernicus DSM for move2 tracks

This repository provides an R helper function `get_dsm_for_move2()` to:

- Expand the spatial extent of a `move2` trajectory object.
- Request Copernicus DEM GLO‑30 (DSM, 30 m) from the Copernicus Data Space Ecosystem (CDSE) via the Sentinel Hub Process API.
- Extract elevation values at each animal location (nearest neighbour or bilinear interpolation).
- Optionally plot an overview map (DSM + points colored by track ID).

The function is intended for exploratory work with movement data and elevation.

***

## Installation

Clone this repository and source the main R script:

```r
source("R/get_dsm_for_move2.R")
```

The script will check for required packages and install any that are missing:

- keyring
- CDSE
- httr2
- terra
- sf
- move2
- dplyr

You can also install them manually beforehand if you prefer.

***

## Getting CDSE / Sentinel Hub credentials

To use this function you need an OAuth client in the Copernicus Data Space Ecosystem with access to Sentinel Hub APIs. The function uses that client’s ID and secret to obtain access tokens.

### 1. Create an account and log in

1. Go to the Copernicus Data Space Ecosystem portal (dataspace.copernicus.eu).
2. Create a user account if you do not have one yet, and sign in.

### 2. Open the Sentinel Hub dashboard

1. In the CDSE portal, open the **Sentinel Hub** dashboard (e.g. from the left sidebar or “Dashboards” menu).
2. This opens the Sentinel Hub user interface in a new tab.

### 3. Create an OAuth client

In the Sentinel Hub UI:

1. Go to **User settings** or **Account** (top‑right user menu).
2. Find the **OAuth clients** (or **API clients**) section.
3. Click **Create new client** (or similar), and provide:
    - A descriptive name (e.g. `move2-dsm-client`).
    - Optional description.
    - Redirect URI can be a dummy value if you use the client‑credentials flow only (for script use).
4. After creation, the interface will show:
    - **Client ID**
    - **Client Secret**

Copy both immediately and store them securely. You normally only see the secret once.

### 4. Store credentials in the R keyring

The function expects the client ID and secret to be stored in your system keyring under the service name `"SentinelHUB"` with keys:

- `"sentinelhub_client_id"`
- `"sentinelhub_client_secret"`

In R:

```r
library(keyring)

key_set("sentinelhub_client_id",     service = "SentinelHUB")
key_set("sentinelhub_client_secret", service = "SentinelHUB")
```

You will be prompted to paste each value; they are then held by the OS credential store, not in your code.

You only need to do this once per machine/user.

***

## Function overview

```r
get_dsm_for_move2(
  x,
  extent_factor = 1.1,
  method        = c("nearest", "bilinear"),
  plot_overview = TRUE,
  dem_file      = NULL,
  key_service   = "SentinelHUB"
)
```

- `x` – a `move2` object.
- `extent_factor` – factor to expand the bounding box before requesting DEM (1.1 = 10% larger in each dimension).
- `method` – `"nearest"` (nearest neighbour) or `"bilinear"` (bilinear interpolation) for elevation extraction.
- `plot_overview` – if `TRUE`, plots a quick base R map of DSM and track points.
- `dem_file` – optional path to write the DEM GeoTIFF; if `NULL`, a temporary file is used.
- `key_service` – keyring service name; default `"SentinelHUB"`.

Return value (invisible list):

- `$move2` – input `move2` object with new column `dem_copernicus30` (elevation in metres).
- `$dem` – `terra::SpatRaster` DSM.
- `$demfile` – path to the DSM GeoTIFF written to disk.

***

## Example usage

```r
library(move2)

# Example data from move2
fishers <- mt_read(mt_example())

# Get DSM and annotate locations
res <- get_dsm_for_move2(
  x             = fishers,
  extent_factor = 1.2,
  method        = "bilinear",
  plot_overview = TRUE,
  key_service   = "SentinelHUB"
)

fishers_annotated <- res$move2
head(fishers_annotated$dem_copernicus30)

# The DEM SpatRaster
res$dem
```


***

## Notes and limitations

- This script is aimed at small to medium spatial extents; for large areas, consider tiling and caching strategies.
- The Process API and DEM access are subject to CDSE/Sentinel Hub quotas and terms of use; monitor your usage.
- Elevation extraction returns `NA` where:
    - the original `move2` record has no geometry, or
    - the point falls outside the requested DEM bbox.

