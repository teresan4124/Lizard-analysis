# Lizard-analysis
This is the code that I used to analyze different skull measurements of 102 different lizards in the Sceloporus species.
# Lizard Cranial Morphometrics Analysis Using Geometric Morphometrics in R

This project analyzes 3D landmark data of **Sceloporus** lizard skulls using geometric morphometric techniques in R. The analysis involves cleaning, aligning, and statistically comparing morphological traits between arboreal and non-arboreal species.

## Overview

The code performs the following key steps:

- Loads 3D landmark datasets (`.nts` files) from CT scans
- Performs **Generalized Procrustes Analysis (GPA)** to align shapes
- Extracts Procrustes coordinates, centroid size, and shape data
- Identifies potential outliers
- Merges landmark data with metadata (`scelinfoV2.csv`)
- Calculates linear distances between landmarks (e.g., snout length, jaw width)
- Normalizes measurements by body size (SVL)
- Performs group-wise visualizations and statistical comparisons (e.g., t-tests between arboreal and non-arboreal lizards)

## Tools and Packages

This analysis uses the following R packages:

- `geomorph` — for shape analysis and GPA
- `mvMORPH`, `PCDimension` — for multivariate morphometrics
- `ips`, `phytools`, `geiger` — for phylogenetic context (if needed)
- `scatterplot3d`, `car` — for visualization and statistical tools

## Features Measured

The following cranial traits are measured as distances between key 3D landmarks:

| Code | Trait Description             |
|------|-------------------------------|
| A    | Snout length                  |
| B    | Height at coronoid            |
| C    | Snout width                   |
| D    | Closing-in lever (jaw)        |
| E    | Opening-in lever (jaw)        |
| F    | Jaw width                     |
| G    | Lower jaw length (A + D)      |

These measurements are computed both as raw distances and size-corrected values (e.g., A/SVL).

## Statistical Analysis

- **Outlier detection**: Using Procrustes coordinates
- **T-tests**: Comparing scaled trait values between arboreal and non-arboreal species
- **Quartile analysis**: Visual inspection of shape variation

## Output Files

- `nosize.csv` — Raw trait measurements per species
- `wsize.csv` — Size-corrected (SVL-scaled) measurements per species

## Notes

- Ensure all `.nts` landmark files and the metadata file `scelinfoV2.csv` are in the same working directory.
- Update species names in metadata to match those extracted from the 3D data files exactly.

## Example Output

```r
plot(AS ~ arb) # Visualize snout length (scaled) by habitat type
t.test(AS ~ arb) # Statistical comparison between arboreal vs non-arboreal species
