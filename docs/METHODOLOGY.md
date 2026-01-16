# Methodology: Bayesian Groundwater Storage Estimation

## Overview

GWLevel_USA implements a Bayesian framework for estimating groundwater storage changes across the conterminous United States (CONUS) using in-situ well observations. The methodology combines:

1. **Historical baseline establishment** (pre-development conditions)
2. **Spatial bias correction** (IDW interpolation)
3. **Monte Carlo uncertainty quantification** (literature-based Sy priors)

---

## 1. Data Sources and Harmonization

### Input Data
- **USGS NWIS**: National Water Information System (~80% of data)
- **CA DWR**: California Department of Water Resources
- **TX TWDB**: Texas Water Development Board

### Harmonization
All sources converted to common format:
- Coordinates: WGS84 (EPSG:4326)
- Depth units: meters below land surface
- Temporal resolution: annual (Winter/Spring measurements preferred)

### Quality Filters
```
- Latitude: 24°N to 50°N (CONUS bounds)
- Longitude: 125°W to 66°W
- Depth: 0 < depth < 1000 m
- Minimum observations per well: 5
- Valid seasons: Winter, Spring (pre-irrigation)
```

---

## 2. Baseline Calculation

### Rationale
The baseline period (default: 1940-1955) represents pre-major development conditions:
- Before widespread irrigation expansion (High Plains)
- Before major municipal pumping increases
- Consistent with USGS historical assessments

### Method
For each well with observations in the baseline period:
```
Baseline_i = mean(Depth_{i,t}) for t ∈ [1940, 1955]
```

Requirements:
- Minimum 2 years of observations in baseline period
- Only wells meeting this criterion are included in analysis

---

## 3. Anomaly Calculation

### Convention (GRACE-compatible)
```
Anomaly = Baseline - Current_Depth
```

**Interpretation:**
- **Negative anomaly** → Water table is deeper → Depletion
- **Positive anomaly** → Water table is shallower → Recovery
- **Zero anomaly** → No change from baseline

This convention matches GRACE satellite observations where negative TWS anomalies indicate water loss.

---

## 4. Spatial Bias Correction

### Problem
Monitoring wells are not randomly distributed:
- Clustered in areas of concern (heavy pumping)
- Biased toward depleted regions
- Simple averaging overestimates depletion

### Solution: Inverse Distance Weighting (IDW)

IDW interpolation provides spatially-weighted estimates:

```
z_IDW = Σ(w_i × z_i) / Σ(w_i)

where:
  w_i = 1 / d_i^p
  p = 2 (default power parameter)
  d_i = distance from observation i to prediction point
```

### Implementation
1. Project coordinates to EPSG:5070 (NAD83/Conus Albers Equal-Area)
2. Build KD-tree for efficient neighbor search
3. Use k=12 nearest neighbors (configurable)
4. Compute weighted mean at aquifer centroid

### Spatial Correction Impact
Typical correction: 5-10% reduction in estimated depletion, indicating monitoring bias toward more depleted areas.

---

## 5. Storage Estimation

### Fundamental Equation
```
ΔS = Sy × A × Δh

where:
  ΔS = Storage change (km³)
  Sy = Specific yield (dimensionless)
  A  = Aquifer area (km²)
  Δh = Water level change (m)
```

### Unit Conversion
```
ΔS (km³) = Sy × Area (km²) × Δh (m) / 1000
```

---

## 6. Uncertainty Quantification

### Monte Carlo Approach
Generate N samples (default: 10,000) by sampling from:

1. **Specific Yield Distribution**
   ```
   Sy ~ TruncatedNormal(μ_Sy, σ_Sy, 0.0001, 0.5)
   ```
   - Lower bound: 0.0001 (physical minimum)
   - Upper bound: 0.5 (cannot exceed porosity)

2. **Anomaly Uncertainty** (optional)
   ```
   Δh ~ Normal(Δh_mean, Δh_sd)
   ```

3. **Storage Samples**
   ```
   ΔS_i = Sy_i × Area × Δh_i / 1000
   ```

### Summary Statistics
From the Monte Carlo samples:
- Mean: `mean(ΔS_samples)`
- SD: `std(ΔS_samples)`
- 95% CI: `[quantile(0.025), quantile(0.975)]`

---

## 7. Specific Yield Priors

### Literature-Based Values
Each aquifer assigned Sy prior from published USGS studies:

| Aquifer Type | Sy Range | Example Aquifers |
|-------------|----------|------------------|
| Unconfined alluvial | 0.15-0.25 | High Plains, Mississippi Valley |
| Basin fill | 0.12-0.20 | Basin and Range, Central Valley |
| Glacial deposits | 0.15-0.20 | Great Lakes, NE Glaciers |
| Fractured rock | 0.01-0.05 | Piedmont, Columbia Basalt |
| Confined | 0.0001-0.001 | Floridan, Carrizo-Wilcox |

### Key References
- Gutentag et al. (1984): High Plains Sy = 0.15 ± 0.04
- Faunt (2009): Central Valley Sy = 0.12 ± 0.04
- Pool & Coes (1999): Basin and Range Sy = 0.18 ± 0.04
- Clark & Hart (2009): Mississippi Valley Sy = 0.22 ± 0.05

---

## 8. Quality Control

### Outlier Detection
t-mixture model for robust outlier identification:

```
P(outlier | residual) = π₀ × t(ν_out) / [π₀ × t(ν_out) + (1-π₀) × t(ν_in)]

where:
  π₀ = prior outlier probability (0.01)
  ν_in = 7 (inlier degrees of freedom)
  ν_out = 1 (outlier, Cauchy-like)
```

### Source-Specific Priors
```
NWIS: π₀ = 0.005 (highest quality)
CA:   π₀ = 0.012
TX:   π₀ = 0.015
```

### Quality Tiers
Aquifers classified by data quality:

| Tier | CV | R² | Years | Wells |
|------|-----|-----|-------|-------|
| HIGH | <0.5 | >0.5 | ≥70 | ≥100 |
| MEDIUM | <1.0 | >0.2 | ≥50 | ≥30 |
| LOW | other | other | other | other |

---

## 9. Projection System

### EPSG:5070 (NAD83/Conus Albers Equal-Area)
Parameters:
```
Latitude of origin: 23°N
Central meridian: 96°W
Standard parallel 1: 29.5°N
Standard parallel 2: 45.5°N
Ellipsoid: GRS80
```

### Why Equal-Area?
- Accurate area calculations for storage estimation
- Standard for CONUS hydrological analysis
- Consistent with EPA/USGS conventions

---

## 10. Validation

### Comparison with Literature
| Region | This Study | Literature | Reference |
|--------|------------|------------|-----------|
| High Plains | -649 km³ | -410 km³ | McGuire (2017) |
| Basin & Range | -594 km³ | -270 km³ | Pool & Coes (1999) |

Differences explained by:
- Earlier baseline (1940 vs 1950)
- Extended analysis period (through 2025)
- Full aquifer extent vs regional studies

### GRACE Comparison
Recent rates (2000-2024) consistent with GRACE-derived estimates within uncertainty bounds.

---

## References

1. Clark, B.R., & Hart, R.M. (2009). The Mississippi Embayment Regional Aquifer Study (MERAS). USGS SIR 2009-5172.

2. Faunt, C.C. (2009). Groundwater Availability of the Central Valley Aquifer, California. USGS PP 1766.

3. Gutentag, E.D., et al. (1984). Geohydrology of the High Plains aquifer. USGS PP 1400-B.

4. Konikow, L.F. (2013). Groundwater depletion in the United States (1900-2008). USGS SIR 2013-5079.

5. McGuire, V.L. (2017). Water-level and recoverable water in storage changes, High Plains aquifer. USGS SIR 2017-5040.

6. Pool, D.R., & Coes, A.L. (1999). Hydrogeologic Investigations of the Sierra Vista Subwatershed. USGS WRIR 99-4197.

---

*Last updated: January 2026*
