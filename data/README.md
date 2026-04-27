# WCVI-E2E Input Data Documentation

This directory contains the physical and biological forcing datasets used by the WCVI-E2E (West Coast Vancouver Island End-to-End) ecosystem model.

All forcing datasets are preprocessed, formatted, and interpolated using functions in data/functions/ before being passed to the coupled physical–biogeochemical–trophic model.

---

## 1\. Physical Input Datasets

### 1.1 Atmospheric Forcing (NARR) (Folder physics/Forcing)

- **Variables**

|**Variables**|**Description**|**Units**|
|-|-|-|
|Qi\_input|Downward shortwave radiation flux|W m⁻²|
|airtmp\_input|Air temperature at 2 m|°C|
|dewptT\_input|Dew point temperature at 2 m|°C|
|uWndSpd\_input|Zonal wind speed at 10 m|m s⁻¹|
|vWndSpd\_input|Meridional wind speed at 10 m|m s⁻¹|
|p\_input|Precipitation rate|m s⁻¹|

- **Source**

North American Regional Reanalysis (NARR)
https://www.esrl.noaa.gov/psd/data/gridded/data.narr.html

- **Native temporal resolution**: 3-hourly

- **Native spatial resolution**: ~32 km gridded fields

- **Model usage**

Averaged spatially over latitude/longitude falling within WCVI shelf and slope boxes
Averaged temporally over 1992–2017 to construct annual climatologies with a 3 hour time step

- **Transition from _input.m to final input datasets used in model simulations**

See utility functions in folder /functions
Interpolated to ODE timestep (typically 30 minutes)
smooth using a 24-hour moving average (remove daily fluctuations) 

**Notes**

QO: clear sky irradiance calculated in processPhysicalForcing.m (Folder /functions)
WSPD10: wind velocity at 10 m above surface calculated in processPhysicalForcing.m (Folder /functions)

---

### 1.2 Wind stress (Folder physics/Forcing)

- **Variables**

|**Variables**|**Description**|**Units**|
|-|-|-|
|tauy2\_input|Alongshore (N–S) wind stress|N m⁻²|

- **Source**

Provided by collaborator Zelalem Engida (local daily values, 1995–2002)

- **Model usage**

Year 1997 selected for its representativeness
Values assumed at midnight and interpolated to 3-hourly resolution
Used primarily for comparison and validation


- **Transition from _input.m to final input datasets used in model simulations**

See utility functions in folder /functions
Interpolated to ODE timestep (typically 30 minutes)

---

### 1.3 Upwelling Forcing (Folder physics/Forcing)

- **Variables**

|**Variables**|**Description**|**Units**|
|-|-|-|
|xfil\_input|total upwelling index|m s⁻¹|
|xfilphy\_input|Physically filtered total upwelling index|m s⁻¹|

- **Source**

Debby: total upwelling proxy based on local currents ("A two-dimensional nitrogen and carbon flux model in a coastal upwelling region",  https://doi.org/10.1029/2001GB001451)

- **Notes**

xfilphy\_input modified to impose plateaus during upwelling/downwelling seasons for improved temperature and salinity skill

- **Transition from _input.m to final input datasets used in model simulations**

See utility functions in folder /functions
Interpolated to ODE timestep (typically 30 minutes)

---

### 1.4 Vertical Mixing (Folder physics/Forcing)

- **Variables**

|**Variables**|**Description**|**Units**|
|-|-|-|
|mld\_input|Mixed layer depth|m|
|entrnmnt\_input|Entrainment rate at base of mixed layer|m s⁻¹|

- **Source**

Provided by collaborator Zelalem Engida
Local values for the WCVI shelf and slope
Daily values, 1995-2002

- **Model usage**

Climatological averaging across available years
Interpolated between daily midnight values to a 3-hour resolution
Entrainment recalculated from MLD to ensure physical consistency: ent(t)= (MLD(t+1)−MLD(t−1))/2Δt

- **Transition from _input.m to final input datasets used in model simulations**

See utility functions in folder /functions
Interpolated to ODE timestep (typically 30 minutes)

---

### 1.5 Alongshore Currents (Folder physics/Forcing)

- **Variables**

|**Variables**|**Description**|**Units**|
|-|-|-|
|dc\_input|Davidson Current transport proxy|s⁻¹|
|sbc\_input|Shelf Break Current transport proxy|s⁻¹|
|cu\_input|California Undercurrent transport proxy|s⁻¹|

- **Sources**

"Currents along the pacific Coast of Canada"
"Poleward reach of the California Undercurrent extension"

- **Assumptions**

Currents distributed over fixed widths and depths
Fluxes converted to equivalent volumetric exchange rates (s⁻¹)

- **Transition from _input.m to final input datasets used in model simulations**

See utility functions in folder /functions
Interpolated to ODE timestep (typically 30 minutes)

---

### 1.6 Initial Temperature \& Salinity (Folder physics/IC)

- **Variables**

|**Variables**|**Description**|**Units**|
|-|-|-|
|t\_input|Initial temperature profiles|°C|
|s\_input|Initial salinity profiles|psu|

- **Source**

World Ocean Atlas (WOA)

- **Model usage**

January climatology
Vertical averaging mapped onto WCVI box structure
Shelf and slope treated separately according to model geometry

---

### 1.7 Lateral Boundary Temperature \& Salinity (Folder physics/LB)

- **Variables**

|**Variables**|**Description**|**Units**|
|-|-|-|
|LBtmp\_input|Lateral boundary temperatures|°C|
|LBsal\_input|Lateral boundary salinities|psu|

- **Sources**

WOA climatologies
NARR (rain temperature)
BC stream temperature datasets --> "Stream Temperature Patterns in British Columbia, Canada, based on routine spot measurements"

- **Temporal resolution**

Monthly climatologies interpolated between mid-month values
Converted to sub-daily resolution for model integration

- **Transition from _input.m to final input datasets used in model simulations**

See utility functions in folder /functions
Interpolated to ODE timestep (typically 30 minutes)

---

## 2\. Biological Forcing Datasets (Non-Boundary)

### 2.1 Dissolved Oxygen (Folder bio/Forcing)

- **Variables**

|**Variable**|**Description**|**Units**|
|-|-|-|
|o2\_input|Dissolved oxygen concentration|mol O₂ m⁻³|

- **Source**

World Ocean Atlas (WOA)

- **Model usage**

Converted from µmol kg⁻¹ to mol m⁻³ using seawater density (1026 kg m⁻³)
Shelf and slope treated separately
Annual climatology

O2 is o2\_input formatted as a structure and serves as the data input for model simulations

---

### 2.2 Fisheries Forcing (Folder bio/Forcing)

- **Variables**

|**Variable**|**Description**|**Units**|
|-|-|-|
|fisheriesTS\_input|Fishery mortality and migratory biomass forcing|Mortality rates: s⁻¹; Biomass: t km⁻²|

- **Source**

Ecopath outputs
DFO literature
Survey data and CPUE-based extrapolations

- **Notes**
BIOTS is fisheriesTS\_input formatted as a structure and serves as the data input for model simulations.

---

### 2.3 NEMURO Parameters \& Initial Conditions (Folder bio/IC)

- **Variables**

|**Variables**|**Description**|**Units**|
|-|-|-|
|NemParam|NEMURO parameter set||
|bnem0|Initial NEMURO state variables|mol N. m-3|

- **Source**

NemParam: Calibrated standalone NEMURO model (Chapter 2, PhD thesis)
Use nemuroParamSets.m to generate structure NemParam: NemParam = nemuroParamSets('Kishi et al.A7')
bnem0: see sections 3.4 and 3.5

---

### 2.4 Additional Zooplankton Initial Conditions (Folder bio/IC)

- **Variables**

|**Variable**|**Description**|**Units**|
|-|-|-|
|bz0|Extra-zooplankton initial conditions|mol N. m-3|

- **Source**

see sections 3.4 and 3.5

---

### 2.5 Trophic \& mortality Parameters (Folder bio/Parameters)

- **Variables**

|**Variables**|**Description**|**Units**|
|-|-|-|
|x, d, theta|Foraging arena parameters|no unit|
|grmax|Maximum grazing rates at 0 °C for planktonic groups|s⁻¹|
|thresh|Ivlev threshold value for grazing for planktonic groups|mol N.m-3|
|m0exp|Non-predatory mortality exponent|no unit|

- **Sources**

Foraging arena parameters used for trophic links involving nektonic groups, default values and can be calibrated or adjusted during end-to-end model calibration
grmax: calibrated standalone NEMURO model (Chapter 2, PhD thesis) for NEMURO groups; calibrated for extrazooplankton groups during end-to-end model calibration
thresh: calibrated standalone NEMURO model (Chapter 2, PhD thesis) for NEMURO groups; calibrated for extrazooplankton groups during end-to-end model calibration
m0exp: linear or quadratic mortality - tested during end-to-end model calibration

---

### 2.6 Ecopath model structure (EM)(Folder bio/Parameters)

- **Variables**

|**Variables**|**Description**|**Units**|
|-|-|-|
|EM|MATLAB structure containing the Ecopath representation of the WCVI ecosystem used to parametrize higher trophic levels in the WCVI-E2E model||

- **Description**

It defines:

- living and non-living groups,
- biomasses,
- production and consumption rates,
- diet composition,
- import/export terms,
- trophic linkages

- **Role in the WCVI-E2E model**

- Provides the trophic backbone for nekton and extra zooplankton groups
- Supplies parameters used to compute feeding, mortality, and initial biomass terms
- Ensures consistency with a mass-balanced food-web representation

- **Source**

- Ecopath model developed for the WCVI ecosystem
- Described in Chapter 3 of the PhD thesis
- Generated using the ecopath\_matlab package https://github.com/kakearney/ecopath\_matlab-pkg

- **Temporal resolution**

Static (time-invariant)

- **Spatial resolution**

Implicit (applied uniformly within WCVI model boxes)

- **Notes**

The Ecopath model is used only when isnem = false
When running NEMURO-only configurations, EM is not required

---

### 2.7 Functional Group Classification (types)(Folder bio/Parameters)

- **Variables**

|**Variables**|**Description**|**Units**|
|-|-|-|
|types|cell array defining how each biological state variable in the WCVI-E2E model is treated dynamically.||

- **Description**

Each entry specifies whether a group is:

- a NEMURO planktonic variable,
- an extra zooplankton group,
- or a nekton group.

- **Allowed values**

- NEMURO variable name (e.g. 'ps', 'pl', 'zs', 'zl', 'zp', 'pon')
- 'z' : extra zooplankton group (Ecopath-based, planktonic)
- 'n' : nekton group (Ecopath-based, not physically mixed)

- **Notes**

types must be consistent with the ordering of state variables in the Ecopath model

---

## 3\. Lateral Boundary Biogeochemical Forcing (Folder bio/LB)

This folder contains prescribed lateral boundary concentrations for biogeochemical variables used by the WCVI-E2E model. These forcings represent external nutrient and plankton inputs entering the model domain via open-ocean boundaries, riverine inflows, and coastal currents.
All concentrations are expressed in moles N m⁻³ (or moles Si m⁻³ for silicate) and are interpolated to the model timestep during preprocessing.

### Boundary Geometry and vertical Structure

Boundary concentrations are specified separately for multiple source waters and vertical layers, consistent with model geometry:

|Source|Vertical extent|
|-|-|
|Open ocean (UL)|Surface --> MLD|
|Open ocean (LL)|200-300 m|
|VICC|Surface --> MLD (UL), MLD --> 50 m (LL)|
|DC/SBC|Surface --> 100 m|
|CU|150 -  300 m|
|Juan de Fuca (JdF)|Surface --> MLD|
|River runoff|Surface layer|
|Rain|Surface layer|

---

### 3.1 Inorganic nutrients (NO3, NH4, SiOH4)

Multiple data sources were evaluated. The current implementation uses a hybrid approach, combining World Ocean Atlas (WOA) climatologies with literature-derived estimates for nearshore and boundary currents.

**Open ocean boundaries**

- Source: World Ocean Atlas (WOA)
- Coordinates: ~48.5°N, -126.5°W (offshore edge of slope box)
- Units: originally µmol kg⁻¹, converted to mol m⁻³ using ρ = 1024 kg m⁻³

**Boundary currents**

- DC/CU: single WOA location (~47.5°N, -125.5°W)
- SBC: sparse WOA coverage (selected coastal coordinates)
- VICC: limited WOA coverage, nearest location to the Juan de Fuca Strait

**Literature-derived estimates**

Used where WOA coverage is insufficient or absent:

- NH₄:

UL open ocean ≈ 0.4 µM ("Seasonal variability in nitrogenous nutrition of phytoplankton assemblages in the northeastern subarctic Pacific Ocean")

LL open ocean = UL / 8.5 ("The marine nitrogen cycle: overview and challenges")

DC/SBC assumed same as UL open ocean (no clear longitudinal or seasonal gradient according to literature)

VICC assumed same as UL open ocean (no clear longitudinal or seasonal gradient according to literature)

CU assumed same as LL open ocean

- Rainwater:

NO₃ = 0.037 mol m⁻³ ("Rainwater as a chemical agent of geologic processes")

NH₄ = 0.024 mol m⁻³ ("Rainwater as a chemical agent of geologic processes")

SiOH₄ ≈ 0

- River runoff (http://aquatic.pyr.ec.gc.ca/webdataonlinenational/?lang=en):

Fraser River (Hope station)

Converted from mg L⁻¹ to mol m⁻³

---

### 3.2 Dissolved Organic Nitrogen (DON)

**Source**

"Seasonal changes in the distribution of dissolved organic nitrogen in coastal and open-ocean waters in the Northeast Pacific"

---

### 3.3 Particulate Organic Nitrogen (PON) and Opal

**Source**

"²³⁴Th as a tracer of particulate organic carbon export in the subarctic northeast Pacific Ocean"

**Assumptions**

VICC concentrations ≈ 2.5 × open ocean UL

Opal derived as: Opal = 2 × PON

---

### 3.4 Phytoplankton (Diatoms and Flagellates)

**Source** 

Satellite chlorophyll-a (NSERC project)

**Method**

Partition total Chl-a into diatoms and non-diatoms thanks to proportion of diatoms

Convert: Chl-a → C (C:Chl = 50) → N (C molecular weight and Redfield ratio)

Units converted to mol N m⁻³

**Caveats**

Satellite Chl-a near river mouths and in Juan de Fuca Strait is uncertain

Used with caution in nearshore regions

---

### 3.5 Zooplankton 

**Source**

La Perouse sampling program (data provided by Moira at DFO)

Vertical distribution and migratory behavior taken into account (DVM and diapause)

Night:Day (N:D) ratios applied to correct daytime samples for species undertaking DVM

Unit conversions: mg DW m⁻³ → mol N m⁻³ using calibrated conversion factor 9.4\*10^-6

---

### Transition from _input.m to final input datasets used in model simulations

See utility functions in folder /functions
Interpolated to ODE timestep (typically 30 minutes)
































































