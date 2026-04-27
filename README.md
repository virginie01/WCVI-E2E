# WCVI-E2E: End-to-End Ecosystem Model for the WCVI Coastal Upwelling System

## Overview
This repository contains the WCVI-E2E model, a fully coupled end-to-end (E2E) ecosystem model developed for the West Coast of Vancouver Island (WCVI), Canada.

The model integrates:
- A 2D physical box model
- A NEMURO-based biogeochemical model (lower trophic levels)
- An Ecopath with Ecosim (EwE)-like food-web model (upper trophic levels)

These components are dynamically linked through a two-way coupling approach, allowing feedbacks between physical processes, nutrient dynamics, plankton, and higher trophic levels.

This framework enables the exploration of how long-term environmental change and trophic interactions shape ecosystem structure and function in coastal upwelling systems.

---

## Scientific Background
Understanding how changes propagate across trophic levels is critical for predicting ecosystem responses to long-term climate and anthropogenic changes.

Traditional models often:
- Represent only lower trophic levels (biogeochemical models), or
- Focus on higher trophic levels (food-web models)

The WCVI-E2E model bridges this gap by explicitly coupling these components and resolving two-way feedbacks between them.

This work builds on previous modeling frameworks and extends them to represent a coastal upwelling ecosystem.

---

## Model Description

### Overview of Development
The WCVI-E2E model was developed through a multi-stage process, combining independently developed sub-components into a unified end-to-end framework.

1. **Physical–Biogeochemical Model (NEMURO-based)**  
   A lower trophic level model was first implemented and calibrated independently.  
   - Parameter sensitivity analysis was conducted to identify key drivers  
   - Model parameters were optimized using observational datasets and MATLAB’s surrogate optimization algorithm (based on the Metric Stochastic Response Surface method)    
   - The calibrated model reproduces seasonal and spatial dynamics of nutrients and plankton 

**Notes**

See Section 2. and specifically section 2.2.3. in  `/docs/thesis.pdf` for a description of the calibration/optimization process  

2. **Ecopath-based Food-Web Model**  
   An upper trophic level model was developed to represent the WCVI food web.  
   - Functional groups were defined across trophic levels  
   - Biomass and energy fluxes were balanced using the Ecopath framework  
   - The model provides a static, mass-balanced representation of ecosystem structure

**Notes**

For a complete description of the WCVI Ecopath model, refer to section 3. in `/docs/thesis.pdf`  

3. **Coupled End-to-End Model (WCVI-E2E)**  
   The two components were dynamically coupled based on the EwE approach, using a two-way feedback approach.  
   - Lower trophic levels provide bottom-up forcing to higher trophic levels  
   - Higher trophic levels impose top-down control on plankton dynamics  
   - Alternative coupling strategies were tested and evaluated
   - Heuristic parameter tuning due to model complexity 

**Notes**

Refer to Section 4. of `/docs/thesis.pdf` and `/docs/paper_submitted.docx`.    

---

### Key Features
- Two-way coupling between lower and upper trophic levels  
- Explicit representation of zooplankton as a key interface  
- Integration of physical, chemical, and biological processes  
- Simulation of long-term ecosystem dynamics (multi-decadal) 

---

## Repository Structure  

```
WCVI-E2E/
│
├── data/        # Input datasets and documentation on data collection, cleaning, and formatting
├── src/         # Core model implementation
├── scripts/     # Simulation, analysis, and plotting scripts
├── results/     # Model outputs (including long-term simulations)
├── external/    # External dependencies (e.g., Ecopath MATLAB functions)
├── docs/        # Thesis and manuscript
└── README.md
```

---

## How to Run the Model

### Requirements
This software requires [MATLAB](https://www.mathworks.com/products/matlab.html). It has been run using R2024b.

The full package requires additional Mathworks toolboxes:

- [Statistics and Machine Learning Toolbox](https://www.mathworks.com/products/statistics.html)
- [Image Processing Toolbox](https://www.mathworks.com/products/image-processing.html)

### Downloading the Repository

There are two ways to obtain a local copy of this repository:

**Option 1: Using Git (recommended)**  
git clone https://github.com/virginie01/WCVI-E2E.git

**Option 2: Download as ZIP**  
- Click the green **Code** button at the top of the repository page on GitHub
- Select **Download ZIP**
- Extract the downloaded file to your desired location 

### Running the Model

After downloading the repository:  
1. Open MATLAB
2. Navigate to the project folder
3. Open `scripts/WCVI-E2E_run.m`
4. Follow instructions to add all the necessary subfolders to your Matlab path
5. Follow instructions to run either the physics-only version or the full end-to-end version

### Syntax

- To run the physics-only version:  
`physicalmodel('outputs_physics');`

- To run the full E2E version:  
`physicalmodel('outputs_e2e','biofun',@biomodel);`

---

## Outputs

Simulation outputs are saved in the `/results/` directory and the appropriate subfolder and include a table `archivedata.mat` containing:
- Monthly time series of ecosystem state variables
- Monthly time-series of diagnostic variables (derivatives, intermediate fluxes, additional diagnostic variables)

Refer to `scripts/WCVI-E2E_run.m` and `scripts/WCVIE2E_variables.docx` for plotting simulation outputs

---

## Results

The `/results/` folder contains long-term simulations used in PhD analyses along with a README document

---

## Related Work

**PhD thesis**  
A full description of the model and its development is available in:
- `/docs/thesis.pdf`

**Manuscript**  
A condensed version focusing on coupling strategies and ecosystem dynamics:  
- `/docs/paper_submitted.docx`

---

## Pre-existing Work and References

1. The 2D physical box model builds upon and extends the following framework:  
- Ianson D, Allen SE. A two‐dimensional nitrogen and carbon flux model in a coastal upwelling region. Global Biogeochemical Cycles 2002;16. https://doi.org/10.1029/2001GB001451.

 Key extensions include:
- Reprogrammed in MATLAB (initially coded in FORTRAN)
- Addtitional alongshore currents
- Temperature coded as a prognostic variable

2. The biogechemical model builds upon and extends the NEMURO framework:  
- Kishi MJ, Kashiwai M, Ware DM, Megrey BA, Eslinger DL, Werner FE, et al. NEMURO—a lower trophic level model for the North Pacific marine ecosystem. Ecological Modelling 2007;202:12–25. https://doi.org/10.1016/j.ecolmodel.2006.08.021.

 Key modifications of the NEMURO model include:
- Modifications of the nitrification and decomposition formulations (see `/docs/paper_submitted.docx`)

3. The E2E framework and coupling strategies draw from the WCE model:  
- Kearney KA, Stock C, Aydin K, Sarmiento JL. Coupling planktonic ecosystem and fisheries food web models for a pelagic ecosystem: Description and validation for the subarctic Pacific. Ecological Modelling 2012;237–238:43–62. https://doi.org/10.1016/j.ecolmodel.2012.04.006.

 Key adaptations include:
- Model refactoring for coastal upwelling systems
- Testing of one-way forcing vs. two-way coupling on ecosystem structure and function

---

## Notes and Limitations

- This codebase is research-oriented and not optimized for general-purpose use
- The model is tailored to the WCVI ecosystem
- Large simulations (e.g., 100-year runs with monthly output) require ~48 hours on a standard local machine, reflecting a balanced trade-off between model granularity and computational efficiency

---

## Author

Virginie C. Bornarel
PhD, Institute for the Oceans and Fisheries
University of British COlumbia

---

## Acknowledgements
This model builds on prior work and code contributions from collaborators and previous studies in physical, biogeochemical, and ecosystem modeling.

This work also builds on contributions from multiple collaborators and data providers.

The author gratefully acknowledges:

- Debby Ianson (Fisheries and Oceans Canada, Institute of Ocean Sciences) for her support and guidance in developing the physical–biogeochemical component of the model  
- Fisheries and Oceans Canada (DFO) staff for providing key datasets used in model development, including zooplankton, groundfish, sardine, herring, and salmon data  
- Specific contributions from:
  - Moira Galbraith (zooplankton data, LaPerouse Zooplankton Monitoring Program)  
  - Maria Surry (groundfish data)  
  - Linnea Flostrand (Pacific sardine data)  
  - Matthew Grinnell (herring assessment outputs)  
  - John Davidson (salmon data)  
  - Andrew Edwards (rockfish assessment outputs)  
- Sarah Z. Rosengard (University of British Columbia) for satellite-based phytoplankton composition data  

These contributions were essential for model development, calibration, and validation.

---

## Citation
If you use this model, please cite:

Bornarel, V. C. (2025). *A newly-developed end-to-end model for coastal upwelling systems: Insights into the mechanistic linkages between trophic levels under varying environmental conditions*. PhD Thesis, University of British Columbia.

---

## Usage Note
This code is provided for research and educational purposes.  
Please cite the associated thesis or publication when using this model.

---

## License
This project is licensed under the MIT License. See the LICENSE file for details.




