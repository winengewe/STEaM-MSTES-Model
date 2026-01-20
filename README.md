# STEaM-MSTES-Model: Mine Shaft Thermal Energy Storage Simulation

**Techno-economic modeling of GigaWatt-hour Subsurface Thermal Energy Storage in engineered structures and legacy mine shafts.**

This repository contains Python models developed by the Energy System Research Unit (ESRU) at the University of Strathclyde as part of the EPSRC grant **EP/W027763/1 (STEaM project)**. 

The model investigates the feasibility and potential for Mine Shaft Thermal Energy Storage (MSTES) to deliver low-cost renewable district heating and support the balancing of future renewable electricity grids.

## 📂 Repository Structure

The project consists of three main Python files:

### 1. `shaftstore_1d_0i.py` (The Physics Core)
This is the core class file that models the thermodynamics of the mine shaft. It utilizes a 1D finite volume/difference approach to simulate:
* **Thermal Nodes:** Calculates heat transfer between water layers, the shaft wall, insulation, and surrounding ground rings.
* **Stratification:** Models buoyancy and mixing within the water column.
* **Heat Loss:** Calculates radial and vertical heat diffusion into the surrounding geology.

### 2. `MSTES-HP.py` (Heat Pump Scenario)
A control and economic wrapper that simulates a **Heat Pump (HP)** integrated system.
* **Logic:** Prioritizes using wind surplus (low tariff) electricity to charge the store or supply demand.
* **Simulation:** Runs year-long simulations (hourly or half-hourly) over a defined lifetime (e.g., 20 years).
* **Outputs:** Calculates Levelized Cost of Heat (LCOH), Coefficient of Performance (COP), and system efficiency.

### 3. `MSTES-CHP.py` (Combined Heat & Power Scenario)
A control and economic wrapper that simulates a **Combined Heat and Power (CHP)** integrated system.
* **Logic:** Models CHP operation modes (power-led vs. heat-led) and integrates Gas Boilers (GB) for peak loads.
* **Economics:** Accounts for fuel costs, electricity export income, and grid tariffs.

## ⚙️ Features

* **Grid Integration:** Models interaction with grid constraints, specifically targeting wind surplus utilization.
* **Sensitivity Analysis:** Capable of looping through various parameters (Store Height, Tariffs, HP/CHP sizing, Insulation thickness).
* **Ground Heating Envelope:** Tracks the long-term temperature evolution of the ground surrounding the shaft.
* **Automated Reporting:** Generates detailed CSV logs (Key Performance Indicators, Time-series results) and visual plots (Temperature profiles, Energy balance graphs).

## 🚀 Getting Started

### Prerequisites
The code requires Python 3.x and the following libraries:
* `numpy`
* `pandas`
* `matplotlib`
* `scipy`
* `tqdm` (for progress bars)
* `IPython`

### Installation
1. **Clone the repository:**
   ```bash
   git clone [https://github.com/winengewe/STEaM-MSTES-Model.git](https://github.com/winengewe/STEaM-MSTES-Model.git)
   ```
2.  **Install dependencies:**
    ```bash
    pip install numpy pandas matplotlib scipy tqdm ipython
    ```
### Configuration & Input Data
Note: The scripts currently point to local file paths (e.g., C:\Users\...). You must update the data_loc variable in MSTES-HP.py and MSTES-CHP.py to point to your local data directory.
The models expect the following CSV input files in your data directory:
gsp_data_1.csv: Electricity demand, Wind generation, and Heat demand profiles.
geo_win.csv: Thermal conductivity data for ground, concrete, and insulation.
geo_win2.csv: Heat capacity, density, and porosity data for the ground.

### Running the Model
Run the desired scenario script directly:
   ```bash
   python MSTES-HP.py
    # OR
    python MSTES-CHP.py
   ```
### Outputs
The simulation automatically creates results folders based on the simulation parameters. Outputs include:
* **kpi.csv:** Lifetime system performance metrics (LCOH, COP, Heat Loss %).
* **RES.csv / tRES.csv:** Time-series energy balance data.
* **f_n_t.csv:** Full nodal temperature history.
* **Plots:** PNG images showing node temperatures, seasonal performance (Summer/Winter), and ground heating envelopes.

## Authors & Acknowledgments
* **Authors:** Ewe Win Eng, Graeme Flett, Paul Tuohy
* **Affiliation:** Energy System Research Unit (ESRU), Dept. of Mechanical and Aerospace Engineering, University of Strathclyde, Glasgow, UK.
* **Funding:** Engineering and Physical Sciences Research Council (EPSRC) grant EP/W027763/1.
   
