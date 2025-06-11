"""
MSTES-CHP.py (Mine Shaft Thermal Energy Storage - Combined Heat & Power)
Model intro: EPSRC grant EP/W027763/1 STEaM project
GigaWatt-Hour Subsurface Thermal Energy storAge: Engineered structures and legacy Mine shafts
A techno-economic model to investigates the feasibility and potential for MSTES 
to deliver low-cost renewable district heating 
and support balancing of the future renewable electricity grid  
by Ewe Win Eng*, Graeme Flett, Paul Tuohy
Energy System Research Unit (ESRU)
Department of Mechanical and Aerospace Engineering
University of Strathclyde, Glasgow G1 1XJ, UK
*engwinewe@gmail.com 
"""

"Import python libraries and packages"
from IPython import get_ipython
ipython = get_ipython()
if ipython is not None:  # Ensure the script runs inside IPython
    ipython.run_line_magic('reset', '-sf')  # Reset all variables
import os
import sys
def stop(error_message):
    """
    Terminates the script with a non-zero exit code and prints an error message.
    """
    print(f"{error_message}")
    sys.exit(1)
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from itertools import product
from tqdm import trange # see looping progress bar %
import datetime
current_time = datetime.datetime.now().strftime("%Y%m%d_%H%M%S") # get current time
import shaftstore_1d_0i
DegC = '°C'
GBP = '£'
combine_tRESnkpi = 0 # [1] combine all tRES.csv AND kpi.csv files

#%% Inputs
ts       = 3 # Timesteps per hourly base data
Mod_Year = 10 # Run years
HOURS_PER_YEAR = 8760  # 1 data point/hour * 24 hours/day * 365 days/year
nt = int(HOURS_PER_YEAR * ts * Mod_Year)
print(f"nt (Total Time Steps): {nt}") 

data_loc = r'C:\Users\cmb22235\OneDrive - University of Strathclyde\Desktop\STEaM WP4 team\MSTES-insul' # folder path
Dat = np.loadtxt(f'{data_loc}\\gsp_data_1.csv', delimiter=',', skiprows=1) # file path
Dat = Dat[:, 1:] # Skip the first column
# [0] Current Coylton GSP Electricity Demand (MWh)
# [1] Future Coylton Electricity Demand (MWh)
# [2] Coylton 33kV Wind (MWh)
# [3] Coylton GSP Wind (MWh)
# [4] Domestic Heat per kWh
# [5] Air Temp (°C)
# [6] West Whitlawburn District Heating (kWh)

Temp = Dat[:,5] # Temperature data equals to Dat file column 6th (python begins at 0)
Dat_row = len(Dat[:,1]) # check the import yearly data row
def validate_Dat(Dat, Dat_row,):
    """
    Validates the Yearly Data Input File.
    Terminates the script if the data length is invalid.
    """
    if Dat_row not in (8760, 17520):
        stop (f"Invalid Dat Row: {Dat_row}. Expected 8760 (hourly) or 17520 (half-hourly) Yearly Data.")
    if Dat_row == 17520:
        Dat = Dat.reshape(8760, 2, len(Dat[1,:])).sum(axis=1) # Sum of every two rows
    return Dat
if HOURS_PER_YEAR == 8760: # if want to use 8760 yearly hours
    Dat = validate_Dat(Dat, Dat_row)  
    MTS = int(3600 / ts) # Timestep in seconds (needs to be run at hourly)
else:
    MTS = int(1800 / ts) # Timestep in seconds (needs to be run at half-hourly)
print(f"MTS (Model Time Step): {MTS} seconds")
Dat = np.repeat(Dat,ts,0) / ts # Set demand and surplus data for sub-half-hour timesteps
Temp = np.repeat(Temp,ts,0) # Ambient temperature for ASHP COP
HD_Option = 1 # [1]=WestWhilatwburn Yearly HD; [2]=Timestep HD

#%%% Store
RLW_A         = [38] # Ground and Water layers height (m) # [0,8,18,28,38,48] = NoStore case [28] = match old case 280m
number_nodes  = 15  # nodes for 8 ground rings, shaftwall, insulation and air/minewater
number_layers = 14  # put [56] if create 52 water layers, must add 4 air layers (Monktonhall basis)
top_wat       = 5   # top layer of heated water section (3 dummy air layers not included in the thermal analysis)
RLA           = 20. # Consolidated air Layer (m) (modelled as a single 20m zone)
mu            = 50. # Buoyancy model factor (set lower (20-50) if 'overflow encountered in exp' warnings occur)
r             = 3.5 # radius of TS (m)
# ithk = 0.9 # insulation thickness (m)
ithk = 0.1 # insulation thickness (m)
print(f"nodes={number_nodes} layer={number_layers}, radius={r}m")
icp  = 880 # insulation heat capacity (J/kgK)
iden = 1600 # insulation density (kg/m3)
cden = iden # concrete wall density (kg/m)
Rx = np.array([256.+ithk, # Ground ring 12 
               128.+ithk, # Ground ring 11
               64.+ithk, # Ground ring 10
               48.+ithk, # Ground ring 9
               32.+ithk, # Ground ring 8
               24.+ithk, # Ground ring 7
               20.+ithk, # Ground ring 6
               16.+ithk, # Ground ring 5
               12.+ithk, # Ground ring 4
               8.+ithk, # Ground ring 3
               4.+ithk, # Ground ring 2
               2.+ithk, # Ground ring 1
               1.+ithk, # Mineshaft wall
               ithk+ithk, # Insulation layer
               ithk, # Air/Minewater 
               -r + ithk]) \
               + (r - ithk) # Node outer radii (m); related to the number_nodes
Rx_len = len(Rx)
if Rx_len != number_nodes+1: # if Rx_len not equal to number of nodes + 1
    stop(f"Invalid Rx_len: {Rx_len}. Expected {number_nodes+1}")
Qlimit = 0 # charge limit, 0 <= Qlimit < 1, 0=use store if it has charge, 0.5=use store if store has half-fully charge [0,0.5]

#%%% Ground properties
# XTC = np.loadtxt(f'{data_loc}\\mth_geo_1v2.csv',delimiter=',') 
# Geo = np.loadtxt(f'{data_loc}\\mth_geo2_1.csv',delimiter=',') 
XTC = np.loadtxt(f'{data_loc}\\geo_win.csv',delimiter=',', skiprows=1) # ground, concrete, insulation, water, air thermal conductivity (W/mK)
XTC = XTC[:, 1:]  # Skip the first column
Geo = np.loadtxt(f'{data_loc}\\geo_win2.csv',delimiter=',', skiprows=1) # ground heat capacity (J/kgK), density (kg/m3), porosity (fraction)
Geo = Geo[:, 1:]  # Skip the first column
# Validate
XTC_row = len(XTC[1,:]) # check XTC row
XTC_col = len(XTC[:,1]) # check XTC column
Geo_col = len(Geo[:,1]) # check Geo column
if XTC_row != number_nodes: # if XTC row not equal to number of nodes
    stop(f"Invalid XTC Row: {XTC_row}. Expected {number_nodes}. Check geo_win.csv and number_nodes")
if XTC_col and Geo_col < number_layers:
    stop(f"Invalid XTC {XTC_col} and Geo {Geo_col} Column. Expected {number_layers}. Check geo_win.csv and geo_win2.csv and number_layers")

#%%% Cost
surp_tarA = [0.1] # tariff during wind surplus period (£/kWh e)
grid_tarA = [0.3] # tariff during wind shortfall period (£/kWh e) [0.1,0.3,0.6,1.2] base0.3
lifetime = 20 # lifetime of system (years) for LCOH
CHP_CAPEX = 4000 # CAPEX of CHP (£/kW) 3-5.5k
GB_CAPEX = 100 # CAPEX of Gas Boiler (£/kW) 1.5-3.4k
DR       = 0.05 # Discount Rate
Inf      = 0.02 # Inflation
chp_fuel_costA = [0.075] # (£/kWh) biomass [0.075] = biomass; hydrogen = [0.375] base0.075
gb_fuel_costA = chp_fuel_costA # (£/kWh) # assume gas boiler use same fuel with CHP plant

#%%% Sensitivity Inputs
HDFacA    = [1e7] # Annual Heat Demand (kWh) for HD_Option = [2]
Size_CHPA = [4200] # Main CHP Capacity (kW) - CHP size assumed to be > peak heat demand (4088.44 WWDH*3 peak HD)
Size_GBA  = [3100] # GB Heat Capacity (kW)
disDTA    = [10] # Store discharge DeltaT (C) # fixed value!!!!
OptionA   = [1] # [1] = CHP only (CHP always on to max export income); [2] = CHP + GB (CHP off during wind surplus)
CHPmodeA  = [1] # [1] = Run to its peak every timestep; [2] = Stop if HD satisfied
eff_GB_thermal = 0.88 # Gas boiler efficiency
eff_CHP_choice = [1] # select 'eff_CHP' [1,2,3,4] base [1]
eff_CHP_map = { # CHP output efficiency for elec,heat
1: (0.38, 0.45),  # gas turbine
2: (0.25, 0.58),  # steam / gas turbine
3: (0.13, 0.70),  # steam
4: (0.50, 0.33),  # fuel cell
}     
#%%% Temperature Limits (°C)
initial_node_temp = 12 # Initial node temperature
heat_tempA         = [50] # Heat demand minimum supply temperature [50]
min_tempA          = [50] # Minimum store supply temperature (must <= heat_temp & store_temp, = heat_temp for ASHP cases) [50]
store_tempA        = [55] # Temperature to fill store (no direct flow from store if < heat_temp) [55,70]

if min_tempA > heat_tempA or min_tempA > store_tempA:
    stop(f"Invalid min_tempA: {int(min_tempA)}{DegC}. Expected less than heat_tempA {int(heat_tempA)}{DegC} and store_tempA {int(store_tempA)}{DegC}")
if OptionA == 1 or OptionA == 3 or OptionA == 4:
    if min_tempA != heat_tempA:
        stop(f"Invalid min_tempA:{int(min_tempA)}{DegC}. Expected equal heat_tempA {int(heat_tempA)}{DegC} for ASHP cases.")

#%% Looping
# Flatten the nested loops using itertools.product
for RLW,surp_tar,grid_tar,chp_fuel_cost,gb_fuel_cost,HDFac,Size_CHP,Size_GB,disDT,Option,CHPmode,eff_CHP,heat_temp,min_temp,store_temp in product(
        RLW_A,surp_tarA,grid_tarA,chp_fuel_costA,gb_fuel_costA,HDFacA,Size_CHPA,Size_GBA,disDTA,OptionA,CHPmodeA,eff_CHP_choice,heat_tempA,min_tempA,store_tempA):
    eff_CHP_electrical, eff_CHP_thermal = eff_CHP_map.get(eff_CHP, (None, None)) # read 'eff_CHP_map' string
    h = RLA + (number_layers - (top_wat-1)) * RLW   # TS height (m) 
    TS_volume = np.pi * r **2 * h # Thermal store volume (m3)
    Qmax = (np.pi * r **2 * RLW  * disDT * 4.181 * 1000. * (number_layers - top_wat + 1)) / 3600. # Max store heat charge (fixed delta T)
    if RLW == 0:
        h = 0  # NoStore case
        TS_CAPEX = 0 # Store Capital Cost (£/m3)
        print("NoStore Case")
    else:
        TS_CAPEX = 7982. * TS_volume**-0.483 # Store Capital Cost (£/m3)

    # Data Initialisation
    eWdS = np.zeros((nt,1)) # Available low tariff/wind surplus electricity 
    DemH = np.zeros((nt,1)) # Heat Demand 
    chp2HD = np.zeros((nt,1)) # heat from CHP to HD  
    chp2tes = np.zeros((nt,1)) # heat from CHP to store 
    hS2H = np.zeros((nt,1)) # Store to Heat Demand    
    gb2HD = np.zeros((nt,1)) # heat from GB to HD    
    waste_heat = np.zeros((nt,1)) # chp waste heat
    Un_H = np.zeros((nt,1)) # Unmet heat demand
    CHPrem = np.zeros((nt,1)) # Residual CHP capacity at end of timestep
    GBrem = np.zeros((nt,1)) # Residual GB capacity at end of timestep
    chp_fuel2H = np.zeros((nt,1)) # chp fuel to HD 
    chp_fuel2S = np.zeros((nt,1)) # chp fuel to store
    chp_fuel_used = np.zeros((nt,1)) # chp fuel used 
    gb_fuel2H = np.zeros((nt,1)) # gb fuel to HD
    gb_fuel2S = np.zeros((nt,1)) # gb fuel to store
    gb_fuel_used = np.zeros((nt,1)) # gas boiler fuel used 
    chp2g = np.zeros((nt,1)) # CHP exports electricity to grid
    chp_export_income = np.zeros((nt,1)) # chp export income (£)
    HeatOpex = np.zeros((nt,1)) # Operating (electricity) costs
    THL = np.zeros((nt,1)) # Store heat losses
    
    Qavail = np.zeros((nt+1,1)) # Store available charge
    OpH = np.zeros((lifetime)) # Yearly operating costs
    
    # Nodes Temp Initialisation
    nodes_temp = np.ones(number_layers * number_nodes) * initial_node_temp # Initial node temperatures
    f_n_t = np.zeros((int(nt/ts), number_layers * number_nodes)) # Results setup
    tr = 0
    Qp = np.zeros(number_layers)
    
    #%% Basic Store Analysis (all kWh)
    for t in trange(nt,desc='Timestep'): # timestep
        td = np.remainder(t,HOURS_PER_YEAR*ts)
        
        Av_WS = max(0.,Dat[td,2] - Dat[td,0] - 0) * 1000. # Available Wind Surplus (kWh) (Coylton 33kV Wind(2) - Current Elec Demand(0))
        eWdS[t,0] = np.copy(Av_WS) # Available Wind Surplus
        
        if HD_Option == 1: 
            DemH[t,0] = Dat[td,6] * 3 # West Whitlawburn District Heating x3, ~10GWh (kWh th)
        elif HD_Option == 2:
            DemH[t,0] = Dat[td,4] * HDFac # Timestep Heat Demand (kWh)
        
        Res_CHP = (Size_CHP / (3600. / MTS)) # Timestep CHP output available (kWh th)
        Res_GB = (Size_GB / (3600. / MTS)) # Timestep GB output available (kWh th)
        Res_HD = np.copy(DemH[t,0]) # Timestep heat demand (kWh th)
        
        #%% CHP only case
        if Option == 1: # CHP on during wind surplus/low tariff
            Size_GB = 0; Res_GB = 0
            
            #%%% Wind shortfall/high tariff
            if Av_WS == 0: # if wind shortfall/high tariff (CHP on)
            
                #%%%% Priority 1: CHP supply heat to HD
                if Res_CHP > 0 and Res_HD > 0: # if residual CHP and HD exists
                    chp2HD[t,0] = min(Res_CHP, Res_HD) # heat from CHP to HD (kWh th)
                    Res_HD -= chp2HD[t,0] # Residual HD (kWh th)
                    Res_CHP -= chp2HD[t,0] # Remaining CHP heatout (kWh th)
                    chp_fuel2H[t,0] += chp2HD[t,0] / eff_CHP_thermal # chp fuel to HD (kWh)
                    
                #%%%% Priority 2: Store supply heat to HD
                if h > 0 and Res_HD > 0: # if store and residual HD exists
                    if Qavail[t,0] > Qlimit * Qmax: # Use store based on store charge
                        # if top Store > HD Temp & bottom > min_temp - 10 (prevents return temperature >> bottom temp, prevent mixing, stratification)
                        if nodes_temp[top_wat*number_nodes-1] > heat_temp and nodes_temp[number_layers*number_nodes-1] > (min_temp - disDT): 
                            # Store direct supply to HD
                            hS2H[t,0] += np.copy(Res_HD) # Heatout direct from Store to HD
                            Res_HD -= hS2H[t,0] # Residual Heat Demand
                        
                #%%%% Priority 3: CHP supply heat to Store (Charge)  
                if h > 0 and Res_CHP > 0.: # if store and residual CHP exists
                    if nodes_temp[number_layers*number_nodes-1] < store_temp - 2.: # Only charge if bottom of Store is less than a limit (prevents > node volume flow per timestep)
                        chp2tes[t,0] = Res_CHP # heat from CHP to store (kWh th)
                        Res_CHP -= chp2tes[t,0] # Remaining CHP heatout (kWh th)
                        chp_fuel2S[t,0] += chp2tes[t,0] / eff_CHP_thermal # chp fuel to Store (kWh)
                            
            #%%% Wind surplus/low tariff    
            elif Av_WS > 0: # if wind surplus/low tariff (CHP on)
            
                #%%%% Priority 1: Store supply heat to HD
                if h > 0 and Res_HD > 0: # if store and residual HD exists
                    if Qavail[t,0] > Qlimit * Qmax: # Use store based on store charge
                        # if top Store > min_temp & bottom > min_temp - disDT (prevents return temperature >> bottom temp, prevent mixing, stratification)
                        if nodes_temp[top_wat*number_nodes-1] > min_temp and nodes_temp[number_layers*number_nodes-1] > (min_temp - disDT): 
                            # Store direct supply to HD
                            if nodes_temp[top_wat*number_nodes-1] > heat_temp: # top Store temp > HD Temp
                                hS2H[t,0] += np.copy(Res_HD) # Heatout direct from Store to HD
                                Res_HD -= hS2H[t,0] # Residual Heat Demand
                
                            #%%%% Priority 2: Store + CHP supply heat to HD
                            elif Res_CHP > 0 and (heat_temp - disDT) < nodes_temp[top_wat*number_nodes-1] < heat_temp: # if Residual CHP and HD exists and return HD < top store < supply HD temp      
                                hS2H[t,0] = (nodes_temp[top_wat*number_nodes-1] - (heat_temp - disDT)) / (heat_temp - (heat_temp - disDT)) * Res_HD # heat from store to HD (kWh th)
                                chp2HD[t,0] = min(Res_CHP, Res_HD - hS2H[t,0]) # CHP boost heat from store to HD (kWh th)                                      
                                Res_HD -= (hS2H[t,0] + chp2HD[t,0]) # Residual HD (kWh th)
                                Res_CHP -= chp2HD[t,0]
                                chp_fuel2S[t,0] += chp2HD[t,0] / eff_CHP_thermal # chp fuel to Store (kWh)
                                
                #%%%% Priority 3: CHP supply heat to HD
                if Res_CHP > 0 and Res_HD > 0 : # if Residual CHP and HD exists
                    chp2HD[t,0] = min(Res_CHP, Res_HD) # heat from CHP to HD (kWh th)
                    Res_HD -= chp2HD[t,0] # Residual HD (kWh th)
                    Res_CHP -= chp2HD[t,0] # Remaining CHP Qout (kWh th) 
                    chp_fuel2H[t,0] += chp2HD[t,0] / eff_CHP_thermal # chp fuel to HD (kWh)

        #%% CHP + post GB case 
        elif Option == 2: # CHP off during wind surplus/low tariff
            
            #%%% Wind shortfall/high tariff
            if Av_WS == 0: # if wind shortfall/high tariff (CHP on)
            
                #%%%% Priority 1: CHP supply heat to HD
                if Res_CHP > 0 and Res_HD > 0: # if Residual CHP and HD exists
                    chp2HD[t,0] = min(Res_CHP, Res_HD) # heat from CHP to HD (kWh th)
                    Res_HD -= chp2HD[t,0] # Residual HD (kWh th)
                    Res_CHP -= chp2HD[t,0] # Remaining CHP heatout (kWh th)
                    chp_fuel2H[t,0] += chp2HD[t,0] / eff_CHP_thermal # chp fuel to HD (kWh)
                    
                #%%%% Priority 2: Store supply heat to HD
                if h > 0 and Res_HD > 0: # if store and residual HD exists
                    if Qavail[t,0] > 0 * Qmax: # Use store based on store charge (charge limit can be varied, set to > 0 to ignore)
                        # if top Store > min_temp & bottom > min_temp - disDT (prevents return temperature >> bottom temp, prevent mixing, stratification)
                        if nodes_temp[top_wat*number_nodes-1] > min_temp and nodes_temp[number_layers*number_nodes-1] > (min_temp - disDT): 
                            # Store direct supply to HD
                            if nodes_temp[top_wat*number_nodes-1] > heat_temp: # top Store temp > HD Temp
                                hS2H[t,0] += np.copy(Res_HD) # Heatout direct from Store to HD
                                Res_HD -= hS2H[t,0] # Residual Heat Demand
                    
                            #%%%% Priority 3: Store + GB supply heat to HD
                            elif Res_GB > 0 and (heat_temp - disDT) < nodes_temp[top_wat*number_nodes-1] < heat_temp: # if Residual GB and HD exists and return HD < top store < supply HD temp      
                                hS2H[t,0] = (nodes_temp[top_wat*number_nodes-1] - (heat_temp - disDT)) / (heat_temp - (heat_temp - disDT)) * Res_HD # heat from store to HD (kWh th)
                                gb2HD[t,0] = min(Res_GB, Res_HD - hS2H[t,0]) # GB boost heat from store to HD (kWh th)                                     
                                Res_HD -= (hS2H[t,0] + gb2HD[t,0]) # Residual HD (kWh th)
                                Res_GB -= gb2HD[t,0]
                                gb_fuel2S[t,0] += gb2HD[t,0] / eff_GB_thermal # GB fuel to Store (kWh)
                                    
                #%%%% Priority 4: CHP supply heat to Store (charge)
                if h > 0 and Res_CHP > 0.: # if store and residual CHP exists
                    if nodes_temp[number_layers*number_nodes-1] < store_temp - 2.: # Only charge if bottom of Store is less than a limit (prevents > node volume flow per timestep)
                        chp2tes[t,0] = Res_CHP # heat from CHP to store (kWh th)
                        Res_CHP -= chp2tes[t,0] # Remaining CHP heatout (kWh th)
                        chp_fuel2S[t,0] += chp2tes[t,0] / eff_CHP_thermal # chp fuel to Store (kWh)
                            
            #%%% Wind surplus/low tariff     
            elif Av_WS > 0: # if wind surplus (CHP off)
            
                #%%%% Priority 1: Store supply heat to HD
                if h > 0 and Res_HD > 0: # if store and residual HD exists
                    if Qavail[t,0] > Qlimit * Qmax: # Use store based on store charge
                        # if top Store > HD Temp & bottom > min_temp - 10 (prevents return temperature >> bottom temp, prevent mixing, stratification)
                        if nodes_temp[top_wat*number_nodes-1] > min_temp and nodes_temp[number_layers*number_nodes-1] > (min_temp - disDT): 
                            # Store direct supply to HD
                            if nodes_temp[top_wat*number_nodes-1] > heat_temp: # top Store temp > HD Temp
                                hS2H[t,0] = Res_HD # heat from store to HD (kWh th)
                                Res_HD -= hS2H[t,0] # Residual HD (kWh th)
                            
                            #%%%% Priority 2: Store + GB supply heat to HD
                            elif Res_GB > 0 and (heat_temp - disDT) < nodes_temp[top_wat*number_nodes-1] < heat_temp: # if Residual GB exists and return HD < top store < supply HD temp      
                                hS2H[t,0] = (nodes_temp[top_wat*number_nodes-1] - (heat_temp - disDT)) / (heat_temp - (heat_temp - disDT)) * Res_HD # heat from store to HD (kWh th)
                                gb2HD[t,0] = min(Res_GB, Res_HD - hS2H[t,0]) # GB boost heat from store to HD (kWh th)                                      
                                Res_HD -= (hS2H[t,0] + gb2HD[t,0]) # Residual HD (kWh th)
                                Res_GB -= gb2HD[t,0]
                                gb_fuel2S[t,0] += gb2HD[t,0] / eff_GB_thermal # GB fuel to Store (kWh)
                
                #%%%% Priority 3: GB supply heat to HD
                if Res_GB > 0 and Res_HD > 0 : # if Residual GB and HD exists
                    gb2HD[t,0] = min(Res_GB, Res_HD) # heat from GB to HD (kWh th)
                    Res_HD -= gb2HD[t,0] # Residual HD (kWh th)
                    Res_GB -= gb2HD[t,0] # Remaining GB Qout (kWh th) 
                    gb_fuel2H[t,0] += gb2HD[t,0] / eff_GB_thermal # GB fuel to HD (kWh)
        
        chp_fuel_used[t,0] = chp_fuel2H[t,0] + chp_fuel2S[t,0] # total chp fuel used (kWh)
        gb_fuel_used[t,0] = gb_fuel2H[t,0] + gb_fuel2S[t,0] # total gb fuel used (kWh)
        
        if CHPmode == 1 and Res_CHP > 0: # Power mode and residual CHP heat exists    
            chp_fuel_used[t,0] += Res_CHP / eff_CHP_thermal # Runs at 100% in power priority mode, more waste heat
            waste_heat[t,0] = Res_CHP # chp waste heat (kWh th) to optimize CHP size for heat priority
                                                                                                     
        Un_H[t,0] += np.copy(Res_HD) # Unmet heat demand (kWh th)  
        if Un_H[t,0] > 1e-10:
            print(t, Un_H[t,0]) # check unmet heat demand                                
    
        #%%% Other Priorities 
        CHPrem[t,0] = np.copy(Res_CHP) # Residual CHP Capacity (optimize size for heat priority)
        GBrem[t,0] = np.copy(Res_GB) # Residual GB Capacity (reduce or increase size for heat priority)
              
        chp2g[t,0] = eff_CHP_electrical * chp_fuel_used[t,0] # chp export electricity to grid (kWh e)
        chp_export_income[t,0] = grid_tar * chp2g[t,0] # chp export income (£)        
        HeatOpex[t,0] = (chp_fuel_used[t,0] * chp_fuel_cost + gb_fuel_used[t,0] * gb_fuel_cost) - chp_export_income[t,0]  # Heat Opex
        
        if h > 0:
            if chp2tes[t,0] - hS2H[t,0] > 0.: # Charge > Discharge
                charge = chp2tes[t,0] - hS2H[t,0] # Net Heat Added to Store
                discharge = 0.
            elif chp2tes[t,0] - hS2H[t,0] < 0.: # Discharge > Charge
                charge = 0.
                discharge = -(chp2tes[t,0] - hS2H[t,0]) # Net Heat Removed from Store
            else:
                charge = 0.
                discharge = 0.
                
            return_temp = nodes_temp[top_wat*number_nodes-1] - disDT # Return temperature from Heat Demand loop (top Store temp minus deltaT)
            # print(discharge, return_temp)
            next_nodes_temp, Hloss, Qp, Cxy, hms = shaftstore_1d_0i.ShaftStore(Rx, number_layers, number_nodes, top_wat, RLA, RLW, XTC, Geo, mu, ithk, icp, iden, cden).new_nodes_temp(nodes_temp, store_temp, return_temp, charge, discharge, disDT, Qp, MTS) # Calculate new store temperatures at end of timestep
            if np.isnan(nodes_temp).any(): # check NaN values generate due to ts too small
                stop("ts too small. Increase the ts")
            nodes_temp = next_nodes_temp[1]

            THL[t,0] = Hloss # Heat loss from Store to ground
            
            # Available heat from store (estimate based on all layers above min_temp)
            fac = 0
            cutoff = -99
            for lyr in range(top_wat-1,number_layers-1):
                if nodes_temp[lyr*number_nodes+number_nodes-1] > min_temp and nodes_temp[(lyr+1)*number_nodes+number_nodes-1] > min_temp:
                    Qavail[t+1,0] += np.pi * r**2 * RLW * disDT * 4.181 * 1000. / 3600. #kWh (fixed delta T)
                elif nodes_temp[lyr*number_nodes+number_nodes-1] > min_temp: # final layer above min_temp
                    fac = (0.5 + ((nodes_temp[lyr*number_nodes+number_nodes-1] - min_temp) / (nodes_temp[lyr*number_nodes+number_nodes-1] - nodes_temp[(lyr+1)*number_nodes+number_nodes-1]))) # proportion of layer and layer below above min_temp assuming linear temperature drop between layer temperatures
                    Qavail[t+1,0] += fac * np.pi * r**2 * RLW * disDT * 4.181 * 1000. / 3600. #kWh (fixed delta T)
                    cutoff = lyr - top_wat + 2 # final layer above min_temp
                    
            if nodes_temp[(number_layers-1)*number_nodes+number_nodes-1] > min_temp:
                Qavail[t+1,0] += np.pi * r**2 * RLW * disDT * 4.181 * 1000. / 3600. #kWh (fixed delta T)
                        
            if np.remainder(t,ts) == 0:
                f_n_t[tr,:] = nodes_temp # Full Temperature Results      
                tr = tr + 1
    
    # %% Generate outputs
    #%%% Create folder in specific location
    def process_folder(location, folder_name):
        folder_path = os.path.join(location, folder_name)
        # Check if the folder exists
        folder_exists = os.path.exists(folder_path)
        if not folder_exists:
            os.makedirs(folder_path)
    
    # Set base location and define folder names
    # base_location = "/Users/cmb22235/OneDrive - University of Strathclyde/Desktop/STEaM WP4 team/MSTES-insul/results/"
    base_location = "/Users/cmb22235/OneDrive - University of Strathclyde/Desktop/"
    folder_name = f'y{Mod_Year},£{surp_tar}&{grid_tar},{int(h)}m,O{Option},M{CHPmode},dT{int(disDT)},hT{int(heat_temp)},sT{int(store_temp)},mT{int(min_temp)}'
    location_1 = os.path.join(base_location, folder_name)
    process_folder(base_location, folder_name)
    
    folder_full = 'full'
    full_location = os.path.join(location_1, folder_full)
    process_folder(location_1, folder_full)
    
    folder_year = 'first_year'
    year_location = os.path.join(location_1, folder_year)
    process_folder(location_1, folder_year)
    
    folder_year = 'final_year'
    year_location = os.path.join(location_1, folder_year)
    process_folder(location_1, folder_year)
    
    folder_summer = 'summer'
    summer_location = os.path.join(location_1, folder_summer)
    process_folder(location_1, folder_summer)
    
    folder_winter = 'winter'
    winter_location = os.path.join(location_1, folder_winter)
    process_folder(location_1, folder_winter)
    
    folder_summer3d = 'summer3d'
    summer3d_location = os.path.join(location_1, folder_summer3d)
    process_folder(location_1, folder_summer3d)
    
    folder_winter3d = 'winter3d'
    winter3d_location = os.path.join(location_1, folder_winter3d)
    process_folder(location_1, folder_winter3d)
    
    results_location = base_location + folder_name
        
    #%% Results
    #%%% Print node temp
    f_s_t = f_n_t[:,number_nodes-1:number_layers*number_nodes+1:number_nodes] # Store Air & Water Temperature Results
    f_s_t = f_s_t[:,top_wat-1:] # Store Water Temperature Results
    f_i_t = f_n_t[:,number_nodes-2:number_layers*number_nodes+1:number_nodes] # Insulation Temperature Results
    f_i_t = f_i_t[:,top_wat-1:] # Insulation Temperature Results
    f_c_t = f_n_t[:,number_nodes-3:number_layers*number_nodes+1:number_nodes] # Concrete Shaft Wall Temperature Results
    f_c_t = f_c_t[:,top_wat-1:] # Concrete Shaft Wall Temperature Results
    f_g_t = f_n_t[:,number_nodes-4:number_layers*number_nodes+1:number_nodes] # Concrete Shaft Wall Temperature Results
    f_g_t = f_g_t[:,top_wat-1:] # Concrete Shaft Wall Temperature Results 

    if RLW > 0:    
        header_node = [f"L{layer}_N{node}" for layer in range(number_layers) for node in range(number_nodes)]
        header = np.array(header_node)
        # filename = f"f_n_t {current_time}.csv"
        filename = "f_n_t.csv"
        np.savetxt(os.path.join(location_1, filename), 
                    f_n_t, delimiter=',', header=', '.join(header), fmt='%f', comments='')

    #%%% RES - Modelled years only (kWh basis)
    # Reshape the variables into hourly data (8760 data per year)
    eWdS = eWdS.reshape(-1, ts).sum(axis=1).reshape(-1, 1)
    DemH = DemH.reshape(-1, ts).sum(axis=1).reshape(-1, 1)
    chp2HD = chp2HD.reshape(-1, ts).sum(axis=1).reshape(-1, 1)
    chp2tes = chp2tes.reshape(-1, ts).sum(axis=1).reshape(-1, 1)
    hS2H = hS2H.reshape(-1, ts).sum(axis=1).reshape(-1, 1) 
    gb2HD = gb2HD.reshape(-1, ts).sum(axis=1).reshape(-1, 1)
    waste_heat = waste_heat.reshape(-1, ts).sum(axis=1).reshape(-1, 1)
    Un_H = Un_H.reshape(-1, ts).sum(axis=1).reshape(-1, 1)    
    CHPrem = CHPrem.reshape(-1, ts).sum(axis=1).reshape(-1, 1) 
    GBrem = GBrem.reshape(-1, ts).sum(axis=1).reshape(-1, 1)
    chp_fuel2H = chp_fuel2H.reshape(-1, ts).sum(axis=1).reshape(-1, 1)
    chp_fuel2S = chp_fuel2S.reshape(-1, ts).sum(axis=1).reshape(-1, 1)
    chp_fuel_used = chp_fuel_used.reshape(-1, ts).sum(axis=1).reshape(-1, 1)
    gb_fuel2H = gb_fuel2H.reshape(-1, ts).sum(axis=1).reshape(-1, 1)
    gb_fuel2S = gb_fuel2S.reshape(-1, ts).sum(axis=1).reshape(-1, 1)
    gb_fuel_used = gb_fuel_used.reshape(-1, ts).sum(axis=1).reshape(-1, 1)
    chp2g = chp2g.reshape(-1, ts).sum(axis=1).reshape(-1, 1)
    chp_export_income = chp_export_income.reshape(-1, ts).sum(axis=1).reshape(-1, 1) 
    HeatOpex = HeatOpex.reshape(-1, ts).sum(axis=1).reshape(-1, 1) 
    THL = THL.reshape(-1, ts).sum(axis=1).reshape(-1, 1)
    
    # List of source arrays # no THL
    RES_arrays = [eWdS, DemH, chp2HD, chp2tes, hS2H, gb2HD, waste_heat, Un_H, CHPrem, GBrem,
                  chp_fuel2H, chp_fuel2S, chp_fuel_used, gb_fuel2H, gb_fuel2S, gb_fuel_used,
                  chp2g, chp_export_income, HeatOpex, THL] 
    
    RES = np.zeros((HOURS_PER_YEAR * Mod_Year,len(RES_arrays)))
    
    # Assign each column in RES
    for i, source in enumerate(RES_arrays):
        RES[:, i] = source[:, 0]
        
    header = ['Available Wind Surplus (kWh)',
              'Heat Demand',
              'Heat CHP to HD',
              'Heat CHP to Store',
              'Heat Store to HD',
              'Heat GB to HD',
              'CHP waste heat',
              'Unmet HD',
              'Residual CHP',
              'Residual GB',
              'CHP fuel to HD',
              'CHP fuel to Store',
              'CHP fuel used',
              'GB fuel to HD',
              'GB fuel to Store',
              'GB fuel used',
              'Elec CHP to grid',
              'CHP export income (£)',              
              'Heat Opex (£)',
              'Total store heat loss'
              ]
    if len(RES[0,:]) != len(header):
        stop(f'Invalid column RES {len(RES[0,:])}. Expected {len(header)}.')
    filename = "RES.csv"
    np.savetxt(os.path.join(location_1, filename), 
                RES, delimiter=',', header=', '.join(header), encoding='utf-8-sig', fmt='%f', comments='')    
    
    #%%% tRES - Modelled years only (MWh basis)
    tRES_arrays = [Mod_Year, h, surp_tar, grid_tar, chp_fuel_cost, gb_fuel_cost,
                   Size_CHP, Size_GB, disDT, Option, CHPmode, 
                   eff_GB_thermal, eff_CHP_electrical, eff_CHP_thermal,
                   heat_temp, min_temp, store_temp,
            (np.sum(DemH) / 1000),  # Total heat demand (MWh)
            (np.sum(chp2g) / 1000),  # Total CHP exports elec2grid
            ((np.sum(chp_fuel_used) + np.sum(gb_fuel_used)) / 1000),  # Total fuel used
            (np.sum(hS2H) / 1000),  # Total store2HD
            ((np.sum(chp_fuel2S) + np.sum(gb_fuel2S)) / 1000),  # Total fuel2store
            (np.sum(chp2HD) / 1000),  # Total chp2HD
            ((np.sum(chp_fuel2H) + np.sum(gb_fuel2H)) / 1000),  # Total fuel2heat
            (np.sum(THL) / 1000),  # Total heat losses from store
            (np.sum(chp2tes) / 1000),  # Total heat charge to store
            (np.sum(Un_H) / 1000)  # Total unmet heat demand         
             ]
    
    tRES = np.zeros((1,len(tRES_arrays)))
    
    # Assign each column in tRES
    for i in range (len(tRES_arrays)):
        tRES[0, i] = tRES_arrays[i]
        
    header = ['Modelling years',
              'Thermal store height (m)',
              f'Surplus Tariff ({GBP}/kWh)',
              f'Non-surplus Tariff ({GBP}/kWh)',
              f'chpfuelcost ({GBP}/kWh)',
              f'gbfuelcost ({GBP}/kWh)',
              'CHP size (kW)',
              'GB size (kW)',
              f'Store Discharge DeltaT {DegC}',
              'System option',
              'CHPmode',
              'eff GB thermal (%)',
              'eff chp elec (%)',
              'eff chp thermal (%)',
              f'Heat Demand supply temp {DegC}',
              f'Minimum store supply temp {DegC}',
              f'Heat pump flow temp to fill store {DegC}',
        
              'Total heat demand (MWh)',
              'Total CHP exports elec2grid',
              'Total fuel used',
              'Total store2HD',
              'Total fuel2store',
              'Total chp2HD',
              'Total fuel2heat',
              'Total heat losses from store',
              'Total heat charge to store',
              'Total unmet heat demand'
              ]
    
    if len(tRES[0,:]) != len(header):
        stop(f'Invalid column tRES {len(tRES[0,:])}. Expected {len(header)}.')
        
    filename = "tRES.csv"
    np.savetxt(os.path.join(location_1, filename), 
                tRES, delimiter=',', header=', '.join(header), encoding='utf-8-sig', fmt='%f', comments='')    
    
    #%%% LCOH
    for yr1 in range(Mod_Year):
        OpH[yr1] = np.sum(HeatOpex[yr1*HOURS_PER_YEAR*ts:(yr1+1)*HOURS_PER_YEAR*ts,0]) # Op cost for modelled years
    for yr2 in range(Mod_Year,lifetime):
        OpH[yr2] = OpH[Mod_Year-1] # Copy last op cost for unmodelled years
    
    CAPEX = Size_CHP * CHP_CAPEX + Size_GB * GB_CAPEX + TS_volume * TS_CAPEX # Total Capital Expenses (£)
    CHPEX = Size_CHP * CHP_CAPEX
    
    Maintenance = CAPEX * 0.02 # O&M costs usually 2-5% of CAPEX (£/yr.)
    
    # Initialize arrays for lifetime
    CP, OP, OM, ET = [np.zeros(lifetime) for _ in range(4)]
    CapI = CAPEX # Capital Costs
    if lifetime > 20:
        CP[20] = CHPEX # replacement heat pumps after year 20
    if lifetime > 40:
        CP[40] = CHPEX # replacement heat pumps after year 40
    
    for yr in range(lifetime):
        CP[yr] = CP[yr] / ((1+DR)**(yr+1)) # Additional capital after operation begins
        operating_cost = OpH[yr]  # Extract the scalar value from the 2D array
        OP[yr] = operating_cost / ((1 + DR) ** (yr + 1))  # Discounted operating cost
        OM[yr] = Maintenance / ((1+DR)**(yr+1)) # Maintenance
        ET[yr] = np.sum(DemH[HOURS_PER_YEAR * yr : HOURS_PER_YEAR * (yr + 1)]) / ((1 + DR) ** (yr + 1)) # Yearly HD
    LCOH = (np.sum(CP) + np.sum(OP) + np.sum(OM)) / np.sum(ET) # £/kWh
    
    #%%% KPIs - Lifetime basis      
    # List of arrays to be extended
    arrays_to_extend = [eWdS, DemH, chp2HD, chp_fuel2H, chp2tes,chp_fuel2S, hS2H, gb2HD, waste_heat, Un_H, CHPrem, GBrem, 
                  chp_fuel_used, chp2tes, chp2g, HeatOpex, chp_export_income, THL] 
    
    # Loop through each year in the extended lifetime
    for exyr in range(lifetime - Mod_Year):
        # Calculate the starting index for slicing
        start_idx = int(nt - (nt / Mod_Year))
        # Extend each array by appending the last year's data
        for i, array in enumerate(arrays_to_extend):
            arrays_to_extend[i] = np.append(array, array[start_idx:nt])
            
    # Unpack the extended arrays
    eWdS, DemH, chp2HD, chp_fuel2H, chp2tes, chp_fuel2S, hS2H, gb2HD, waste_heat, Un_H, CHPrem, GBrem, \
    chp_fuel_used, chp2tes, chp2g, HeatOpex, chp_export_income, THL = arrays_to_extend        

    KPI_arrays = [lifetime, h, surp_tar, grid_tar, chp_fuel_cost, gb_fuel_cost,
                   Size_CHP, Size_GB, disDT, Option, CHPmode, 
                   eff_GB_thermal, eff_CHP_electrical, eff_CHP_thermal,
                   heat_temp, min_temp, store_temp,
                   
    (np.sum(DemH)+np.sum(chp2g)) / (np.sum(chp_fuel_used) + np.sum(gb_fuel_used)), # Overall COP (kWh th/kWh)
    np.sum(hS2H) / (np.sum(chp_fuel2S) + np.sum(gb_fuel2S)), # Store Heat COP (kWh th/kWh)
    (np.sum(chp2HD) + np.sum(gb2HD)) / (np.sum(chp_fuel2H) + np.sum(gb_fuel2H)), # HD Direct Supply COP (kWh th/kWh)
    ((np.sum(hS2H) / np.sum(DemH)) * 100), # Proportion of heat via store (%)
    ((np.sum(chp2HD) / np.sum(DemH)) * 100), # Proportion of heat via chp (%)
    ((np.sum(chp2g) / np.sum(chp_fuel_used)) * 100), # Proportion of elec via chp (%)
    ((np.sum(THL) / np.sum(chp2tes)) * 100), # Store heat loss proportion (%)
    LCOH # LCOH (£/kWh)
   ]
    
    KPI = np.zeros((1,len(KPI_arrays)))
    
    # Assign each column in KPI
    for i in range (len(KPI_arrays)):
        KPI[0, i] = KPI_arrays[i]
        
    header = ['Lifetime years',
              'Thermal store height (m)',
              f'Surplus Tariff ({GBP}/kWh)',
              f'Non-surplus Tariff ({GBP}/kWh)',
              f'chpfuelcost ({GBP}/kWh)',
              f'gbfuelcost ({GBP}/kWh)',
              'CHP size (kW)',
              'GB size (kW)',
              f'Store Discharge DeltaT {DegC}',
              'System option',
              'CHPmode',
              'eff GB thermal (%)',
              'eff chp elec (%)',
              'eff chp thermal (%)',
              f'Heat Demand supply temp {DegC}',
              f'Minimum store supply temp {DegC}',
              f'Heat pump flow temp to fill store {DegC}',

              'Overall COP (kWh/kWh)',
              'Store Heat COP (kWh th/kWh)',
              'HD Direct Supply COP (kWh th/kWh)',
              'Proportion of heat via store (%)',
              'Proportion of heat via chp (%)',
              'Proportion of elec via chp fuel (%)',
              'Store heat loss proportion (%)',
              f'LCOH ({GBP}/kWh)'
              ]
    
    if len(KPI[0,:]) != len(header):
        stop(f'Invalid column KPI {len(KPI[0,:])}. Expected {len(header)}.')    
        
    filename = "kpi.csv"
    np.savetxt(os.path.join(location_1, filename), 
                KPI, delimiter=',', header=', '.join(header), encoding='utf-8-sig', fmt='%f', comments='')
        
    #%%% Plot graphs
    plt.style.use('ggplot')
    plt.rcParams.update({'font.size': 6})   
        
    # plot colour   
    c = np.array([
    "#c1272d", # red
    "#0000a7", # blue 
    "#008176", # green
    "#ba03af", # purple
    "#eecc16", # yellow
    "#000000", # black                              
    "#854802", # brown
    "#b3b3b3"  # grey
    ])
           
    def pick_time(start_hour, end_hour):
        """
        Returns a list of time indices for the specified range.
        """
        return list(range(start_hour, end_hour))
    
    # Main plotting function
    def plot_all_graph(start_hour, end_hour, x_ticks, x_labels, duration, xaxis):
        """
        Plots the data for the specified time range.
        """
        # Extract data for the specified time range
        time_indices = pick_time(start_hour, end_hour)
        plot_f_at_t(time_indices, x_ticks, x_labels, duration, xaxis) # all top node layer
        plot_f_ab_t(time_indices, x_ticks, x_labels, duration, xaxis) # all top node layer
        plot_f_s_t(time_indices, x_ticks, x_labels, duration, xaxis) # shaftstore water layer
        # plot_f_i_t(time_indices, x_ticks, x_labels, duration, xaxis) # ins layer
        # plot_f_w_t(time_indices, x_ticks, x_labels, duration, xaxis) # wall layer
        # plot_f_g_t(time_indices, x_ticks, x_labels, duration, xaxis) # ground layer
        plot_heat(time_indices, x_ticks, x_labels, duration, xaxis) # Ways to supply heat demand
        plot_store(time_indices, x_ticks, x_labels, duration, xaxis) # Store heat
        plot_elec(time_indices, x_ticks, x_labels, duration, xaxis) # Available Wind Surplus
    
    # Function to plot all top layer temperature
    def plot_f_at_t(time_indices, x_ticks, x_labels, duration, xaxis):
        """
        Plots the all top layer temperature for the specified time range.
        """
        if h > 0:  # If store exists
            fig = plt.figure(figsize=(10,8))      
            legend_labels = []
            for i in range(number_nodes):
                f_at_t = f_n_t[:,i:number_layers*number_nodes+1:number_nodes]
                f_at_t = f_at_t[:,top_wat-1:]
                plt.plot(f_at_t[time_indices,0], ls = "-", lw = "0.5")     
                legend_labels.append(f"N{i} {Rx[i+1]}-{Rx[i]}m")
                plt.legend(legend_labels,title='Node',bbox_to_anchor=(1, 1),loc='upper left', prop={"size": 6}) # Place legend outside
            plt.xlabel(f"Time {xaxis}")
            plt.xticks(x_ticks, x_labels)
            plt.ylabel("Node temperature ({DegC})") 
            plt.title(f"All top layer temp for {int(store_temp)}{DegC} Store({int(h)}m{int(TS_volume)}m3)")

            # Save the plot
            filename = f"{results_location}/{duration}/f_at_t_{duration}.png"
            plt.savefig(filename, format="png", dpi=300, bbox_inches="tight")
            fig.clear()
            plt.close(fig)     
            
    # Function to plot all btm layer temperature
    def plot_f_ab_t(time_indices, x_ticks, x_labels, duration, xaxis):
        """
        Plots the all btm layer temperature for the specified time range.
        """
        if h > 0:  # If store exists
            fig = plt.figure(figsize=(10,8))  
            legend_labels = []
            for i in range(number_nodes):
                f_ab_t = f_n_t[:,i:number_layers*number_nodes+1:number_nodes]
                f_ab_t = f_ab_t[:,top_wat-1:]
                plt.plot(f_ab_t[time_indices,-1], ls = "-", lw = "0.5")     
                legend_labels.append(f"N{i} {Rx[i+1]}-{Rx[i]}m")
                plt.legend(legend_labels,title='Node',bbox_to_anchor=(1, 1),loc='upper left', prop={"size": 6}) # Place legend outside
            plt.xlabel(f"Time {xaxis}")
            plt.xticks(x_ticks, x_labels)
            plt.ylabel("Node temperature ({DegC})") 
            plt.title(f"All btm layer temp for {int(store_temp)}{DegC} Store({int(h)}m{int(TS_volume)}m3)")

            # Save the plot
            filename = f"{results_location}/{duration}/f_ab_t_{duration}.png"
            plt.savefig(filename, format="png", dpi=300, bbox_inches="tight")
            fig.clear()
            plt.close(fig)     

    # Function to plot water temperature
    def plot_f_s_t(time_indices, x_ticks, x_labels, duration, xaxis):
        """
        Plots the water temperature for the specified time range.
        """
        if h > 0:  # If store exists
            fig = plt.figure(figsize=(10,8))
            plt.plot(f_s_t[time_indices], ls="-", lw="0.5")
            plt.title(f"Water({Rx[-1]}-{Rx[-2]}m) Temp for {int(store_temp)}{DegC} Store({int(h)}m{int(TS_volume)}m3)")
            plt.xticks(x_ticks, x_labels)
            plt.xlabel(f"Time {xaxis}")
            plt.ylabel(f"Node temperature for each layer ({DegC})")
    
            # Add legend
            legend_labels = [f"L{i+1}" for i in range(number_layers-(top_wat-1))]  # L1 to L10
            plt.legend(legend_labels, bbox_to_anchor=(1, 1), loc='upper left', prop={"size": 6})  # Place legend outside
    
            # Save the plot
            filename = f"{results_location}/{duration}/f_s_t_{duration}.png"
            plt.savefig(filename, format="png", dpi=300, bbox_inches="tight")
            fig.clear()
            plt.close(fig) 
            
    # Function to plot insulation temperature
    def plot_f_i_t(time_indices, x_ticks, x_labels, duration, xaxis):
        """
        Plots the insulation temperature for the specified time range.
        """
        if h > 0:  # If store exists
            fig = plt.figure(figsize=(10,8))
            plt.plot(f_i_t[time_indices], ls="-", lw="0.5")
            plt.title(f"Insulation({Rx[-2]}-{Rx[-3]}m) Temp for {int(store_temp)}{DegC} Store({int(h)}m{int(TS_volume)}m3)")
            plt.xticks(x_ticks, x_labels)
            plt.xlabel(f"Time {xaxis}")
            plt.ylabel(f"Node temperature for each layer ({DegC})")
    
            # Add legend
            legend_labels = [f"L{i+1}" for i in range(number_layers-4)]  # L1 to L10
            plt.legend(legend_labels, bbox_to_anchor=(1, 1), loc='upper left', prop={"size": 6})  # Place legend outside
    
            # Save the plot
            filename = f"{results_location}/{duration}/f_i_t_{duration}.png"
            plt.savefig(filename, format="png", dpi=300, bbox_inches="tight")
            fig.clear()
            plt.close(fig) 
                
    # Function to plot wall temperature
    def plot_f_w_t(time_indices, x_ticks, x_labels, duration, xaxis):
        """
        Plots the wall temperature for the specified time range.
        """
        if h > 0:  # If store exists
            fig = plt.figure(figsize=(10,8))
            plt.plot(f_c_t[time_indices], ls="-", lw="0.5")
            plt.title(f"Wall({Rx[-3]}-{Rx[-4]}m) Temp for {int(store_temp)}{DegC} Store({int(h)}m{int(TS_volume)}m3)")
            plt.xticks(x_ticks, x_labels)
            plt.xlabel(f"Time {xaxis}")
            plt.ylabel(f"Node temperature for each layer ({DegC})")
    
            # Add legend
            legend_labels = [f"L{i+1}" for i in range(number_layers-4)]  # L1 to L10
            plt.legend(legend_labels, bbox_to_anchor=(1, 1), loc='upper left', prop={"size": 6})  # Place legend outside
    
            # Save the plot
            filename = f"{results_location}/{duration}/f_w_t_{duration}.png"
            plt.savefig(filename, format="png", dpi=300, bbox_inches="tight")
            fig.clear()
            plt.close(fig)     
            
    # Function to plot ground temperature
    def plot_f_g_t(time_indices, x_ticks, x_labels, duration, xaxis):
        """
        Plots the ground temperature for the specified time range.
        """
        if h > 0:  # If store exists
            fig = plt.figure(figsize=(10,8))
            plt.plot(f_g_t[time_indices], ls="-", lw="0.5")
            plt.title(f"Ground({Rx[-4]}-{Rx[-5]}m) Temp for {int(store_temp)}{DegC} Store({int(h)}m{int(TS_volume)}m3)")
            plt.xticks(x_ticks, x_labels)
            plt.xlabel(f"Time {xaxis}")
            plt.ylabel(f"Node temperature for each layer ({DegC})")
    
            # Add legend
            legend_labels = [f"L{i+1}" for i in range(number_layers-4)]  # L1 to L10
            plt.legend(legend_labels, bbox_to_anchor=(1, 1), loc='upper left', prop={"size": 6})  # Place legend outside
    
            # Save the plot
            filename = f"{results_location}/{duration}/f_g_t_{duration}.png"
            plt.savefig(filename, format="png", dpi=300, bbox_inches="tight")
            fig.clear()
            plt.close(fig) 

    def plot_heat(time_indices, x_ticks, x_labels, duration, xaxis): 
        fig = plt.figure(figsize=(10,8))
        
        plt.subplot(4,1,1)
        plt.plot(DemH[time_indices], label='HD', ls = '-', lw = '0.5', c=c[0])
        plt.title(f"Heat demand {duration}")
        plt.ylabel("Power (kW th)")
        plt.xticks(ticks=x_ticks, labels="")
        
        plt.subplot(4,1,2)
        plt.plot(chp2HD[time_indices], label='HD direct from chp', ls = '-', lw = '0.5', c=c[1])
        plt.title("HD direct from chp")
        plt.ylabel("Power (kW th)")
        plt.xticks(ticks=x_ticks, labels="")
        
        plt.subplot(4,1,3)
        plt.plot(hS2H[time_indices], label='Store heat to HD', ls = '-', lw = '0.5', c=c[2])
        plt.title("Store heat to HD")
        plt.ylabel("Power (kW th)")
        plt.xticks(ticks=x_ticks, labels="")
        
        plt.subplot(4,1,4)
        plt.plot(gb2HD[time_indices], label='HD direct from gb', ls = '-', lw = '0.5', c=c[3])
        plt.title("HD direct from gb")
        plt.ylabel("Power (kW th)")
        plt.xticks(ticks=x_ticks, labels=x_labels)
        plt.xlabel(f"Time {xaxis}")

        filename = f'{results_location}/{duration}/heat {duration}.png'
        plt.savefig(filename, format = 'png',dpi=300, bbox_inches='tight')
        fig.clear()
        plt.close(fig)  
        
    def plot_store(time_indices, x_ticks, x_labels, duration, xaxis): 
        fig = plt.figure(figsize=(10,8))

        plt.subplot(3,1,1)
        plt.plot(chp2tes[time_indices], label='Heat charge to store', ls = '-', lw = '0.5', c=c[0])
        plt.title("Heat charge to store")
        plt.ylabel("Power (kW th)")
        plt.xticks(ticks=x_ticks, labels="")

        plt.subplot(3,1,2)
        plt.plot(hS2H[time_indices], label='Store heat to HD', ls = '-', lw = '0.5', c=c[2])
        plt.title("Store heat to HD")
        plt.ylabel("Power (kW th)")
        plt.xticks(ticks=x_ticks, labels="")        
        
        plt.subplot(3,1,3)
        plt.plot(THL[time_indices], label='Heat losses from store', ls = '-', lw = '0.5', c=c[1])
        plt.title("Heat losses from store")
        plt.ylabel("Power (kW th)")
        plt.xticks(ticks=x_ticks, labels=x_labels)
        plt.xlabel(f"Time {xaxis}")
                
        filename = f'{results_location}/{duration}/store {duration}.png'
        plt.savefig(filename, format = 'png',dpi=300, bbox_inches='tight')
        fig.clear()
        plt.close(fig) 
        
    def plot_elec(time_indices, x_ticks, x_labels, duration, xaxis): 
        fig = plt.figure(figsize=(10,8))
        
        plt.subplot(4,1,1)
        plt.plot(eWdS[time_indices], label='Wind surplus', ls = '-', lw = '0.5', c = c[0])
        plt.title(f"Wind surplus {duration}")
        plt.ylabel("Power (kW e)")
        plt.xticks(ticks=x_ticks, labels="")
        
        plt.subplot(4,1,2)
        plt.plot(chp2g[time_indices], label='CHP exports elec2grid', ls = '-', lw = '0.5', c = c[1])
        plt.title("CHP exports elec2grid")
        plt.ylabel("Power (kW e)")
        plt.xticks(ticks=x_ticks, labels="")        
        
        plt.subplot(4,1,3)
        plt.plot(chp_fuel_used[time_indices], label='CHP fuel used', ls = '-', lw = '0.5', c = c[2])
        plt.title("CHP fuel used")        
        plt.ylabel("Power (kW e)")
        plt.xticks(ticks=x_ticks, labels="")
        
        plt.subplot(4,1,4)
        plt.plot(gb_fuel_used[time_indices], label='GB fuel used', ls = '-', lw = '0.5', c = c[3])
        plt.title("GB fuel used")        
        plt.ylabel("Power (kW e)")
        plt.xticks(ticks=x_ticks, labels=x_labels)
        plt.xlabel(f"Time {xaxis}")
        
        filename = f'{results_location}/{duration}/elec {duration}.png'
        plt.savefig(filename, format = 'png',dpi=300, bbox_inches='tight')
        fig.clear()
        plt.close(fig)    
 
    #%%% Function to pick time indices 
    # Time For Each Season
    # 1 day  = 24h * data/hr = 24
    # 3 day  = 24*3 = 72  
    # 1 week = 24*7 = 168
    # 1 year = 365 days * 24hr = 8760 hr * 1 data per hour -> 8760
    # data start from [0]
    if HOURS_PER_YEAR == 8760:
        SEASON_RANGES = {
        "spring": (1416, 3624),  # March 1 to May 31
        "summer": (3624, 6552),  # June 1 to August 31
        "fall": (6552, 8016),    # September 1 to November 30
        "winter": [(8016, 8760),(0,1416)]  # December 1 to February 28/29
        }
        # Constants
        MONTH_TICKS = [0, 744, 1416, 2160, 2880, 3624, 4344, 5088, 5832, 6552, 7296, 8016]  # Approximate hours per month
        # MONTH_TICKS = [tick * ts for tick in MONTH_TICKS]
        MONTH_LABELS = ['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec']    
    elif HOURS_PER_YEAR == 17520:
        SEASON_RANGES = {
        "spring": (2832, 7248),    # March 1 to May 31 (1416*2 to 3624*2)
        "summer": (7248, 13104),   # June 1 to August 31 (3624*2 to 6552*2)
        "fall": (13104, 16032),    # September 1 to November 30 (6552*2 to 8016*2)
        "winter": [(16032, 17520), (0, 2832)]  # December 1 to February 28/29 (8016*2 to 8760*2 and 0 to 1416*2)
        }
        # Constants
        MONTH_TICKS = [0, 1488, 2832, 4320, 5760, 7248, 8688, 10176, 11664, 13104, 14592, 16032] # Half-hourly intervals
        # MONTH_TICKS = [tick * ts for tick in MONTH_TICKS]
        MONTH_LABELS = ['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec']
    
    # Extend MONTH_TICKS and MONTH_LABELS for multiple years
    extended_month_ticks = []; extended_month_labels = []
    for year in range(Mod_Year+1):
        # Add the offset for each year to MONTH_TICKS
        extended_month_ticks.extend([tick + year * HOURS_PER_YEAR for tick in MONTH_TICKS])
        # Repeat the MONTH_LABELS for each year
        extended_month_labels.extend(MONTH_LABELS)
    full_month_ticks = []; full_month_labels = []
    for year in range(Mod_Year+1):     
        full_month_ticks.append(HOURS_PER_YEAR * year)
        full_month_labels.append(str(year))
        
    # Helper function to generate ticks and labels
    def generate_ticks_labels(total_steps, step_size, label_prefix=""):
        """
        Generates ticks and labels for the x-axis.
        """
        ticks = [step_size * i for i in range(total_steps + 1)]
        labels = [f"{label_prefix}{i}" for i in range(total_steps + 1)]
        return ticks, labels
    
    def plot_full():
        """
        Plots the entire dataset over all years.
        """
        # Calculate start and end hours for all years
        start_hour = 0
        end_hour = HOURS_PER_YEAR * Mod_Year
            
        # Use predefined monthly ticks and labels
        if Mod_Year <= 2:
            x_ticks = extended_month_ticks
            x_labels = extended_month_labels
        else:
            x_ticks = full_month_ticks
            x_labels = full_month_labels
            
        # Plot the data
        plot_all_graph(start_hour, end_hour, x_ticks, x_labels, 'full', '(years)')     
            
    # Function to plot one year of data
    def plot_first_year():
        """
        Plots the data for first full year.
        """
        # Calculate start and end hours for the specified year
        start_hour = 0
        end_hour = HOURS_PER_YEAR

        # Use predefined monthly ticks and labels
        x_ticks = MONTH_TICKS
        x_labels = MONTH_LABELS

        # Plot the data
        plot_all_graph(start_hour, end_hour, x_ticks, x_labels, "first_year", "(months)")
        
    # Function to plot one year of data
    def plot_final_year():
        """
        Plots the data for a full year.
        """
        # Calculate start and end hours for the specified year
        start_hour = HOURS_PER_YEAR * (Mod_Year - 1)
        end_hour = HOURS_PER_YEAR * Mod_Year

        # Use predefined monthly ticks and labels
        x_ticks = MONTH_TICKS
        x_labels = MONTH_LABELS

        # Plot the data
        plot_all_graph(start_hour, end_hour, x_ticks, x_labels, "final_year", "(months)")

    def plot_summer():
        """
        Plots the data for the summer period of the final year.
        """
        start_hour = HOURS_PER_YEAR * (Mod_Year - 1) + 3624  # Start of summer (June)
        end_hour = HOURS_PER_YEAR * (Mod_Year - 1) + 5832   # End of summer (August)
        
        # Generate ticks and labels    
        total_weeks = int((end_hour - start_hour) / (24 * 7))  # 1 data point/hour * 24 hours/day * 7 days/week
        summer_ticks, summer_labels = generate_ticks_labels(total_weeks, 24 * 7, "W")
        
        # Plot the data
        plot_all_graph(start_hour, end_hour, summer_ticks, summer_labels, 'summer', '(weeks)')
    
    def plot_winter():
        """
        Plots the data for the winter period of the second-to-last year.
        """
        if Mod_Year > 1:
            start_hour = HOURS_PER_YEAR * (Mod_Year - 2) + 8016  # Start of winter (December)
            end_hour = HOURS_PER_YEAR * (Mod_Year - 2) + 8016 + 1416  # End of winter (February)
        elif Mod_Year == 1:
            start_hour = 8016  # Start of winter (December)
            end_hour = HOURS_PER_YEAR  # End of winter (February)
        
        # Generate ticks and labels
        total_weeks = int((end_hour - start_hour) / (24 * 7))
        winter_ticks, winter_labels = generate_ticks_labels(total_weeks, 24 * 7, "W")
    
        # Plot the data
        plot_all_graph(start_hour, end_hour, winter_ticks, winter_labels, 'winter', '(weeks)')
        
    # Function to plot summer 3D data
    def plot_summer_3d():
        """
        Plots the data for a 3-day period in last summer.
        """
        start_hour = HOURS_PER_YEAR * (Mod_Year - 1) + 4344 + 24 * 28  # Start of 3-day period
        end_hour = HOURS_PER_YEAR * (Mod_Year - 1) + 4344 + 24 * (28 + 3)  # End of 3-day period

        # Generate ticks and labels
        total_days = int((end_hour - start_hour) / 24)
        summer_3d_ticks, summer_3d_labels = generate_ticks_labels(total_days, 24, "D")

        # Plot the data
        plot_all_graph(start_hour, end_hour, summer_3d_ticks, summer_3d_labels, "summer3d", "(days)")
    
    # Function to plot winter 3D data
    def plot_winter_3d():
        """
        Plots the data for a 3-day period in last winter.
        """
        start_hour = HOURS_PER_YEAR * (Mod_Year - 1) + 8016 + 24 * 27  # Start of 3-day period
        end_hour = HOURS_PER_YEAR * (Mod_Year - 1) + 8016 + 24 * (27 + 3)  # End of 3-day period

        # Generate ticks and labels
        total_days = int((end_hour - start_hour) / 24)
        winter_3d_ticks, winter_3d_labels = generate_ticks_labels(total_days, 24, "D")

        # Plot the data
        plot_all_graph(start_hour, end_hour, winter_3d_ticks, winter_3d_labels, "winter3d", "(days)")
            
    #%%% Generate Plots for speific duration
    plot_full()
    plot_first_year()
    plot_final_year()
    plot_summer()
    plot_winter()    
    plot_summer_3d()
    plot_winter_3d()
    
    #%%% Plot ground heating envelope
    if h > 0:  # If store exists
        your_array = f_n_t  # Ensure this is a 2D array (rows x columns)
    
        def calculate_column_averages_with_similar_names_and_rows(array, start_row, end_row):
            column_sums = {}
            column_counts = {}
            # Iterate through headers to find columns with similar names after '_'
            headers = header_node  # Ensure headers are defined
            for header in headers:
                if '_' in header:
                    _, suffix = header.split('_')
                    if suffix not in column_sums:
                        column_sums[suffix] = 0
                        column_counts[suffix] = 0
    
            # Iterate through rows within the specified range
            for current_row in range(start_row, end_row):
                # Iterate through columns with similar names after '_'
                for col, header in enumerate(headers):
                    if '_' in header:
                        prefix, suffix = header.split('_')
                        if suffix in column_sums:
                            try:
                                column_sums[suffix] += float(array[current_row][col])
                                column_counts[suffix] += 1
                            except ValueError:
                                pass  # Skip non-numeric values
    
            # Calculate averages for each column with similar names
            column_averages = {}
            for suffix, sum in column_sums.items():
                count = column_counts[suffix]
                column_averages[suffix] = sum / count if count != 0 else 0
            return column_averages
    
        # Calculate averages for each year
        averages_all = []
        for year in range(Mod_Year):
            start_row = year * 8759  # Starting row (inclusive)
            end_row = (year + 1) * 8759  # Ending row (exclusive)
            averages = calculate_column_averages_with_similar_names_and_rows(your_array, start_row, end_row)
            averages_all.append(averages)
    
        # Custom x-axis array (reversed)
        x_axis_values = Rx[::-1]  # Ensure Rx is defined
    
        # Plotting
        fig = plt.figure(figsize=(10,8))
        for i, average in enumerate(averages_all):
            suffixes = list(average.keys())[::-1]  # Reverse the order of the keys
            values = [average[suffix] for suffix in suffixes]
            
            # Ensure x_axis_values and values have the same length
            if len(x_axis_values) != len(values):
                x_axis_values = x_axis_values[:len(values)]  # Trim x_axis_values to match
            
            plt.plot(x_axis_values, values, label=f'Year {i + 1}')
    
        plt.xlabel('Distance from center of MSTES (m)')
        plt.ylabel(f'Mean temperature ({DegC})')
        plt.title(f'Ground heating envelope over {Mod_Year}yrs {int(store_temp)}{DegC} Store({int(h)}m{int(TS_volume)}m3)')
        # plt.legend(loc='best', ncol=4, fancybox=True, prop={'size': 8})
        plt.legend(loc='best', fancybox=True, prop={'size': 8})
        plt.xlim(right=100)  # Adjust x-axis limit if necessary
    
        # Save the plot
        filename = f'{results_location}/Ghe{Mod_Year}yrs{int(store_temp)}{DegC}.png'
        plt.savefig(filename, format='png', dpi=300, bbox_inches='tight')
        plt.close(fig)
        
#%% Combine kpi_last      
if combine_tRESnkpi == 1:                                
  def combine_csv_files(root_dir, file_name):
    combined_data = pd.DataFrame()
    for subdir, dirs, files in os.walk(root_dir):
      for file in files:
        if file_name in file:
          temp_data = pd.read_csv(os.path.join(subdir, file), encoding='ISO-8859-1')
          combined_data = pd.concat([combined_data, temp_data])
    return combined_data
  root_dir = '/Users/cmb22235/OneDrive - University of Strathclyde/Desktop/STEaM WP4 team/MSTES-insul/results'
  file_name = 'kpi'  # replace with your file name
  combined_data = combine_csv_files(root_dir, file_name)
  combined_data.to_csv('combined_kpi.csv', index=False, encoding='ISO-8859-1')
  file_name = 'tRES'  # replace with your file name
  combined_data = combine_csv_files(root_dir, file_name)
  combined_data.to_csv('combined_tRES.csv', index=False, encoding='ISO-8859-1')
  
