# How to use
1) Download the all the files.
2) Open all '.csv' files via Microsoft Excel:
   
    i) 'gsp_data_1.csv' to modify the input data: Current Coylton GSP Electricity Demand (MWh), Future Coylton Electricity Demand (MWh), Coylton 33kV Wind (MWh), Coylton GSP Wind (MWh), Domestic Heat per           kWh, Air Temp (DegC), West Whitlawburn District Heating (kWh).
   
    ii) 'geo_win.csv' to modify the ground, concrete, insulation, air water thermal conductivity (W/mK) for each node and layer.
   
    iii) 'geo_win2.csv' to modify the ground heat capacity (J/kgK), ground density (kg/m3), ground porosity (fraction) for each node and layer.
   
4) Open all '.py' files in an Integrated Development Environment (IDE) for Python (etc. Spyder, Pycharm, Vscode).
   
    i) 'shaftstore_1d_0i.py' to modify the code for the mineshaft thermal store.
   
    ii) 'MSTES-HP.py' for HP integration energy system.
   
6) Install all python libraries/packages (numpy==1.26.4, pandas==2.2.2, matplotlib==3.9.2).
7) Change the directory in the 'MSTES-HP.py' for the inputs (lines 46,47,101,103) and outputs (lines 35,1146).
8) Change the variables and parameters according to your system design. Main inputs are from line 39 to 136. 
9) Run the code.
10) Check the results (graphs in png. format and data in csv. format) in the 'results' folder.

# Cite
Ewe WE, Tuohy P, Flett G. STEaM source code [python], Github. August 20 2024. https://github.com/winengewe/STEaM
