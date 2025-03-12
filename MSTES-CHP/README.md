# How to use
1) Download the all the files.
2) Open all '.csv' files via Microsoft Excel:

   i) 'MSTES-CHP/MSTES-CHP.py' to modify the input data: Current Coylton GSP Electricity Demand (MWh), Coylton 33kV Wind (MWh), West Whitlawburn District Heating (kWh), Air Temp (DegC)

4) Open all '.py' files in an Integrated Development Environment (IDE) for Python (etc. Spyder, Pycharm, Vscode).
   
    i) 'shaftstore_5a.py' to modify the code for the mineshaft thermal store.
   
    ii) 'MSTES-CHP.py' for CHP integration energy system.
   
5) Install all python libraries/packages (numpy==1.26.4, pandas==2.2.2, matplotlib==3.9.2).
6) Change the directory in the 'MSTES-CHP.py' for the inputs (line 13) and outputs (lines 484,1175,178).
7) Change the variables and parameters according to your system design. Main inputs are from line 22 to 79. 
8) Run the code.
9) Check the results (graphs in png. format and data in csv. format) in the 'results' folder.

# Cite
Ewe WE, Tuohy P, Flett G. STEaM source code [python], Github. August 20 2024. https://github.com/winengewe/STEaM
