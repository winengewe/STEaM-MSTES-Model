import numpy as np

from scipy.integrate import odeint

class ShaftStore(object):

    def __init__(self, Rx, number_layers, number_nodes, top_wat, RLA, RLW, XTC, Geo, mu, ithk, icp, iden, cden):
        
        self.Rx = Rx
    
        self.XTC = XTC
        self.Geo = Geo
    
        self.number_nodes = number_nodes
        self.number_layers = number_layers
        self.top_wat = top_wat - 1 # python -1 conversion
        self.RLA = RLA # Consolidated Air Layer (m)
        self.RLW = RLW # Ground and Water layers height (m)
        self.mu = mu
        self.alp = 0.143e-6
        
        self.icp = icp
        self.iden = iden
        self.cden = cden
        
    def HTC_S(self,layer,Tf,QBest,ds): # Store wall heat transfer coefficients

        DH = ds
        
        if layer < self.top_wat:
            rho = (1.1881 * 293.15) / (273.15 + Tf) # Density
            Rmu = (1.8462 + (4.576e-3 * (Tf-27.))) * 1e-5 # Dynamic viscosity
            k = 0.024 # Air thermal conductivity
            tdf = 18.24e-6 # Air thermal diffusivity
            tec = 0.0034 # volumetric thermal expansion coefficient
            
        else:
            rho = 1001.1 - (0.082 * Tf) - (0.0036 * Tf**2) # Density
            Rmu = (1.825372e-7 * Tf**2) - (3.1195e-5 * Tf) + 0.001639325 # Dynamic viscosity
            k = 0.579 # Water thermal conductivity
            tdf = 1.41e-7 # Water thermal diffusivity
            tec = 0.000214 # volumetric thermal expansion coefficient

        RaMS = (9.81 * tec * abs(QBest) * DH**4) / (k * (Rmu/rho) * tdf) # Rayleigh number - mineshaft
        NuMS = 0.2 * RaMS**0.25 # Nusselt number - mineshaft
        hms = (NuMS * k) / DH # Mineshaft heat transfer coefficient
            
        return hms
    
    def coefficient_base(self, layer, node, nodes_temp, Qp):
        
        Rx = self.Rx
        XTC = self.XTC
        if layer == -1 and node == self.number_nodes - 1:
            RL = self.RLA
        else:
            RL = self.RLW
        
        layer1 = layer
            
        if node == self.number_nodes - 1:
            nnode = (layer * self.number_nodes) + node
            hms = max(1e-10, self.HTC_S(layer1, nodes_temp[nnode], Qp[layer], self.Rx[self.number_nodes - 1]))

        elif node == self.number_nodes - 2:
            nnode = (layer * self.number_nodes) + node + 1
            hms = max(1e-10, self.HTC_S(layer1, nodes_temp[nnode], Qp[layer], self.Rx[self.number_nodes - 1]))
                      
        else:
            hms = 0.
            
        hms = np.copy(hms)
        
        if node == 0: # Outer Ring      
            rf2 = 0.5 * (Rx[node] + Rx[node+1])
            rf3 = 0.5 * (Rx[node+1] + Rx[node+2])
            RTsx = np.log(Rx[node] / rf2)
            RTsy = np.log(rf2 / Rx[node+1])
            RTsz = np.log(Rx[node+1] / rf3)
            Cxy = ((2. * np.pi * XTC[layer1,node] * RL) / RTsx)
            Cyz = ((2. * np.pi * RL) / ((RTsy / XTC[layer1,node]) + (RTsz / XTC[layer1,node+1])))
            
        elif node == self.number_nodes - 2: # Insulation
            rf1 = 0.5 * (Rx[node-1] + Rx[node])
            rf2 = 0.5 * (Rx[node] + Rx[node+1])
            RTsx = np.log(rf1 / Rx[node])
            RTsy = np.log(Rx[node] / rf2)
            RTsz = np.log(rf2 / Rx[node+1])
            Cxy = ((2. * np.pi * RL) / ((RTsx / XTC[layer1,node-1]) + (RTsy / XTC[layer1,node])))           
            Cyz = 1. / ((1./ (hms * 2. * np.pi * Rx[node+1] * RL)) + (RTsz / (2. * np.pi * XTC[layer1,node] * RL)))
 
            # rf1 = 0.5 * (Rx[node-1] + Rx[node])
            # rf2 = 0.5 * (Rx[node] + Rx[node+1])
            # RTsx = np.log(rf1 / Rx[node])
            # RTsy = np.log(Rx[node] / rf2)
            # RTsz = np.log(rf2 / Rx[node+1])
            # Cxy = ((2. * np.pi * RL) / ((RTsx / XTC[layer1,node-1]) + (RTsy / XTC[layer1,node])))
            # Cyz = ((2. * np.pi * RL) / ((RTsy / XTC[layer1,node]) + (RTsz / XTC[layer1,node+1])))
            
        elif node == self.number_nodes - 1: # Fluid                   
            rf1 = 0.5 * (Rx[node-1] + Rx[node])
            RTs = np.log(rf1 / Rx[node])
            Cxy = 1. / ((1./ (hms * 2. * np.pi * Rx[node] * RL)) + (RTs / (2. * np.pi * XTC[layer1,node-1] * RL)))
            Cyz = 0.
            
            # rf1 = 0.5 * (Rx[node-1] + Rx[node])
            # RTsx = np.log(Rx[node-1] / rf1)
            # RTsy = np.log(rf1 / Rx[node])
            # Cxy = ((2. * np.pi * RL) / ((RTsx / XTC[layer1,node-1]) + (RTsy / XTC[layer1,node])))
            # Cyz = 0.
            
        elif node > 0 and node < self.number_nodes - 2: # Intermediate Rings       
            rf1 = 0.5 * (Rx[node-1] + Rx[node])
            rf2 = 0.5 * (Rx[node] + Rx[node+1])
            rf3 = 0.5 * (Rx[node+1] + Rx[node+2])
            RTsw = np.log(rf1 / Rx[node])
            RTsx = np.log(Rx[node] / rf2)
            RTsy = np.log(rf2 / Rx[node+1])
            RTsz = np.log(Rx[node+1] / rf3)
            Cxy = ((2. * np.pi * RL) / ((RTsw / XTC[layer1,node-1]) + (RTsx / XTC[layer1,node])))
            Cyz = ((2. * np.pi * RL) / ((RTsy / XTC[layer1,node]) + (RTsz / XTC[layer1,node+1])))
    
        else:
            Cxy = 0.
            Cyz = 0.
            
        return Cxy, Cyz, hms
    
    def coefficient_A(self, layer, node, Cxy, Cyz, mdotu, mdotd):
    
        if node < self.number_nodes - 3: # Ground
            cp = self.Geo[layer,0] # J/kg Soil Heat Capacity
            Por = self.Geo[layer,2] # Soil porosity
            Den = (self.Geo[layer,1] * (1. - Por)) + (1000. * Por) # Node density
            alpha = 0.
        elif node == self.number_nodes - 3: # Concrete
            cp = 880. # J/kg Heat Capacity
            # Den = 2300. # Node density
            Den = self.cden
            alpha = 0.
        elif node == self.number_nodes - 2: # Insulation
            cp = self.icp # J/kg Heat Capacity
            Den = self.iden # Node density
            alpha = 0.
        elif layer < self.top_wat:
            cp = 1005.
            Den = 1.23
            alpha = 0.
        else:
            cp = 4184.
            Den = 1000.
            alpha = self.alp
                
        if layer < self.top_wat and node == self.number_nodes -1:
            RL = self.RLA
        else:
            RL = self.RLW
        
        node_mass = RL * Den * np.pi * (self.Rx[node]**2 - self.Rx[node+1]**2)
    
        if node < self.number_nodes-1 or layer < self.top_wat:
            A =  (-Cxy - Cyz) / (node_mass * cp)
        else:
            if layer == self.top_wat or layer == self.number_layers - 1:
                A = (-Cxy - Cyz - ((mdotu + mdotd) * cp) - ((alpha * node_mass * cp) / RL**2)) / (node_mass * cp)
            else:
                A = (-Cxy - Cyz - ((mdotu + mdotd) * cp) - ((2. * alpha * node_mass * cp) / RL**2)) / (node_mass * cp)
            
        return A
    
    
    def coefficient_B(self, layer, node, Cxy):
        
        if node < self.number_nodes - 3: # Ground
            cp = self.Geo[layer,0] # J/kg Soil Heat Capacity
            Por = self.Geo[layer,2] # Soil porosity
            Den = (self.Geo[layer,1] * (1. - Por)) + (1000. * Por) # Node density
            # alpha = 0.
        elif node == self.number_nodes - 3: # Concrete
            cp = 880. # J/kg Heat Capacity
            # Den = 2300. # Node density
            Den = self.cden
            # alpha = 0.
        elif node == self.number_nodes - 2: # Insulation
            cp = self.icp # J/kg Heat Capacity
            Den = self.iden # Node density
            # alpha = 0.
        elif layer < self.top_wat:
            cp = 1005.
            Den = 1.23
            # alpha = 0.
        else:
            cp = 4184.
            Den = 1000.
            # alpha = self.alp
                
        if layer < self.top_wat and node == self.number_nodes - 1:
            RL = self.RLA
        else:
            RL = self.RLW
        
        node_mass = RL * Den * np.pi * (self.Rx[node]**2 - self.Rx[node+1]**2)
            
        B = Cxy / (node_mass * cp)
        
        return B
    
    def coefficient_C(self, layer, node, Cyz):
        
        if node < self.number_nodes - 3: # Ground
            cp = self.Geo[layer,0] # J/kg Soil Heat Capacity
            Por = self.Geo[layer,2] # Soil porosity
            Den = (self.Geo[layer,1] * (1. - Por)) + (1000. * Por) # Node density
            # alpha = 0.
        elif node == self.number_nodes - 3: # Concrete
            cp = 880. # J/kg Heat Capacity
            # Den = 2300. # Node density
            Den = self.cden
            # alpha = 0.
        elif node == self.number_nodes - 2: # Insulation
            cp = self.icp # J/kg Heat Capacity
            Den = self.iden # Node density
            # alpha = 0.
        elif layer < self.top_wat:
            cp = 1005.
            Den = 1.23
            # alpha = 0.
        else:
            cp = 4184.
            Den = 1000.
            # alpha = self.alp
                
        if layer == 0 and node == self.number_nodes - 1:
            RL = self.RLA
        else:
            RL = self.RLW
        
        node_mass = RL * Den * np.pi * (self.Rx[node]**2 - self.Rx[node+1]**2)
        
        C = Cyz / (node_mass * cp)
        
        return C
    
    def coefficient_D(self, layer, node, mdotd):
    
        if node == self.number_nodes - 1:
            cp = 4184.
            Den = 1000.
            
            if layer > 0:
                alpha = self.alp
            else:
                alpha = 0.
                
            node_mass = self.RLW * Den * np.pi * (self.Rx[node]**2 - self.Rx[node+1]**2)
        
            if layer < self.top_wat + 1:
                D = 0
            else:
                D = ((mdotd * cp) / (node_mass * cp)) + (alpha / self.RLW**2)
                
        else:
            D = 0
                
        return D
    
    def coefficient_E(self, layer, node, mdotu):
    
        if node == self.number_nodes - 1:
            cp = 4184.
            Den = 1000.
            if layer > 0:
                alpha = self.alp
            else:
                alpha = 0.
            
            node_mass = self.RLW * Den * np.pi * (self.Rx[node]**2 - self.Rx[node+1]**2)
        
            if layer < self.top_wat or layer == self.number_layers - 1:
                E = 0
            else:
                E = ((mdotu * cp) / (node_mass * cp)) + (alpha / self.RLW**2)
                
        else:
            E = 0
            
        return E
    
    def coefficient_F(self, layer, node): # Furthest ground node
        
        if layer == 0 and node == self.number_nodes - 1:
            RL = self.RLA
        else:
            RL = self.RLW
        
        layer1 = layer
        
        cp = self.Geo[layer,0] # J/kg Soil Heat Capacity
        Por = self.Geo[layer,2] # Soil porosity
        Den = (self.Geo[layer,1] * (1. - Por)) + (1000. * Por) # Node density
        
        node_mass = RL * Den * np.pi * (self.Rx[node]**2 - self.Rx[node+1]**2)
    
        if node == 0:
            rf2 = 0.5 * (self.Rx[node] + self.Rx[node+1])
            RTsx = np.log(self.Rx[node] / rf2)
            Cxy = ((2. * np.pi * self.XTC[layer1,node] * RL) / RTsx)
            
            F = (10. * Cxy) / (node_mass * cp)
        else:
            F = 0
    
        return F
    
    def coefficient_G(self, layer, node): # Axial heat transfer to lower layer
        
        if layer == 0 and node == self.number_nodes - 1:
            RL = self.RLA
        else:
            RL = self.RLW
        
        layer1 = layer
        layer2 = layer + 1
        
        if node < self.number_nodes - 3: # Ground
            cp = self.Geo[layer,0] # J/kg Soil Heat Capacity
            Por = self.Geo[layer,2] # Soil porosity
            Den = (self.Geo[layer,1] * (1. - Por)) + (1000. * Por) # Node density
            # alpha = 0.
        elif node == self.number_nodes - 3: # Concrete
            cp = 880. # J/kg Heat Capacity
            # Den = 2300. # Node density
            Den = self.cden
            # alpha = 0.
        else: # Insulation
            cp = self.icp # J/kg Heat Capacity
            Den = self.iden # Node density
            # alpha = 0.
        
        node_mass = RL * Den * np.pi * (self.Rx[node]**2 - self.Rx[node+1]**2)
    
        if layer < (self.number_layers - 1):
            G1 = (self.XTC[layer1,node] * np.pi * (self.Rx[node]**2 - self.Rx[node+1]**2)) / (RL / 2.)
            G2 = (self.XTC[layer2,node] * np.pi * (self.Rx[node]**2 - self.Rx[node+1]**2)) / (RL / 2.)            
            G = (G1 + G2) / (node_mass * cp)
        else:
            G = 0
    
        return G
    
    def coefficient_H(self, layer, node): # Axial heat transfer to upper layer
        
        if layer == 0 and node == self.number_nodes - 1:
            RL = self.RLA
        else:
            RL = self.RLW
        
        layer1 = layer - 1
        layer2 = layer
        
        if node < self.number_nodes - 3: # Ground
            cp = self.Geo[layer,0] # J/kg Soil Heat Capacity
            Por = self.Geo[layer,2] # Soil porosity
            Den = (self.Geo[layer,1] * (1. - Por)) + (1000. * Por) # Node density
            # alpha = 0.
        elif node == self.number_nodes - 3: # Concrete
            cp = 880. # J/kg Heat Capacity
            # Den = 2300. # Node density
            Den = self.cden
            # alpha = 0.
        else: # Insulation
            cp = self.icp # J/kg Heat Capacity
            Den = self.iden # Node density
            # alpha = 0.
        
        node_mass = RL * Den * np.pi * (self.Rx[node]**2 - self.Rx[node+1]**2)
    
        if layer > 0:
            H1 = (self.XTC[layer1,node] * np.pi * (self.Rx[node]**2 - self.Rx[node+1]**2)) / (RL / 2.)
            H2 = (self.XTC[layer2,node] * np.pi * (self.Rx[node]**2 - self.Rx[node+1]**2)) / (RL / 2.)            
            H = (H1 + H2) / (node_mass * cp)
        else:
            H = 0
    
        return H
    
    def coefficient_Z(self, layer, node, mdotu, mdotd, source_temp, return_temp): # Fluid movement between water layers (direct charge/discharge)
        
        if node == self.number_nodes - 1:
            cp = 4184.
            Den = 1000.
            
            node_mass = self.RLW * Den * np.pi * (self.Rx[node]**2 - self.Rx[node+1]**2)
        
            if layer < self.top_wat:
                Z = 0
            elif layer == self.top_wat:
                Z = (mdotd * cp * source_temp) / (node_mass * cp)
            elif layer == self.number_layers - 1:
                Z = (mdotu * cp * return_temp) / (node_mass * cp)
            else:
                Z = 0
                
        else:
            Z = 0
            
        return Z
    
    def set_of_coefficients(self, nodes_temp, mdotu, mdotd, source_temp, return_temp, Qp):
        
        H_Ms = np.zeros(self.number_layers)
    
        c = [[] for i in range(self.number_layers)]
        
        for layer in range(self.number_layers):
            
            for node in range(self.number_nodes):
                
                Cxy, Cyz, H_Ms[layer] = self.coefficient_base(layer, node, nodes_temp, Qp)
                
                coefficients = {'A': self.coefficient_A(layer, node, Cxy, Cyz, mdotu, mdotd),
                                'B': self.coefficient_B(layer, node, Cxy),
                                'C': self.coefficient_C(layer, node, Cyz),
                                'D': self.coefficient_D(layer, node, mdotd),
                                'E': self.coefficient_E(layer, node, mdotu),
                                'F': self.coefficient_F(layer, node),
                                'G': self.coefficient_G(layer, node),
                                'H': self.coefficient_H(layer, node),
                                'Z': self.coefficient_Z(layer, node, mdotu, mdotd, source_temp, return_temp),
                                }
                
                c[layer].append(coefficients)
        
        return c, H_Ms
    
    def new_nodes_temp(self, nodes_temp, source_temp, return_temp, charge, discharge, disDT, Qp, MTS):
    
        def model_temp(z, t, c):
            dzdt = []
            
            for layer in range(self.number_layers):
                            
                for lnode in range(self.number_nodes):
        
                    node = (layer * self.number_nodes) + lnode            
        
                    if lnode == 0: # Furthest earth 'ring'
                        Ti = nodes_temp[node]
                        Ti_b = nodes_temp[node + 1]
                        if layer < (self.number_layers - 1):
                            Tld = nodes_temp[node + self.number_nodes]
                        else:
                            Tld = 0.
                        if layer > 0:
                            Tlu = nodes_temp[node - self.number_nodes]
                        else:
                            Tlu = 0.
        
                        dTdt = MTS * ((c[layer][lnode]['A'] - c[layer][lnode]['G'] - c[layer][lnode]['H']) * Ti +
                                c[layer][lnode]['C'] * Ti_b +
                                c[layer][lnode]['G'] * Tld + 
                                c[layer][lnode]['H'] * Tlu + 
                                c[layer][lnode]['F'])
        
                        dzdt.append(dTdt)
        
                    elif lnode == (self.number_nodes - 1): # Store Fluid
                    
                        if layer == self.top_wat: # Top Water Layer
                            Ti = nodes_temp[node]
                            Ti_a = nodes_temp[node - 1]
                            Td = nodes_temp[node + self.number_nodes]
                            Tu = nodes_temp[self.number_nodes - 1]
                            Twa = nodes_temp[node - self.number_nodes - 1]
                            
                            cp = 4184.
                            Den = 1000.
                            
                            node_mass = self.RLW * Den * np.pi * self.Rx[self.number_nodes - 1]**2
                            
                            es = 10. * 0.61078 * np.exp((17.27*Ti)/(Ti+237.7)) # mbar
                            ea = 10. * 0.61078 * np.exp((17.27*Tu)/(Tu+237.7))
                            
                            Tsv = Ti / (1. - (0.378 * (es / 1013.25))) #mbar
                            Tav = Tu / (1. - (0.378 * (ea / 1013.25)))
                            diff = Tsv - Tav

                            if es != ea:
                                evap = max(0., ((((2.7 * np.cbrt(diff) * (es - ea)) * (1. + (0.61 * ((Td - Ti) / (es - ea))))) * np.pi * self.Rx[self.number_nodes - 1]**2) / (node_mass * cp)))
                            else:
                                evap = 0.
                            evap = 0.
                                
                            if Ti > Twa:
                                rad = (5.7e-8 * 0.9134 * (Ti - Twa)**4 * np.pi * self.Rx[self.number_nodes - 1]**2) / (node_mass * cp) # Radiation
                            else:
                                rad = 0.
                            rad = 0.

                            slw = (0.5 * (1./self.mu) * np.log(np.exp(0.) + np.exp(self.mu * (Td - Ti)))) # Slow buoyancy
                            slw = 0
                            
                            dTdt = MTS * (c[layer][lnode]['A'] * Ti +
                                          c[layer][lnode]['B'] * Ti_a + 
                                          c[layer][lnode]['E'] * Td +
                                          c[layer][lnode]['Z'] -
                                          rad -
                                          evap
                                          ) + slw
        
                            dzdt.append(dTdt)
                        
                        elif layer == self.number_layers - 1: # Bottom Water Layer
                            Ti = nodes_temp[node]
                            Ti_a = nodes_temp[node - 1]
                            Tu = nodes_temp[node - self.number_nodes]
                            
                            slw = (-0.5 * (1./self.mu) * np.log(np.exp(0.) + np.exp(self.mu * (Ti - Tu))))
                            slw = 0
            
                            dTdt = MTS * (c[layer][lnode]['A'] * Ti +
                                          c[layer][lnode]['B'] * Ti_a + 
                                          c[layer][lnode]['D'] * Tu +
                                          c[layer][lnode]['Z']
                                          ) + slw
                            
                            dzdt.append(dTdt)
                            
                        elif layer > self.top_wat and layer < self.number_layers - 1: # Intermediate Water Layer
                            Ti = nodes_temp[node]
                            Ti_a = nodes_temp[node - 1]
                            Tu = nodes_temp[node - self.number_nodes]
                            Td = nodes_temp[node + self.number_nodes]
                            
                            slw1 = (0.5 * (1./self.mu) * np.log(np.exp(0.) + np.exp(self.mu * (Td - Ti))))
                                        
                            slw2 = (-0.5 * (1./self.mu) * np.log(np.exp(0.) + np.exp(self.mu * (Ti - Tu))))
                                
                            slw = slw1 + slw2
                            slw = 0
            
                            dTdt = MTS * (c[layer][lnode]['A'] * Ti +
                                          c[layer][lnode]['B'] * Ti_a + 
                                          c[layer][lnode]['D'] * Tu +
                                          c[layer][lnode]['E'] * Td
                                          ) + slw
            
                            dzdt.append(dTdt)
                            
                        else:
                            if layer == 0:
                                Ti = nodes_temp[node]
                                Ti_a0 = nodes_temp[node - 1]
                                Ti_a1 = nodes_temp[node - 1 + self.number_nodes]
                                Ti_a2 = nodes_temp[node - 1 + 2 * self.number_nodes]
                                Ti_a3 = nodes_temp[node - 1 + 3 * self.number_nodes]
                                
                                Td = nodes_temp[node + 4 * self.number_nodes]
                                
                                cp = 1005.
                                Den = 1.23
                        
                                node_mass = self.RLA * Den * np.pi * 3.5**2
                                
                                es = 10. * 0.61078 * np.exp((17.27 * Td)/(Td + 237.7)) # mbar
                                ea = 10. * 0.61078 * np.exp((17.27 * Ti)/(Ti + 237.7))
                                
                                Tsv = Td / (1. - (0.378 * (es / 1013.25))) #mbar
                                Tav = Ti / (1. - (0.378 * (ea / 1013.25)))
                                diff = Tsv - Tav

                                if es != ea:    
                                    # evap = max(0., ((((2.7 * np.cbrt(diff) * (es - ea)) * (1. + (0.61 * ((Td - Ti) / (es - ea))))) * np.pi * self.Rx[9]**2) / (node_mass * cp)))
                                    evap = max(0., ((((2.7 * np.cbrt(diff) * (es - ea)) * (1. + (0.61 * ((Td - Ti) / (es - ea))))) * np.pi * self.Rx[10]**2) / (node_mass * cp)))
                                else:
                                    evap = 0.
                                evap = 0
                                    
                                dTdt = MTS * (c[layer][lnode]['A'] * Ti +
                                              c[layer+1][lnode]['A'] * Ti +
                                              c[layer+2][lnode]['A'] * Ti +
                                              c[layer+3][lnode]['A'] * Ti +  
                                              c[layer][lnode]['B'] * Ti_a0 +
                                              c[layer+1][lnode]['B'] * Ti_a1 +
                                              c[layer+2][lnode]['B'] * Ti_a2 +
                                              c[layer+3][lnode]['B'] * Ti_a3 +
                                              evap
                                              ) 
                
                                dzdt.append(dTdt)
                                
                            else:
                                dTdt = 0.
                                dzdt.append(dTdt)
        
                    else:
                        if layer >= self.top_wat or lnode != self.number_nodes - 2:
                            Ti = nodes_temp[node]
                            Ti_b = nodes_temp[node + 1]
                            Ti_a = nodes_temp[node - 1]
                            if layer < (self.number_layers - 1):
                                Tld = nodes_temp[node + self.number_nodes]
                            else:
                                Tld = 0.
                            if layer > 0:
                                Tlu = nodes_temp[node - self.number_nodes]
                            else:
                                Tlu = 0.
                            
                            dTdt = MTS * ((c[layer][lnode]['A'] - c[layer][lnode]['G'] - c[layer][lnode]['H']) * Ti +
                                    c[layer][lnode]['B'] * Ti_a +
                                    c[layer][lnode]['C'] * Ti_b +
                                    c[layer][lnode]['G'] * Tld + 
                                    c[layer][lnode]['H'] * Tlu + 
                                    c[layer][lnode]['F'])
            
                            dzdt.append(dTdt)
                            
                        else:
                            Ti = nodes_temp[node]
                            Ti_b = nodes_temp[self.number_nodes - 1] # consolidated air temp
                            Ti_a = nodes_temp[node - 1]
                            if layer < (self.number_layers - 1):
                                Tld = nodes_temp[node + self.number_nodes]
                            else:
                                Tld = 0.
                            if layer > 0:
                                Tlu = nodes_temp[node - self.number_nodes]
                            else:
                                Tlu = 0.
                            
                            dTdt = MTS * ((c[layer][lnode]['A'] - c[layer][lnode]['G'] - c[layer][lnode]['H']) * Ti +
                                    c[layer][lnode]['B'] * Ti_a +
                                    c[layer][lnode]['C'] * Ti_b +
                                    c[layer][lnode]['G'] * Tld + 
                                    c[layer][lnode]['H'] * Tlu + 
                                    c[layer][lnode]['F'])
            
                            dzdt.append(dTdt)
    
            return dzdt
    
        node_temp_list = []
        node_temp_list.append(nodes_temp)
        
        if charge > 0:
            mdotd = (charge * 3600.) / (MTS * 4.181 * (source_temp - nodes_temp[self.number_nodes * self.number_layers - 1])) # charge in kWh to kg/s, deltaT based on tank bottom temperature
        else:
            mdotd = 0.
            
        if discharge > 0:
            # mdotu = (discharge * 3600.) / (MTS * 4.181 * disDT) # discharge in kWh to kg/s, 10C deltaT
            mdotu = (discharge * 3600.) / (MTS * 4.181 * (nodes_temp[(self.top_wat + 1) * self.number_nodes - 1] - return_temp)) # discharge in kWh to kg/s
            # print(mdotu, discharge, return_temp, nodes_temp[(self.top_wat + 1) * self.number_nodes - 1] - return_temp)
        else:
            mdotu = 0.    
    
        # solve ODE
        # span for next time step
        tspan = [0,1]
        # solve for next step
        
        # new coefficients
        coefficients, H_Ms = self.set_of_coefficients(nodes_temp, mdotu, mdotd, source_temp, return_temp, Qp)
        
        z = odeint(
            model_temp, nodes_temp, tspan,
            args=(coefficients,))

        nodes_temp = z[1]
            
        node_temp_list.append(nodes_temp)
            
        # Calculate heat losses
        Hloss = 0.
        for layer in range(self.top_wat,self.number_layers):
            rf1 = 0.5 * (self.Rx[9] + self.Rx[10])
            RTs = np.log(rf1 / self.Rx[10])
            Cxy = 1. / ((1. / (H_Ms[layer] * 2. * np.pi * self.Rx[self.number_nodes - 1] * self.RLW)) + (RTs / (2. * np.pi * self.XTC[layer,self.number_nodes - 2] * self.RLW)))
            
            Hloss += Cxy * (nodes_temp[(layer + 1) * self.number_nodes - 1] - nodes_temp[(layer + 1) * self.number_nodes - 2]) * (MTS/3600.) * 0.001 # kWh
            
        # Heat transfer for heat transfer coefficient
        Qp = np.zeros(self.number_layers)
        for layer in range(self.number_layers):
            if layer == -1:
                RL = self.RLA
            else:
                RL = self.RLW
            
            nst = int((layer + 1) * self.number_nodes - 1)
            
            rf1 = 0.5 * (self.Rx[9] + self.Rx[self.number_nodes - 1])
            RTs = np.log(rf1 / self.Rx[self.number_nodes - 1])
            P_UAb = 1. / ((1./ (H_Ms[layer] * 2. * np.pi * self.Rx[self.number_nodes - 1] * RL)) + (RTs / (2. * np.pi * self.XTC[layer,self.number_nodes - 2] * RL)))
            Qp[layer] = max(1e-3, P_UAb * (nodes_temp[nst] - nodes_temp[nst-1]))  
    
        return node_temp_list, Hloss, Qp, Cxy, H_Ms[self.top_wat]