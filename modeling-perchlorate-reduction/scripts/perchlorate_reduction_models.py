"""
THEORETICAL YIELD AND STOICHIOMETRY

Derived from reduction potential as in Rittmann and McCarty, 2001.

Reactions for acceptor, donor, and cell synthesis are balanced from energetics using the equation:

R = (fe x Ra) + (fs x Rc) - Rd

For modeling, the important values to obtain are the the mol acetate / mol electron acceptor and g cells / mol electron acceptor (biomass: 113 g / mol C5H7O2N)

Perchlorate to chlorate (ClO4- > ClO3-)
- Ra: ½ ClO4- + H+ + e- > ½ ClO3- + ½ H2O
- Rd: 1/4 HCO3- + 9/8 H+ + e- = 1/8 CH3COO- + 1/2 H2O
- Rc: 1/4 HCO3- + 1/20 NH4+ + 6/5 H+ + e- = 1/20 C5H7O2N + 13/20 H2O

Chlorate to chloride (ClO3- > Cl-)
- Ra: 1/6 ClO3- + H+ + e- → 1/6 Cl- + 3/6 H2O  
- Rd: 1/4 HCO3- + 9/8 H+ + e- = 1/8 CH3COO- + 1/2 H2O
- Rc: 1/4 HCO3- + 1/20 NH4+ + 6/5 H+ + e- = 1/20 C5H7O2N + 13/20 H2O
"""

from energetics import *
import pandas as pd

# Perchlorate reduction to chlorate

# Input
ne_perc = 2 # number of electrons
E_perc = 0.813 # V, redox potential
fs = energy_to_fractions(redox_to_Ga(ne_perc,E_perc))[2] 
fe = energy_to_fractions(redox_to_Ga(ne_perc,E_perc))[3] 
stoich_clo4 = pd.DataFrame({"molecule" : ["e-","H+","H2O","ClO4-","ClO3-","HCO3-","CH3COO-","NH4+","C5H7O2N"],
                            "Ra" : [fe * x for x in [1, 1, 0.5, 0.5, 0.5, 0, 0, 0, 0]],
                            "Rd" : [1, 1.125, 0.5, 0, 0, 0.25, 0.125, 0, 0],
                            "Rc" : [fs * x for x in [1, 1.2, 0.65, 0, 0, 0.25, 0, 0.05, 0.05]]})


stoich_clo4 = stoich_clo4.set_index('molecule')
stoich_clo4_acet = stoich_clo4.loc['CH3COO-','Rd'] / stoich_clo4.loc['ClO4-','Ra'] # mol acetate / mol perchlorate
yield_clo4_acet = stoich_clo4.loc['C5H7O2N','Rc']  / stoich_clo4.loc['ClO4-','Ra'] # mol biomass / mol perchlorate
ypc = [yield_clo4_acet, 0.0] # Yield coefficients for PRB and CRB
ypa = [stoich_clo4_acet, 0.0] # Stoichiometric ratios for PRB and CRB

# Chlorate reduction to chloride

# Input
E_chlor = 0.792 # V, redox potential
ne_chlor = 6 # number of electrons
fs = energy_to_fractions(redox_to_Ga(ne_chlor,E_chlor))[2]
fe = energy_to_fractions(redox_to_Ga(ne_chlor,E_chlor))[3]
stoich_clo3 = pd.DataFrame({"molecule" : ["e-","H+","H2O","ClO3-","Cl-","HCO3-","CH3COO-","NH4+","C5H7O2N"],
                        "Ra" : [fe * x for x in [1, 1,  0.5, 0.16667, 0.16667, 0, 0, 0, 0]],
                        "Rd" : [1, 1.125, 0.5, 0, 0, 0.25, 0.125, 0, 0],
                        "Rc" : [fs * x for x in [1, 1.2, 0.65, 0, 0, 0.25, 0, 0.05, 0.05]]})

# Calculation of yield coefficients and stoichiometric ratios
stoich_clo3 = stoich_clo3.set_index('molecule')
stoich_clo3_acet = stoich_clo3.loc['CH3COO-','Rd'] / stoich_clo3.loc['ClO3-','Ra'] # mol acetate / mol chlorate
yield_clo3_acet = stoich_clo3.loc['C5H7O2N','Rc']  / stoich_clo3.loc['ClO3-','Ra'] # mol biomass / mol chlorate
yc = [yield_clo3_acet, yield_clo3_acet] # Yield coefficients for PRB and CRB
yca = [stoich_clo3_acet, stoich_clo3_acet] # Stoichiometric ratios for PRB and CRB


"""
KINETIC PARAMETERS

"""

import numpy as np

# Substrate affinity [PRB,CRB], Ks (M) 
ks_clo4 = np.array([0.0060,10000])*10**-3 # Ks perchlorate, 0.0019 (Dudley) super high value for CRB because it does not catalyze
ks_clo3 = np.array([0.0074,0.0074])*10**-3 # Ks chorate  0.0007 M - Dudley #0.159 CRB
ks_acet = np.array([1,1])*10**-3 # Ks acetate

# Substrate inhibition factor, Haldane kinetics [PRB,CRB], Ki (M)
# Inhibition only occurs for PRB with perchlorate
# Youngblut et. al 2016 http://www.jbc.org/content/291/17/9190.full.pdf
ki_clo4 = np.array([10000, 10000])*10**-3   # Ki perchlorate, only low for ClO4- and PRB
ki_clo3 = np.array([10000, 10000])*10**-3   # Ki chlorate
ki_acet = np.array([10000, 10000])*10**-3   # Ki acetate

# Growth and death rate
mu = np.array([0.5,0.5]) # Maximum growth rate
m = np.array([0.0,0.0])# Death rate


"""
Equations tailored to modeling PRB and CRB

"""

from kinetics import *
from scipy.integrate import odeint

# Set indexes of state variables
id_prb = 0 # PRB, perchlorate-reducing bacteria
id_crb = 1 # CRB, chlorate-reducing bacteria
id_acet = 2 # Acetate
id_clo4 = 3 # Perchlorate (ClO4-)
id_clo3 = 4 # Chlorate (ClO3-)
id_prb_clo3 = 5 # Chlorate consumed by PRB
id_crb_clo3 = 6 # Chlorate consumed by CRB

# Name values
time = 'Time (h)'
values = ["PRB","CRB","C2H3O2-","ClO4-","ClO3-", "ClO3- to PRB", "ClO3- to CRB", time, "CRB/PRB", "Total Cells"]   

def mm_kinetics(x,t): # Concentrations, Time
    
    # Empty list for reactions
    rxn1_acc = [0] * len([id_prb,id_crb])
    rxn2_acc = [0] * len([id_prb,id_crb])
    rxn_don = [0] * len([id_prb,id_crb])
  
    # Define terms for reactions
    # PRB
    rxn1_acc[id_prb] = mm(x[id_clo4], ks_clo4[id_prb]) # ClO4- > ClO3-
    rxn2_acc[id_prb] = mm(x[id_clo3], ks_clo3[id_prb]) # ClO3- > Cl-
    rxn_don[id_prb] = mm(x[id_acet], ks_acet[id_prb]) # acetate > co2
    
    # CRB
    rxn2_acc[id_crb] = mm(x[id_clo3],ks_clo3[id_crb]) # ClO3- > Cl-
    rxn_don[id_crb] = mm(x[id_acet], ks_acet[id_crb]) # acetate > co2
    
    dxdt = [0] * 7 # Empty list

    # Michaelis-menton changes in [PRB] and [CRB]
    dxdt[id_prb] =  ypc[id_prb] * mu[id_prb] * x[id_prb] * rxn1_acc[id_prb] * rxn_don[id_prb] + yc[id_prb] * mu[id_prb] * x[id_prb] * rxn2_acc[id_prb] * rxn_don[id_prb] - m[id_prb] * x[id_prb]
    dxdt[id_crb] =  yc[id_crb] * mu[id_crb] * x[id_crb] * rxn2_acc[id_crb] * rxn_don[id_crb]
    
    # Michaelis-menton changes in [acetate], [perchlorate], and [chlorate]
    dxdt[id_acet] = - ypa[id_prb] * mu[id_prb] * x[id_prb] * rxn1_acc[id_prb] * rxn_don[id_prb] - mu[id_prb] * x[id_prb] * rxn2_acc[id_prb] * rxn_don[id_prb] - mu[id_crb] * x[id_crb] * rxn2_acc[id_crb] * rxn_don[id_crb]
    dxdt[id_clo4] = - mu[id_prb] * x[id_prb] * rxn1_acc[id_prb] * rxn_don[id_prb]
    dxdt[id_clo3] = mu[id_prb] * x[id_prb] * rxn1_acc[id_prb] * rxn_don[id_prb] - mu[id_prb] * x[id_prb] * rxn2_acc[id_prb] * rxn_don[id_prb] - mu[id_crb] * x[id_crb] * rxn2_acc[id_crb] * rxn_don[id_crb]

    # Calculate ClO3- consumed by PRB or CRB at each timestep
    dxdt[id_prb_clo3] = mu[id_prb] * x[id_prb] * rxn2_acc[id_prb] * rxn_don[id_prb]
    dxdt[id_crb_clo3] = mu[id_crb] * x[id_crb] * rxn2_acc[id_crb] * rxn_don[id_crb]

    return dxdt

def ci_kinetics(x,t): # Concentrations, Time
    
    # Empty list for reacts
    rxn1_acc = [0] * len([id_prb,id_crb])
    rxn2_acc = [0] * len([id_prb,id_crb])
    rxn_don = [0] * len([id_prb,id_crb])
  
    rxn1_acc[id_prb] = ci([x[id_clo4], x[id_clo3]], [ks_clo4[id_prb], ks_clo3[id_prb]]) # ClO4- > ClO3-
    rxn2_acc[id_prb] = ci([x[id_clo3], x[id_clo4]], [ks_clo3[id_prb], ks_clo4[id_prb]]) # ClO3- > Cl-
    rxn_don[id_prb] = ci(x[id_acet], ks_acet[id_prb]) # acetate > co2
    
    rxn2_acc[id_crb] = ci(x[id_clo3], ks_clo3[id_crb]) # ClO3- > Cl-
    rxn_don[id_crb] = ci(x[id_acet], ks_acet[id_crb]) # acetate > co2
    
    dxdt = [0] * 7 # Empty list

    # Michaelis-menton changes in [PRB] and [CRB]
    dxdt[id_prb] =  ypc[id_prb] * mu[id_prb] * x[id_prb] * rxn1_acc[id_prb] * rxn_don[id_prb] + yc[id_prb] *mu[id_prb] * x[id_prb] * rxn2_acc[id_prb] * rxn_don[id_prb] - m[id_prb] * x[id_prb]
    dxdt[id_crb] =  yc[id_crb] * mu[id_crb] * x[id_crb] * rxn2_acc[id_crb] * rxn_don[id_crb]
    
    # Michaelis-menton changes in [acetate], [perchlorate], and [chlorate]
    dxdt[id_acet] = - ypa[id_prb] * mu[id_prb] * x[id_prb] * rxn1_acc[id_prb] * rxn_don[id_prb] - yca[id_prb] * mu[id_prb] * x[id_prb] * rxn2_acc[id_prb] * rxn_don[id_prb] - yca[id_crb] * mu[id_crb] * x[id_crb] * rxn2_acc[id_crb] * rxn_don[id_crb]
    dxdt[id_clo4] = - mu[id_prb] * x[id_prb] * rxn1_acc[id_prb] * rxn_don[id_prb]
    dxdt[id_clo3] = mu[id_prb] * x[id_prb] * rxn1_acc[id_prb] * rxn_don[id_prb] - mu[id_prb] * x[id_prb] * rxn2_acc[id_prb] * rxn_don[id_prb] - mu[id_crb] * x[id_crb] * rxn2_acc[id_crb] * rxn_don[id_crb]
    
    # Calculate ClO3- consumed by PRB or CRB at each timestep
    dxdt[id_prb_clo3] = mu[id_prb] * x[id_prb] * rxn2_acc[id_prb] * rxn_don[id_prb]
    dxdt[id_crb_clo3] = mu[id_crb] * x[id_crb] * rxn2_acc[id_crb] * rxn_don[id_crb]

    return dxdt

# Approximation of Michaelis-Menton kinetics using ECA
# Tang and Riley 2013

# ECA: individual terms for model
def eca(ss,ki,ks,kkd,eed):

    # ss = Substrate concentration, perchlorate and/or chlorate
    # ki = Inhition factor, Haldane kinetics, perchlorate and/or chlorate
    # ks = Substrate affinity, perchlorate and/or chlorate
    # kkd = Substrate affinity, PRB and/or CRB
    # eed = Cell concentration, PRB and/or CRB

    # Allow lists to be one item
    if type(ss) is not list: ss = [ss] 
    if type(ki) is not list: ki = [ki]
    if type(ks) is not list: ks = [ks] 
    if type(kkd) is not list: kkd = [kkd] 
    if type(eed) is not list: eed = [eed] 
    
    # Creates one or two terms for [Substrate] / Ko
    if len(ss) > 1:
        dmn1 = ss[0] / (ks[0] / (ss[0]/ki[0] + 1)) # ClO4- or ClO3- [Substrate] / Ko
        dmn2 = ss[1] / (ks[1] / (ss[1]/ki[1] + 1)) # ClO4- or ClO3- [Substrate] / Ko
    else:
        dmn1 = ss[0] / (ks[0] / (ss[0]/ki[0] + 1)) #  ClO4- or ClO3- [Substrate] / Ko
        dmn2 = 0 # No second substrate
         
    # Creates one or two terms for [Biomass] / Ks_substrate
    if len(kkd) > 1:
        dmn3 = eed[0] / kkd[0] # PRB or CRB [Biomass] / Ks_substrate
        dmn4 = eed[1] / kkd[1] # PRB or CRB [Biomass] / Ks_substrate
    else:
        dmn3 = eed[0] / kkd[0] # PRB or CRB [Biomass] / Ks_substrate
        dmn4 = 0 # No second organism
    
    equation = ss[0] / (ks[0] * (1 + dmn1 + dmn2 + dmn3 + dmn4))
    return equation

# ECA: ODE model
def eca_kinetics(x,t): # Concentrations, Time
    
    # Examples
    #dxdt(id_ac)   = - mu(1) * x(id_prb) * edonor2(1) * eacceptor2(1);    
    #dxdt(id_clo4) = - mu(1) * x(id_prb) * edonor1(1) * eacceptor1(1);    
    #dxdt(id_clo3) = mu(1) * x(id_prb) * edonor1(1) * eacceptor1(1) - mu(1) * x(id_prb) * edonor2(1) * eacceptor2(1);
    #dxdt(id_prb) = mu(1) * x(id_prb) * edonor2(1) * eacceptor2(1)  - m(1) * x(id_prb);
    #dxdt(id_crb) = 0;
    
    rxn1_acc = [0] * len([id_prb,id_crb])
    rxn1_don = [0] * len([id_prb,id_crb])
    rxn2_acc = [0] * len([id_prb,id_crb])
    rxn2_don = [0] * len([id_prb,id_crb])
        
    # Reaction 1: Respiration of perchlorate (PRB ONLY)
    # Acceptor: ClO4- > ClO3-
    # Donor: Acetate > CO2 

    rxn1_acc[id_prb] = eca([x[id_clo4], x[id_clo3]], # PRB reduce both substrates
                           [ki_clo4[id_prb], ki_clo3[id_prb]], 
                           [ks_clo4[id_prb], ks_clo3[id_prb]], 
                           ks_clo4[id_prb], 
                           x[id_prb])
    
    rxn1_don[id_prb] = eca(x[id_acet], 
                           ki_acet[id_prb],
                           ks_acet[id_prb], 
                           [ks_acet[id_prb], ks_acet[id_prb]], 
                           [x[id_prb], x[id_crb]])
    
    # Reaction 2: Respiration of chlorate and on (PRB AND CRB)
    # Acceptor: ClO3- > Cl-
    # Donor: Acetate > CO2 

    rxn2_acc[id_prb] = eca([x[id_clo3], x[id_clo4]], # PRB reduce both substrates
                           [ki_clo3[id_prb], ki_clo4[id_prb]],
                           [ks_clo3[id_prb], ks_clo4[id_prb]], 
                           [ks_clo3[id_prb], ks_clo3[id_crb]], 
                           [x[id_prb], x[id_crb]])
    
    rxn2_don[id_prb] = eca(x[id_acet], 
                           ki_acet[id_prb],
                           ks_acet[id_prb], 
                           [ks_acet[id_prb], ks_acet[id_crb]], 
                           [x[id_prb], x[id_crb]])
    
    rxn2_acc[id_crb] = eca(x[id_clo3], # CRB reduce only ClO3-
                           ki_clo3[id_crb], 
                           ks_clo3[id_crb], 
                           [ks_clo3[id_crb], ks_clo3[id_prb]], 
                           [x[id_crb], x[id_prb]])
    
    rxn2_don[id_crb] = eca(x[id_acet], 
                           ki_acet[id_crb],
                           ks_acet[id_crb], 
                           [ks_acet[id_crb], ks_acet[id_prb]], 
                           [x[id_crb], x[id_prb]])

    dxdt = [0] * 7 # Empty list
    
    # Calculate change in concentration after timestep
    
    dxdt[id_acet] = (- ypa[id_prb] * mu[id_prb] * x[id_prb] * rxn1_acc[id_prb] * rxn1_don[id_prb] # Rxn 1
                     - yca[id_prb] * mu[id_prb] * x[id_prb] * rxn2_acc[id_prb] * rxn2_don[id_prb] # Rxn 2 PRB
                     - yca[id_crb] * mu[id_crb] * x[id_crb] * rxn2_acc[id_crb] * rxn2_don[id_crb]) # Rxn 2 CRB
    
    dxdt[id_clo4] = - mu[id_prb] * x[id_prb] * rxn1_acc[id_prb] * rxn1_don[id_prb] # Rxn 1
    
    dxdt[id_clo3] = (  mu[id_prb] * x[id_prb] * rxn1_acc[id_prb] * rxn1_don[id_prb] # From Rxn 1
                     - mu[id_prb] * x[id_prb] * rxn2_acc[id_prb] * rxn2_don[id_prb] # Rxn 2 PRB
                     - mu[id_crb] * x[id_crb] * rxn2_acc[id_crb] * rxn2_don[id_crb]) # Rxn 2 CRB
    
    dxdt[id_prb] =  (  ypc[id_prb] * mu[id_prb] * x[id_prb] * rxn1_acc[id_prb] * rxn1_don[id_prb] # Rxn 1
                     + yc[id_prb]  * mu[id_prb] * x[id_prb] * rxn2_acc[id_prb] * rxn2_don[id_prb] # Rxn 2
                     - m[id_prb] * x[id_prb]) # Death
    
    dxdt[id_crb] =  (  yc[id_crb] * mu[id_crb] * x[id_crb] * rxn2_acc[id_crb] * rxn2_don[id_crb] # Rxn 2
                     - m[id_crb] * x[id_crb]) # Death
    
    # Calculate ClO3- consumed by PRB or CRB at each timestep
    
    dxdt[id_prb_clo3] = mu[id_prb] * x[id_prb] * rxn2_acc[id_prb] * rxn2_don[id_prb]
    dxdt[id_crb_clo3] = mu[id_crb] * x[id_crb] * rxn2_acc[id_crb] * rxn2_don[id_crb]
    
    return dxdt


