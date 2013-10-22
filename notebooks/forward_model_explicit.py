# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <markdowncell>

# # Forward eco-physiological modelling of $\delta O_{18}$ in tree rings

# <markdowncell>

# ### Import the libraries we will need

# <codecell>

%pylab inline

# <codecell>

import os
from math import exp
from itertools import product
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt

# <markdowncell>

# ### Define the path where to find the csv files (from your home directory) 

# <codecell>

dpath = 'research/NIWA/paleo/model-isotope/excel'

dpath = os.path.join(os.environ['HOME'], dpath)

# <markdowncell>

# ### Define the path where to save the figures (from your home directory)

# <codecell>

fpath = 'research/NIWA/paleo/model-isotope/figures'

fpath = os.path.join(os.environ['HOME'], fpath)

# <markdowncell>

# ### Below is a small function to convert degrees Celsius to Kelvin 

# <codecell>

def C2K(T):
    """
    conversion celsius to Kelvin
    """
    return T + 273.16

# <markdowncell>

# ### Define the parameter space of the model

# <markdowncell>

# **Each parameter entering the calculations of $\delta 18_{O}$ can be defined as :**  
# > + `np.linspace(min, max, steps)` creates steps values between min (included) and max (not included)  
# > + `[value,]` A real value (e.g. 0.015), brackets around and comma **are important**  
# > + `[value1,value2,...] A list of specific real values to test 

# <codecell>

### ===========================================================================
### define the parameter space here 
#gs_l = [0.25,]
gs_l = np.linspace(0.069,0.35,100)
leaf_width_l= [0.015,]
d_source_H2O_l = [-5.17,]
fract_through_stomata_l = [32,]
fract_through_boundary_layer_l = [28,]
# **eff_length** is the effective length at a minimum is double the average distance for a stoma to the central vein. For kauri, the avg distance is 0.0035m.  
# Thus, the minimum for this variable is about 0.0077m.  Untested theory however suggests at least 2 orders of magnitude for the avg distance,
# so 0.35m
eff_length_l = [0.0077,]
C_l = [5.55e4,]
# #### Constants for calculations of $\Delta$ cellulose and $\Delta$ leaf
C_O_fract_l = [27,]
Dcel_Dom_l = [9,]
#prop_exc_l = np.linspace(0.25,0.65,100)
prop_exc_l = [0.45,]
#prop_Xylem_l = np.linspace(0.25,0.65,100)
prop_Xylem_l = [0.56,]
PAR_l = np.linspace(200,1000,100)
### ===========================================================================

# <markdowncell>

# ### Reads the inputs (relative humidity, air temperature, pressure, windspeed)

# <codecell>

inputs = pd.read_csv(os.path.join(dpath,'inputs.csv'), index_col = 0)

# <markdowncell>

# ### Reads the *observed* values of $\delta 18_{0}$ 

# <codecell>

obs = pd.read_csv(os.path.join(dpath,'observed_tree_rings.csv'), index_col=0)

# <markdowncell>

# ### Multiply by the standard deviation and add the mean to get the $\delta 18_{O}$ values

# <codecell>

obs['raw'] = (obs['av']*1.73) + 31.82

# <markdowncell>

# ### Below the main loops (over inputs and over parameter space) are implemented: *Do not* modify anything here

# <codecell>

### ===========================================================================

l = []

for index in xrange(len(inputs)): 

    ### ===========================================================================
    ### create the iterator defining the parameter space for the model 
    parameter_space = product(gs_l,leaf_width_l,d_source_H2O_l,fract_through_stomata_l,fract_through_boundary_layer_l,eff_length_l,C_l,C_O_fract_l,Dcel_Dom_l,prop_exc_l,prop_Xylem_l, PAR_l)

    rh, airtemp, pressure, windspeed = inputs.irow(index)

    param_outputs= [] 

    ### ===========================================================================
    ### start the loop over the parameter space 
    for parameters in parameter_space: 
        gs,leaf_width,d_source_H2O,fract_through_stomata,fract_through_boundary_layer,eff_length,C,C_O_fract,Dcel_Dom,prop_exc,prop_Xylem, PAR = parameters

        ### ===========================================================================
        ### Energy balance calculations 
        rs = 1. / gs

        r_times_b = 3.8 * (leaf_width**0.25)*(windspeed**(-0.5))

        rb = 0.89 * r_times_b

        gr = (4*0.98*(0.000000056703)*(C2K(airtemp)**3))/(29.2)

        rBH = 1./((1./r_times_b)+gr)

        Qtot = (PAR/4.6)*2

        Qabs = 0.5 * Qtot

        ### ===========================================================================
        ### Calculating $\epsilon$

        lesstemp = airtemp - 1.

        estemp = (6.13753 * exp(lesstemp * ((18.564 - (lesstemp/254.4)))/(lesstemp +255.57)))*100

        lesstemp_K = C2K(lesstemp)

        s = (((6.13753 * exp(airtemp * ((18.564 - (airtemp/254.4)))/(airtemp +255.57))))-estemp)/(C2K(airtemp)-C2K(lesstemp))

        smbar = 6.13753*(((airtemp+255.7)*(18.564 - (2*airtemp/254.4)) - airtemp*(18.564 - \
                            (airtemp/254.4)))/((airtemp+255.57)**2))*(exp(airtemp*(18.564 - \
                                                        (airtemp/254.4))/(airtemp + 255.57)))

        epsilon = (smbar*44012)/(29.2*(pressure))

        ### ===========================================================================
        ### Calculating $\frac{EA}{EI}$

        ea = (rh / 100) * (6.13753 * exp(airtemp * ((18.564 - (airtemp/254.4)))/(airtemp +255.57)))

        es = (6.13753 * exp(airtemp * ((18.564 - (airtemp/254.4)))/(airtemp +255.57)))

        D = (((6.13753 * exp(airtemp * ((18.564 - (airtemp/254.4)))\
                             /(airtemp +255.57))))-ea)/pressure

        temp_diff = (rBH*((Qabs*(rs+rb))-(44012*D)))/(29.2*(rs+rb+(epsilon*rBH)))

        leaf_temp = airtemp + temp_diff

        ei = (6.13753 * exp(leaf_temp * ((18.564 - (leaf_temp/254.4)))\
                            /(leaf_temp +255.57)))

        leaf_temp_K = C2K(leaf_temp)

        ea_ei = ea / ei

        ### ===========================================================================
        ### Calculating transpiration

        transpiration = (epsilon * rBH * Qabs / 44012. + D) \
        / (rs + rb + epsilon * rBH)

        ### ===========================================================================
        ### Craig & Gordon parameters 

        d_water_vapour = d_source_H2O + -1*(2.644-3.206*(1000/C2K(airtemp))+\
                                            1.534*(1e6/(C2K(airtemp)*C2K(airtemp))))

        ek = ((fract_through_stomata*1/gs)+(fract_through_boundary_layer*rb))/((1/gs)+rb)

        e_star = 2.644-3.206*((1e3)/leaf_temp_K)+1.534*((10**6)/(leaf_temp_K**2))

        dv = ((d_water_vapour/1e3)*(1+(d_source_H2O/1e3))+(d_source_H2O/1e3))*1e3

        dv = ((d_water_vapour/1e3)*(1+(d_source_H2O/1e3))+(d_source_H2O/1e3))*1e3

        de = ek+e_star+((d_water_vapour-ek)*ea_ei)

        ### ===========================================================================
        ### Estimating the Peclet effect 

        D_Peclet = 0.000000119*(exp(-(637/(leaf_temp_K-137))))

        p_Peclet = (transpiration*eff_length)/(C*D_Peclet)

        DL = (de*(1-exp(-1*p_Peclet)))/p_Peclet

        dL = ((DL/1e3)*(1+(d_source_H2O/1e3))+(d_source_H2O/1e3))*1e3

        ### ===========================================================================
        ### Calculating $\Delta$ cellulose and $\Delta$ leaf

        D_sucrose = DL + C_O_fract

        D_cellulose = (DL*(1-(prop_exc*prop_Xylem)))+C_O_fract

        D_leaf = D_cellulose - Dcel_Dom

        d_sucrose = ((D_sucrose/1e3)*(1+(d_source_H2O/1e3))+(d_source_H2O/1e3))*1e3

        d_leaf = ((D_leaf/1e3)*(1+(d_source_H2O/1e3))+(d_source_H2O/1000))*1e3

        ### ===========================================================================
        ### OUTPUT = $\Delta O_{18}$ in tree-rings cellulose

        OUTPUT = ((D_cellulose/1e3)*(1+(d_source_H2O/1e3))+(d_source_H2O/1e3))*1e3

        param_outputs.append(OUTPUT)

    l.append(param_outputs)

l = np.array(l)

# <markdowncell>

# ### Normalize (subtracts the average, divides by the sample standard deviation) the outputs 

# <codecell>

mean_l = l.mean(0)

std_l = l.std(0)

l_s = (l - mean_l) / std_l

# <markdowncell>

# ### Calculate correlation coefficients and find the max 

# <codecell>

R_coeff = []
for i in np.arange(l_s.shape[1]):                                                                                                                                                                   
    R_coeff.append(np.corrcoef(l_s[:,i], obs['av'].values.flatten())[0,1])
R_coeff = np.array(R_coeff)
maxR = np.argmax(R_coeff)

# <markdowncell>

# ### Now finds the parameters that minimize the root mean square error (${RMSE=\sqrt{\frac{1}{n}\sum_{t=1}^{n}{(observed_{t}-modelled_{t})^2}}}{}$)

# <codecell>

### ===========================================================================
RMSE = []
for i in np.arange(l.shape[1]): 
    RMSE.append(np.sqrt((obs['raw'] - l[:,i])**2).mean())                                                                                                                                                 
RMSE = np.array(RMSE)
minRMSE = np.argmin(RMSE)

# <markdowncell>

# ### Find the corresponding parameters in the parameter space

# <codecell>

parameter_space = product(gs_l,leaf_width_l,d_source_H2O_l,fract_through_stomata_l,fract_through_boundary_layer_l,eff_length_l,C_l,C_O_fract_l,Dcel_Dom_l,prop_exc_l,prop_Xylem_l, PAR_l)

parameter_space = np.array(list(parameter_space))

optimal_params_R = parameter_space[maxR]
optimal_params_RMSE = parameter_space[minRMSE]

# <markdowncell>

# ### Print the results

# <codecell>

results = """
The Maximum R value ({0:<4.2f}) is obtained for: 
gs = {1[0]:<5.2f}
leaf_width = {1[1]:<6.4f}
d_source_H2O = {1[2]:<5.2f}
fract_through_stomata = {1[3]:<5.2f}
fract_through_boundary_layer = {1[4]:<5.2f}
eff_length = {1[5]:<6.4f}
C = {1[6]:<5.2f}
C_O_fract = {1[7]:<5.2f}
Dcel_Dom = {1[8]:<5.2f}
prop_exc = {1[9]:<5.2f}
prop_Xylem = {1[10]:<5.2f}
PAR = {1[11]:<5.2f}
""".format(R_coeff[maxR],optimal_params_R)

print(results)         

# <codecell>

results = """
The Minimum RMSE value ({0:<4.2f}) is obtained for: 
gs = {1[0]:<5.2f}
leaf_width = {1[1]:<6.4f}
d_source_H2O = {1[2]:<5.2f}
fract_through_stomata = {1[3]:<5.2f}
fract_through_boundary_layer = {1[4]:<5.2f}
eff_length = {1[5]:<6.4f}
C = {1[6]:<5.2f}
C_O_fract = {1[7]:<5.2f}
Dcel_Dom = {1[8]:<5.2f}
prop_exc = {1[9]:<5.2f}
prop_Xylem = {1[10]:<5.2f}
PAR = {1[11]:<5.2f}
""".format(RMSE[minRMSE],optimal_params_RMSE)

print(results)         

# <markdowncell>

# ### Creates the figures

# <codecell>

### ===========================================================================
### raw modelled values (possibly P-dimensional)
f, ax = plt.subplots(figsize=(12,8))
#ax.plot(inputs.index, l, color='coral', lw=1.5, zorder=2)
ax.fill_between(inputs.index, l.min(1),l.max(1), color='coral', lw=1.5, zorder=2)

ax.plot(inputs.index, l[:,minRMSE], color='k', lw=3, label='best model (RMSE)')
ax.plot(inputs.index, l[:,maxR], color='g', lw=3, label='best model (R)')
ax.plot(obs.index, obs['raw'].values, color='steelblue', lw=2, label='observations')
ax.errorbar(obs.index, obs['raw'].values, yerr=obs['std'].values, fmt='o', color='steelblue')
ax.set_title('Raw data')
ax.legend(loc=0)
ax.grid('on')
f.savefig(os.path.join(fpath,'raw_modelled_delta18O.png'), bbox_inches='tight', dpi=200)
plt.show()

### ===========================================================================
### normalized modelled values (possibly P-dimensional) and observed values 
f, ax = plt.subplots(figsize=(12,8))
#ax.plot(inputs.index, l_s, color='coral', lw=1.5)
ax.fill_between(inputs.index, l_s.min(1),l_s.max(1), color='coral', zorder=2)
ax.plot(obs.index, obs['av'].values, color='steelblue', lw=2, label='observations')
ax.plot(inputs.index, l_s[:,maxR], color='k', lw=2, label='best model (R)')
ax.plot(inputs.index, l_s[:,minRMSE], color='g', lw=2, label='best model (RMSE)')
ax.errorbar(obs.index, obs['av'].values, yerr=obs['std'].values, fmt='o', color='steelblue')
ax.legend(loc=0)
ax.set_title('normalized data:\n modelled and observed $\delta O_{18}$', fontsize=14)
#ax.text(2006,2.5,'R=%4.2f' % (np.corrcoef(l_s.flatten(),obs['av'].values)[0,1]))
f.savefig(os.path.join(fpath,'normalized_modelled_delta18O.png'), bbox_inches='tight', dpi=200)
plt.show()

# <codecell>

l.shape

# <codecell>

l.min(0).shape

# <codecell>


