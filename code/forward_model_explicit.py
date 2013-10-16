#!/usr/bin/python                                                                                  
# -*- coding: utf-8 -*-
# ==============================================================================
# code forward_model_explicit.py
# Description:  purpose
# TAGS: tags: tags: tags
# created on 2013-10-15
# Nicolas Fauchereau <Nicolas.Fauchereau@gmail.com>
# ==============================================================================

from math import exp
from itertools import product
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt

def C2K(T):
    return T + 273.16

### ===========================================================================
### define the parameter space here 
gs_l = [0.25,]
#gs_l = np.linspace(0.069,0.35,100)
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

### ===========================================================================
inputs = pd.read_csv('../excel/inputs.csv', index_col = 0)

### ===========================================================================
#PAR_l = [925.,]
PAR_l = np.linspace(800,1000,100)
### ===========================================================================

l = []

for index in xrange(len(inputs)): 

    ### ===========================================================================
    ### create the iterator defining the parameter space for the model 
    parameter_space = product(gs_l,leaf_width_l,d_source_H2O_l,fract_through_stomata_l,fract_through_boundary_layer_l,eff_length_l,C_l,C_O_fract_l,Dcel_Dom_l,prop_exc_l,prop_Xylem_l, PAR_l)

    rh, airtemp, pressure, windspeed = inputs.iloc[index,:]

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
                                            1.534*(1000000/(C2K(airtemp)*C2K(airtemp))))

        ek = ((fract_through_stomata*1/gs)+(fract_through_boundary_layer*rb))/((1/gs)+rb)

        e_star = 2.644-3.206*((10**3)/leaf_temp_K)+1.534*((10**6)/(leaf_temp_K**2))

        dv = ((d_water_vapour/1000.)*(1+(d_source_H2O/1000))+(d_source_H2O/1000.))*1000.

        dv = ((d_water_vapour/1000)*(1+(d_source_H2O/1000))+(d_source_H2O/1000))*1000

        de = ek+e_star+((d_water_vapour-ek)*ea_ei)

        ### ===========================================================================
        ### Estimating the Peclet effect 

        D_Peclet = 0.000000119*(exp(-(637/(leaf_temp_K-137))))

        p_Peclet = (transpiration*eff_length)/(C*D_Peclet)

        DL = (de*(1-exp(-1*p_Peclet)))/p_Peclet

        dL = ((DL/1000)*(1+(d_source_H2O/1000))+(d_source_H2O/1000))*1000

        ### ===========================================================================
        ### Calculating $\Delta$ cellulose and $\Delta$ leaf

        D_sucrose = DL + C_O_fract

        D_cellulose = (DL*(1-(prop_exc*prop_Xylem)))+C_O_fract

        D_leaf = D_cellulose - Dcel_Dom

        d_sucrose = ((D_sucrose/1000)*(1+(d_source_H2O/1000))+(d_source_H2O/1000))*1000

        d_leaf = ((D_leaf/1000)*(1+(d_source_H2O/1000))+(d_source_H2O/1000))*1000

        ### ===========================================================================
        ### OUTPUT = $\Delta O_{18}$ in tree-rings cellulose

        OUTPUT = ((D_cellulose/1000)*(1+(d_source_H2O/1000))+(d_source_H2O/1000))*1000

        param_outputs.append(OUTPUT)

    l.append(param_outputs)

l = np.array(l)

mean_l = l.mean(0)
std_l = l.std(0)

l_s = (l - mean_l) / std_l

### ===========================================================================
### load the OBSERVED delta 18O values 

obs = pd.read_csv('../excel/observed_tree_rings.csv', index_col=0)

### ===========================================================================
### plots 
f, ax = plt.subplots(figsize=(10,6))
ax.plot(inputs.index, l, color='r', lw=1.5)
ax.set_title('Raw data')
f.savefig('../figures/raw_modelled_delta18O.png', bbox_inches='tight', dpi=200)
plt.close(f)

f, ax = plt.subplots(figsize=(10,6))
ax.plot(inputs.index, l_s, color='coral', lw=1.5, label='model')
ax.plot(obs.index, obs['av'].values, color='steelblue', lw=2, label='observations')
ax.errorbar(obs.index, obs['av'].values, yerr=obs['std'].values, fmt='o', color='steelblue')
ax.legend(loc=0)
ax.set_title('normalized data:\n modelled and observed $\delta O_{18}$', fontsize=14)
ax.text(2006,2.5,'R=%4.2f' % (np.corrcoef(l_s.flatten(),obs['av'].values)[0,1]))
f.savefig('../figures/normalized_modelled_delta18O.png', bbox_inches='tight', dpi=200)
plt.close(f)

### ===========================================================================
### calculate correlation coefficients 
R_coeff = []

for i in np.arange(l.shape[1]): 
    R_coeff.append(np.corrcoef(l_s[:,i], obs['av'].values.flatten())[0,1])
