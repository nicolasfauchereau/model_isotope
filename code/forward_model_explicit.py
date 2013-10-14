# -*- coding: utf-8 -*-
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

inputs = pd.read_csv('../excel/inputs.csv', index_col = 0)

def C2K(T):
    return T + 273.16

gs = np.linspace(0.069, 0.35, 100)
leaf_width = [0.015,]
d_source_H2O = [-5.17,]
fract_through_stomata = [32,]
fract_through_boundary_layer = [28,]
# **eff_length** is the effective length at a minimum is double the average distance for a stoma to the central vein. For kauri, the avg distance is 0.0035m.  
# Thus, the minimum for this variable is about 0.0077m.  Untested theory however suggests at least 2 orders of magnitude for the avg distance,
# so 0.35m
eff_length = [0.0077,]
C = [5.55e4,]
# #### Constants for calculations of $\Delta$ cellulose and $\Delta$ leaf
C_O_fract = [27,]
Dcel_Dom = [9,]
prop_exc = [0.45,]
prop_Xylem = [0.56,]

### ===========================================================================
# ### Environmental inputs
airtemp = 14.2666666667
windspeed = 5.4
rh = 84.4
pressure = 1015.3333333333
PAR = 925

### ===========================================================================
### create the iterator defining the parameter space for the model 

parameter_space = product(gs,leaf_width,d_source_H2O,fract_through_stomata,fract_through_boundary_layer,eff_length,C,C_O_fract,Dcel_Dom,prop_exc,prop_Xylem)


for index in xrange(len(inputs))


param_outputs= [] 

### ===========================================================================
### start the loop over the parameter space 
for parameters in parameter_space: 
    gs,leaf_width,d_source_H2O,fract_through_stomata,fract_through_boundary_layer,eff_length,C,C_O_fract,Dcel_Dom,prop_exc,prop_Xylem = parameters

    # ### Energy balance calculations 
    rs = 1. / gs

    r_times_b = 3.8 * (leaf_width**0.25)*(windspeed**(-0.5))

    rb = 0.89 * r_times_b

    gr = (4*0.98*(0.000000056703)*(C2K(airtemp)**3))/(29.2)

    rBH = 1./((1./r_times_b)+gr)

    Qtot = (PAR/4.6)*2

    Qabs = 0.5 * Qtot

    # ### Calculating $\epsilon$

    lesstemp = airtemp - 1.

    estemp = (6.13753 * exp(lesstemp * ((18.564 - (lesstemp/254.4)))/(lesstemp +255.57)))*100

    lesstemp_K = C2K(lesstemp)

    s = (((6.13753 * exp(airtemp * ((18.564 - (airtemp/254.4)))/(airtemp +255.57))))-estemp)/(C2K(airtemp)-C2K(lesstemp))

    smbar = 6.13753*(((airtemp+255.7)*(18.564 - (2*airtemp/254.4)) - airtemp*(18.564 - \
                        (airtemp/254.4)))/((airtemp+255.57)**2))*(exp(airtemp*(18.564 - \
                                                    (airtemp/254.4))/(airtemp + 255.57)))

    epsilon = (smbar*44012)/(29.2*(pressure))

    # ### Calculating $\frac{EA}{EI}$

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

    # ### Calculating transpiration

    transpiration = (epsilon * rBH * Qabs / 44012. + D) \
    / (rs + rb + epsilon * rBH)

    # ### Craig / Gordon parameters 

    d_water_vapour = d_source_H2O + -1*(2.644-3.206*(1000/C2K(airtemp))+\
                                        1.534*(1000000/(C2K(airtemp)*C2K(airtemp))))

    ek = ((fract_through_stomata*1/gs)+(fract_through_boundary_layer*rb))/((1/gs)+rb)

    e_star = 2.644-3.206*((10**3)/leaf_temp_K)+1.534*((10**6)/(leaf_temp_K**2))

    dv = ((d_water_vapour/1000.)*(1+(d_source_H2O/1000))+(d_source_H2O/1000.))*1000.

    dv = ((d_water_vapour/1000)*(1+(d_source_H2O/1000))+(d_source_H2O/1000))*1000

    de = ek+e_star+((d_water_vapour-ek)*ea_ei)

    # ### Estimating the Peclet effect 

    D_Peclet = 0.000000119*(exp(-(637/(leaf_temp_K-137))))

    p_Peclet = (transpiration*eff_length)/(C*D_Peclet)

    DL = (de*(1-exp(-1*p_Peclet)))/p_Peclet

    dL = ((DL/1000)*(1+(d_source_H2O/1000))+(d_source_H2O/1000))*1000

    # ### Calculating $\Delta$ cellulose and $\Delta$ leaf

    D_sucrose = DL + C_O_fract

    D_cellulose = (DL*(1-(prop_exc*prop_Xylem)))+C_O_fract

    D_leaf = D_cellulose - Dcel_Dom

    d_sucrose = ((D_sucrose/1000)*(1+(d_source_H2O/1000))+(d_source_H2O/1000))*1000

    d_leaf = ((D_leaf/1000)*(1+(d_source_H2O/1000))+(d_source_H2O/1000))*1000

    # ### OUTPUT = $\Delta O_{18}$ in tree-rings cellulose

    OUTPUT = ((D_cellulose/1000)*(1+(d_source_H2O/1000))+(d_source_H2O/1000))*1000
    
    param_outputs.append(OUTPUT)

