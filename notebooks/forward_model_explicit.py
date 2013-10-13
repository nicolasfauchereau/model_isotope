# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <markdowncell>

# # Forward eco-physiological modelling of $\Delta O_{18}$ in tree rings

# <markdowncell>

# forward version

# <markdowncell>

# ###TODO 
# 
# create a dictionnary with key : value = cell reference : parameter

# <codecell>

d_ex = {}
d_ex['E3']=gs

# <codecell>

def C2K(T):
    """
    conversion celsius to Kelvin
    """
    return T + 273.16

# <markdowncell>

# ### Inputs

# <codecell>

airtemp = 14.2666666667
windspeed = 5.4
rh = 84.4
pressure = 1015.3333333333
gs = 0.25
leaf_width = 0.015
PAR = 925

# <markdowncell>

# ### Energy balance calculations 

# <codecell>

rs = 1. / gs

# <codecell>

r_times_b = 3.8 * (leaf_width**0.25)*(windspeed**(-0.5))

# <codecell>

rb = 0.89 * r_times_b

# <codecell>

gr = (4*0.98*(0.000000056703)*(C2K(airtemp)**3))/(29.2)

# <codecell>

rBH = 1./((1./rb)+gr)

# <codecell>

Qtot = (PAR/4.6)*2

# <codecell>

Qabs = 0.5 * Qtot

# <markdowncell>

# ### Calculating $\epsilon$

# <codecell>

lesstemp = airtemp - 1

# <codecell>

estemp = (6.13753 * exp(lesstemp * ((18.564 - (lesstemp/254.4)))/(lesstemp +255.57)))*100.

# <codecell>

s = (((6.13753 * exp(airtemp * ((18.564 - (airtemp/254.4)))/(airtemp +255.57))))-estemp)/(C2K(airtemp)-C2K(lesstemp))

# <codecell>

smbar = 6.13753*(((airtemp+255.7)*(18.564 - (2*airtemp/254.4)) - airtemp*(18.564 - \
                    (airtemp/254.4)))/((airtemp+255.57)**2))*(exp(airtemp*(18.564 - (airtemp/254.4))/(airtemp + 255.57)))

# <codecell>

epsilon = (smbar*44012)/(29.2*(pressure))

# <codecell>

epsilon

# <markdowncell>

# ### Calculating $\frac{EA}{EI}$

# <codecell>

ea = (rh / 100) * (6.13753 * exp(airtemp * ((18.564 - (airtemp/254.4)))/(airtemp +255.57)))

# <codecell>

D = (((6.13753 * exp(airtemp * ((18.564 - (airtemp/254.4)))/(airtemp +255.57))))-ea)/pressure

# <codecell>

leaf_temp = airtemp + (rBH * ((Qabs * (rs + rb)) - (44012. * D))) / (29.2 * (rs + rb + (epsilon * rBH)))

# <codecell>

leaf_temp_K = C2K(leaf_temp)

# <codecell>

leaf_temp_K

# <markdowncell>

# ### Calculating transpiration

# <codecell>

transpiration = (epsilon * rBH * Qabs / 44012. + D) / (rs + rb + epsilon * rBH)

# <codecell>

transpiration

# <markdowncell>

# ### Craig / Gordon parameters 

# <codecell>


# <codecell>


