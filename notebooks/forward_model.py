# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <markdowncell>

# # Forward eco-physiological modelling of $\Delta O_{18}$ in tree rings

# <markdowncell>

# I am working backward usin the excel spreadsheet as template

# <codecell>


# <codecell>


# <codecell>


# <markdowncell>

# ### first level 

# <codecell>

B16 = 2.644-3.206*((10^3)/K10)+1.534*((10^6)/(K10^2))

# <codecell>

E30 = (E23*(1-(B30*B31)))+B28

# <markdowncell>

# ### output

# <codecell>

output = ((E30/1000)*(1+(B16/1000))+(B16/1000))*1000

