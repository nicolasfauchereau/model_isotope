# Forward eco-physiological modelling of delta 18O in Kauri tree-rings 

This [IPython](www.ipython.org) notebook accompanies the paper:  
    
**Stable oxygen isotope signatures of early season wood in Agathis Australis (*Kauri*) tree-ring earlywood: Preliminary results and prospects for palaeoclimate reconstructions**

By *Andrew M. Lorrey, Tom H. Brookman, Michael N. Evans, Nicolas C. Fauchereau, Alison Criscitiello, Greg Eisheid, Anthony M. Fowler, Travis Horton and Daniel P. Schrag*

Submitted to *Journal of Geophysical Research - Biogeosciences*, December 2013.

---

## Installation instructions 

To read and run this IPython notebook, you need: 

+ Python, available at [www.pythopn.org]()

and the following libraries 

+ [IPython](http://ipython.org/)
+ [Numpy](http://www.numpy.org/)
+ [Scipy](http://www.scipy.org/)
+ [Pandas](http://pandas.pydata.org/)
+ [Matplotlib](http://matplotlib.org/)

All of these are freely available at the link provided and are generally available through the package manager of any modern Linux distribution.

However, we recommend using one of the excellent free, platform independent Scientific Python distributions, such as 

    + Anaconda, from Continuum Analytics: [https://store.continuum.io/cshop/anaconda/](https://store.continuum.io/cshop/anaconda/)
    + Canopy express, from Enthought: [https://www.enthought.com/downloads/](https://www.enthought.com/downloads/)
    + pyzo: [http://www.pyzo.org/](http://www.pyzo.org/)
    
Once you have installed one of these, please test the installation by running (the chevrons indicate a terminal emulator)

>> ipython notebook 

This should open the IPython notebook dashboard. Choose "New Notebook" and in the input cell, try importing the above libraries: 

import numpy
import scipy 
import pandas
import matplotlib 

If all of the above imports smoothly, you should be able to run the notebook provided as supplementary material to our paper.

---

## Running instructions 

The archive for the paper's supplementary material, available on github, is organized into 3 directories: 

+ **notebooks**: Contains the IPython notebook itself, named *forward_model_explicit.ipynb* 
+ **figures**: is where the figures are saved
+ **data**: should contain two csv files 
    1. *inputs.csv* are observed or estimated *environmental parameters* for the model, i.e. annual: 
     - Relative Humidity (%)
     - Air Temperature (degrees Celsius)
     - Pressure (hPa)
     - wind speed (m/s)
    Please follow the template provided in the csv files for the naming of the columns 

    2. *observed_tree_rings.csv* contains the annual *observed* Delta O18 values for Kauri tree-rings, it contains 2 columns: 
     - 'av': annual averages
     - 'std': annual standard deviation 
                
After unpacking the archive, cd to the notebook directory and run 

>> ipython notebook 

You should then be able to run the model, we recommend you run it sequentially, cell by cell, using the '|>' button in the notebook's toolbar

---

For all questions please contact me at: 

Nicolas.Fauchereau@gmail.com
or 
Nicolas.Fauchereau@niwa.co.nz
