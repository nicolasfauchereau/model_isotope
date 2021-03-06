# Mechanistic modeling of delta 18O in Kauri tree-rings

This repository contains the supplementary material accompanying the paper:

**Stable oxygen isotope signatures of early season wood in Agathis Australis (*Kauri*) tree-ring earlywood: Prospects for palaeoclimate reconstructions**

By *Andrew M. Lorrey, Tom H. Brookman, Michael N. Evans, Nicolas C. Fauchereau, Margaret Barbour, Cate Macinnis-Ng, Alison Criscitiello, Greg Eisheid, Anthony M. Fowler, Travis Horton and Daniel P. Schrag*

Submitted to *Dendrochronologia*, 2015.

It basically consists of an [IPython](www.ipython.org) notebook and a couple of csv files containing environmental parameters (Relative Humidity, Air Temperature, Atmospheric Pressure, Wind speed) as well as *observed* (derived by isotopic techniques) delta 18O values for Kauri (Agathis Australis) tree-rings.

The code for this model originates from an Excel file provided by Dr. Margaret Barbour, and follows the BFRE04 model, augmented with a leaf energy balance model

All code and data is available on github at:

[https://github.com/nicolasfauchereau/model_isotope](https://github.com/nicolasfauchereau/model_isotope)

See diagram below, which depicts a flow chart related to Figure 3 in the manuscript. The calculations associated with each major part of the flow chart are found in Supplementary Workbook 1.

<img src='https://raw.githubusercontent.com/nicolasfauchereau/model_isotope/master/figures/Lorreyetal2015SF.jpg' width=800>

IPython notebooks are now rendered as HTML on github, follow [this link](https://github.com/nicolasfauchereau/model_isotope/blob/master/notebooks/forward_model_explicit_widgetsop.ipynb) for a static version.

---

## Installation instructions

To read and run this IPython notebook, you need:

+ Python, available at [www.python.org]()

and the following libraries

+ [IPython](http://ipython.org/)
+ [Numpy](http://www.numpy.org/)
+ [Scipy](http://www.scipy.org/)
+ [Pandas](http://pandas.pydata.org/)
+ [Matplotlib](http://matplotlib.org/)
+ [statsmodels](http://statsmodels.sourceforge.net/)

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
    import statsmodels

If all of the above imports smoothly, you should be able to run the notebook provided as supplementary material to our paper.

---

## Running instructions

The archive for the paper's supplementary material, available on github, is organized into 3 directories:

+ **notebooks**: Contains the IPython notebook itself, named *forward_model_explicit.ipynb*
+ **figures**: is where the figures are saved
+ **excel**: should contain two csv files
    1. *inputs.csv* are observed or estimated *environmental parameters* for the model, i.e. annual:
     - Relative Humidity (%)
     - Air Temperature (degrees Celsius)
     - Pressure (hPa)
     - wind speed (m/s)

    2. *observed_tree_rings.csv* contains the annual *observed* Delta O18 values for Kauri tree-rings, it contains 2 columns:
     - 'av': annual averages
     - 'std': annual standard deviation

Please follow the template provided in the csv files for the naming of the columns.

After unpacking the archive, cd into the *notebooks* directory and run:

    >> ipython notebook

You should then be able to run the model, we recommend you run it sequentially, cell by cell, using the '|>' button in the notebook's toolbar

---

For all questions regarding the code itself or the installation of the required libraries, please contact Dr. Nicolas Fauchereau at:

+ [Nicolas.Fauchereau@gmail.com](mailto:Nicolas.Fauchereau@gmail.com)
+ [Nicolas.Fauchereau@niwa.co.nz](mailto:Nicolas.Fauchereau@niwa.co.nz)

For all questions regarding the model, data and methodology, please contact Dr. Andrew Lorrey at:

+ [Andrew.Lorrey@niwa.co.nz](mailto:Andrew.Lorrey@niwa.co.nz)
