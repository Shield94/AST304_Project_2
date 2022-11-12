White Dwarf
===========

This directory contains scripts and a notebook to calculate the total mass and radius of a white dwarf considering its central pressure. eos.py is a script for calculating the pressure and density of a star. The Structure.py file then utilizes a differential function,initial conditions, and the fourth order runga-kutta method to integrate the initial pressure,radius and mass of a white dwarf. The final goal is to calculate the final mass, pressure and radius of white dwarfs with a mass between 0.1 and 1 solar masses.

Contents
--------

0. `README.md`: this file
1. `astro_const.py`: module containing physical constants, uses `astropy`
2. `eos.py`: starter code for the equation of state
3. `test_eos.py`: unit test for the equation of state; do not modify this file
4. `eos_table.txt`: comparison data for the equation of state, used by `test_eos.py`
5. `structure.py`: starter code to integrate stellar structure equations
6. `observations.py`: module containing class for reading and storing tabulated
    white dwarf masses, radii, and associated uncertainties.
7. `Joyce.txt`: Table 4, Joyce et al. (2018). Data for `observations.py`.
8. `group13_Project2_Final.ipynb`: jupyter notebook containing the code used to complete
    the integration, as well as a report with analysis of the model and a plot of the results.

Report/closeout
---------------

Before beginning the integration process to determine the mass-radius relationship for white dwarfs, we tested our equation of state functions for pressure and density that are written in eos.py and the test_eos.py returned a TRUE value, indicating that the EOS function we constructed was correct and all of our values were within the allowed tolerance.


To get the values of δ, η, and ξ, we had to estimate with our `integrate` function. The procedure to get these values is simple, but it was hard to determine how small these values should be. We settled with δ = 1e-4, η = 5e-7, and ξ = 0.05. These were the values that produced stable results when varying other parameters. These parameters also produced solid results for the table above, boosting our confidence in them. These values would prove useful when calculating the scalings for central pressure and density. As seen in **2.6**, P_guess does a horrible job, giving us a mass 3 times larger than what we started with! We had to scale back P_guess to get sensible consistent results. Using `bisect`, we were able to optimize P_guess and get back the same mass we handed `integrate`. Bisect also had arguments for an lower and upper bounds, these bounds were P_low and P_high, calulated in **2.6**.


Looking at the plot in the notebook, we see that, according to our model, radius decreases at a decreasing rate as the mass of the white dwarf increases. As you can see, our projected mass-radius relationship does not align very well with the observed data from HST and FUSE. It only passes through the error bars of two of the data points. This suggests our model of a white dwarf is not very realistic, and not a reliable predictor of the physical properties of dwarf stars. That being said, it does appear to have the same general shape as the data, so it could be used as a model of white dwarfs with a very low confidence, or if high precision was not needed.

Resources
---------

Coding standards and the grading rubric are posted on D2L.
