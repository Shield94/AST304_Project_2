########################################################################
# Group 13: O Star Students: Nathan Shields, Brenna Chetan, Maya Joyce
# AST 304, Fall 2020
# Michigan State University
########################################################################

"""
<This module contains the neccessary functions to numerically find the Pressure at different radii within a star. First, guess the
central pressure using pressure_guess and use that value in the integrate function with a small delta_m to find the pressure throughout
the star. Adjust delta_m, eta, and xi until the M and R solutions converge. Then adjust Pc to get a range of pressures that find masses
between .1 and 1 solar masses. >
"""

import numpy as np
from eos import *
from ode import *
from astro_const import *

def stellar_derivatives(m,z,mue):
    """
    RHS of Lagrangian differential equations for radius and pressure

    Arguments
        m
            current value of the mass
        z (array)
            current values of (radius, pressure)
        mue
            ratio, nucleons to electrons.  For a carbon-oxygen white dwarf,
            mue = 2.

    Returns
        dzdm (array)
            Lagrangian derivatives dr/dm, dP/dm
    """
    # evaluate dzdm
    r= z[0] #radius
    P= z[1] #pressure
    rho = density(P,mue) #density
    drdm = 1/(ac.fourpi*(r**2)*rho)
    dPdm = (-ac.G*m)/(ac.fourpi*(r**4))

    dzdm = np.zeros_like(z)

    dzdm[0] = drdm
    dzdm[1] = dPdm
    return dzdm

def central_values(Pc,delta_m,mue):
    """
    Constructs the boundary conditions at the edge of a small, constant density
    core of mass delta_m with central pressure P_c

    Arguments
        Pc
            central pressure (Pa)
        delta_m
            core mass (kg)
        mue
            nucleon/electron ratio

    Returns
        z = array([ r, p ])
            central values of radius and pressure (m, Pa)
    """
    z = np.zeros(2)
    # compute initial values of z = [ r, p ]
    P_c = Pc   #central pressure
    rho_c = density(Pc, mue)   #central density
    r_c = (3*delta_m/(ac.fourpi*rho_c))**(1/3)   #central radius

    #create return array
    z[0] = r_c
    z[1] = P_c

    return z

def lengthscales(m,z,mue):
    """
    Computes the radial length scale H_r and the pressure length H_P

    Arguments
        m
            current mass coordinate (kg)
        z (array)
           [ r, p ] (m,Pa)
        mue
            mean electron weight

    Returns
        z/|dzdm| (kg)
    """
    r = z[0] #radius
    P = z[1] #pressure
    rho = density(P,mue)
    H_r = ac.fourpi*(r**3)*rho  #radial length scale
    H_P = ac.fourpi*(r**4)*P/(ac.G*m)  #pressure length scale

    H_z = np.zeros(2)

    #create the return array for the mass overwhich r and P change
    H_z[0] = H_r
    H_z[1] = H_P

    return H_z

def integrate(Pc,delta_m,eta,xi,mue,max_steps=10000):
    """
    Integrates the scaled stellar structure equations

    Arguments
        Pc
            central pressure (Pa)
        delta_m
            initial offset from center (kg)
        eta
            The integration stops when P < eta * Pc
        xi
            The stepsize is set to be xi*min(p/|dp/dm|, r/|dr/dm|)
        mue
            mean electron mass
        max_steps
            solver will quit and throw error if this more than max_steps are
            required (default is 10000)

    Returns
        m_step, r_step, p_step
            arrays containing mass coordinates, radii and pressures during
            integration (kg, m, Pa)
    """
    #initialize arrays
    m_step = np.zeros(max_steps)
    r_step = np.zeros(max_steps)
    p_step = np.zeros(max_steps)

    # set starting conditions using central values
    m = delta_m #initial mass
    z = central_values(Pc,m,mue) #initial r and P

    Nsteps = 0
    for step in range(max_steps):
        radius = z[0]
        pressure = z[1]
        # are we at the surface?
        if (pressure < eta*Pc):
            break

        # store the step
        r_step[step] = radius
        p_step[step] = pressure
        m_step[step] = m

        # set the stepsize
        H_z = lengthscales(m,z,mue)
        H_r = H_z[0]
        H_P = H_z[1]
        h = xi * min(H_r,H_P)

        # take a step
        z = rk4(stellar_derivatives,m,z,h,args=mue)
        m += h

        # increment the counter
        Nsteps += 1
    # if the loop runs to max_steps, then signal an error
    else:
        raise Exception('too many iterations')

    return m_step[0:Nsteps],r_step[0:Nsteps],p_step[0:Nsteps]

def pressure_guess(m,mue):
    """
    Computes a guess for the central pressure based on virial theorem and mass-
    radius relation.

    Arguments
        m
            mass of white dwarf (kg)
        mue
            mean electron mass

    Returns
        P
            guess for pressure (Pa)
    """
    # guess the pressure
    Pguess = (ac.G**5)/(ac.ke**4)*(m*(mue**2))**(10/3)

    return Pguess
