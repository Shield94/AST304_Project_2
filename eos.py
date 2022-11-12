########################################################################
# Group 13: O Star Students: Nathan Shields, Brenna Chetan, Maya Joyce
# AST 304, Fall 2020
# Michigan State University
########################################################################

"""
<This module defines functions that compute the equation of state of a white dwarf, giving us Pressure from mue (baryon/electron ratio) and
rho (mass density), and rho from Pressure and mue.>
"""

import astro_const as ac

def pressure(rho, mue):
    """
    Arguments
        rho
            mass density (kg/m**3)
        mue
            baryon/electron ratio

    Returns
        electron degeneracy pressure (Pascal)
    """

    # replace following lines with body of routine
    #solve for pressure
    p = ac.ke*(rho/mue)**(5/3) #ac.ke defined in astro_const is the constant term

    return p

def density(p, mue):
    """
    Arguments
        p
            electron degeneracy pressure (Pascal)
        mue
            baryon/electron ratio

    Returns
        mass density (kg/m**3)
    """

    # replace following lines with body of routine
    #solve for rho
    rho = mue*(p/ac.ke)**(3/5)
    return rho
