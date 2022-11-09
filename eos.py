########################################################################
# Team 13: Nathan Shields, Brenna Chetan, Maya Joyce
# AST304, Fall 2020
# Michigan State University
########################################################################

"""
<Description of this module goes here: what it does, how it's used.>
"""

import astro_const as ac
from numpy import pi

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
    p = .2*(3/(8*pi))**(2/3)*ac.h**2/ac.m_e*(rho/(mue*ac.m_u))**(5/3)
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
    rho = mue*ac.m_u*(5*p*ac.m_e/ac.h**2*(3/(8*pi))**-(2/3))**(3/5)
    return rho
