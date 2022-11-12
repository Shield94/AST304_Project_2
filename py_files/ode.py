########################################################################
# Group 13: O Star Students: Nathan Shields, Brenna Chetan, Maya Joyce
# AST 304, Fall 2020
# Michigan State University
########################################################################

"""
This module contains 4 methods of numerically solving ODEs (Forward Euler, 2nd order Runga-Kutta and 4th order Runga-Kutta).
Each function advances a numerical integration one step.
"""

# all routines that take a single step should have the same interface
# fEuler is complete, except for documentation. you can use this as a pattern
# for the other two routines.
def fEuler(f,t,z,h,args=()):
    """
    Advances a numerical integration one step with the Forward Euler method: fEuler(f,t,z,h,args = ()

    Arguments
        f(t,z,...)
            function that contains the RHS of the equation dz/dt = f(t,z,...)

        t
            inital time of numerical integration

        z
            position and velocity array: array(x-position,y-position,x-velocity,y-velocity)

        h
            step size for integration

        args (tuple, optional)
            additional arguments to pass to f

    Returns
        znew = z(t+h)
    """

    # The following trick allows us to pass additional parameters to f
    # first we make sure that args is of type tuple; if not, we make it into
    # that form
    if not isinstance(args,tuple):
        args = (args,)

    # when we call f, we use *args to pass it as a list of parameters.
    # for example, if elsewhere we define f like
    # def f(t,z,x,y):
    #    ...
    # then we would call this routine as
    # znew = fEuler(f,t,z,h,args=(x,y))
    #

    return z + h*f(t,z,*args)


# Runge-Kutta step and a fourth order Runge-Kutta step.

def rk2(f,t,z,h,args=()):
    """
    Advances a numerical integration one step with the 2nd order Runga-Kutta method: rk2(f,t,z,h,args = ()

    Arguments
        f(t,z,...)
            function that contains the RHS of the equation dz/dt = f(t,z,...)

        t
            inital time of numerical integration

        z
            position and velocity array: array(x-position,y-position,x-velocity,y-velocity)

        h
            step size for integration

        args (tuple, optional)
            additional arguments to pass to f

    Returns
        znew = z(t+h)
    """

    if not isinstance(args,tuple):
        args = (args,)

    Zmp = z + (h/2)*f(t,z,*args)

    return z + h*f(t + h/2,Zmp,*args)


def rk4(f,t,z,h,args=()):
    """
    Advances a numerical integration one step with the 4th order Runga-Kutta method: rk4(f,t,z,h,args = ())

    Arguments
        f(t,z,...)
            function that contains the RHS of the equation dz/dt = f(t,z,...)

        t
            inital time of numerical integration

        z
            position and velocity array: array(x-position,y-position,x-velocity,y-velocity)

        h
            step size for integration

        args (tuple, optional)
            additional arguments to pass to f

    Returns
        znew = z(t+h)
    """

    if not isinstance(args,tuple):
        args = (args,)

    k1 = f(t,z,*args)

    k2 = f(t + (h/2),z + k1*(h/2),*args)

    k3 = f(t + (h/2),z + k2*(h/2),*args)

    k4 = f(t + h,z + k3*h,*args)

    return z + (h/6)*(k1+2*k2+2*k3+k4)
