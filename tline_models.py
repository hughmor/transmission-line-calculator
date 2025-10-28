"""
Impedance calculator based on the equations from Wadell's 'Transmission line design handbook'

I started by writing a function for a transmission line topology that wasn't in KiCAD's calculator tool or any online calculators I could find.
Once I wrote this, I realized it would be cool to rewrite my own version of this kind of calculator tool.
Some features that would be nice to have:
- Calculation and plotting of impedances for a range of parameters
- Optimization over parameters
- Visualization of waveguides (this is one thing that current tools actually do well), but I could start with a little command line ASCII art version of each one
- Error/tolerance analysis when plotting
"""

import numpy as np
from scipy.special import ellipk as K
from scipy.constants import mu_0, epsilon_0 as e_0

def coplanar_strips(a, b, h, e_r):
    """
    Characteristic impedance of the 'Coplanar Strips' waveguide from Sec. 3.4.6 of Wadell's 'Transmission line design handbook'

    TODO: I'm recomputing every k value a bunch of times; this can be optimized to reduce the number of elliptic integrals being calculated, but it doesn't seem to have an execution time problem right now. Maybe once I add optimization over parameters I'll need to do that.
    TODO: Right now this only computes the characteristic impedance. I should also incorporate the loss equations. 
    """

    # 
    k1 = lambda a, b, h: np.sinh(np.pi*a/4/h) / np.sinh(np.pi*b/4/h)
    k1p = lambda a, b, h: np.sqrt(1-k1(a,b,h)**2)
    k = lambda a,b,h: a/b
    kp = lambda a, b, h: np.sqrt(1-k(a,b,h)**2)
    e_eff = lambda a, b, h, e_r: 1 + ((e_r-1)*K(kp(a,b,h))*K(k1(a,b,h)))/(2*K(k(a,b,h))*K(k1p(a,b,h)))
    eta_0 = np.sqrt(mu_0/e_0)
    Z0 = lambda a,b,h,e_r: eta_0 * K(k(a,b,h)) / np.sqrt(e_eff(a,b,h,e_r)) / K(kp(a,b,h))
    return Z0(a, b, h, e_r)


if __name__ == '__main__':
    a = 0.05e-3
    b = 0.13e-3
    h = 1e-3
    e_r = 1.47
    Z0 = coplanar_strips(a, b, h, e_r)
    print(f'The characteristic impedance of a coplanar strip waveguide with total width {b:.2e}, spacing between strips {a:.2e}, height {h:.2e}, and relative substrate dielectric constant {e_r:.2f} is: {Z0:.4f} Ohms')

