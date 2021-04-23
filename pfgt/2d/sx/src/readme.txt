1. file description:

fgt2d.f     - 2d FGT driver

g2ddirect.f - direct calculation subroutines

g2dhlall.f  - contains form and eval Hermite and local expansion subroutines,
                     Hermite to local, Hermite to SOE and SOE/X, SOE to local,
		     SOE/X to local subroutines, and 1D translation matrices (setup stage)

g2dsoeall.f - contains form and eval SOE expansion subroutines, and SOE shift
		     subroutines

g2dsxall.f  - contains form and eval SOE/X expansion subroutines, and SX shift
		     subroutines

g2dterms.f  - calculates number of terms for various expansions on each level


----------------------------------------------------------------------------------
code to be removed later:

fgt2d2.f, g2dsoeall2.f, g2dsxall2.f - use vertices as SOE expansion center and
               edge midpoints as SX expansion center
--------------------------------------------------------------------------------




	       

2. naming conventions for subroutines/files:

h - Hermite 
l - local, i.e., Taylor
s - SOE
x - plane wave

c - charge
d - dipole

p - potential
g - potential + gradient
h - potential + gradient + hessian



3. order conventions:

hessian order      - xx, xy, yy
h and l expansions - inner loop - x; outer loop - y
SOE expansions     - pp, pm, mp, mm
SX expansions      - px, mx, xp, xm



4. scaling conventions:

The Hermite function h_n(x) is devided by sqrt(n!) to avoid unnecessary overflow.
