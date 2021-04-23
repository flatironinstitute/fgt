1. file description:

fgt3d.f     - 3d FGT driver using plane waves

g3ddirect.f - direct calculation subroutines

g3drouts.f  - contains form and eval Hermite and local expansion subroutines,
                     Hermite to local, Hermite to PW, PW to local,

g3dterms.f  - calculates number of terms for various expansions on each level


2. naming conventions for subroutines/files:

h  - Hermite 
l  - local, i.e., Taylor
pw - PW

c  - charge
d  - dipole

p  - potential
g  - potential + gradient
h  - potential + gradient + hessian



3. order conventions:

hessian order      - xx, yy, zz, xy, xz, yz
h and l expansions - inner loop - x; middle loop - y; outer loop - z


4. scaling conventions:

The Hermite function h_n(x) is devided by sqrt(n!) to avoid unnecessary overflow.
