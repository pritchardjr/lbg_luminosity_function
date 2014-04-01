#
# List of commonly used cosmological constants
#
import math

TINY=1.0e-20
H0=3.24076e-18 # 1/s : Hubble Constant in 1/s without h i.e. H0=100 km/s/Mpc
c=2.99792e8 #m/s : Speed of light in MKS
c_CGS=2.99792e10 #m/s : Speed of light in CGS
cH0_Mpc=2997.92  # Mpc : c/H0 in Mpc
CRITDENMSOLMPC=2.7755e11 #Current critical density, in Msol/Mpc^3, not including h^2!

h=0.7
Yp=0.247
fHe=Yp/4/(1-Yp)
omegar=4.15e-5/h/h
omegacmb=2.47e-5/h/h
Tcmb=2.726 #Kelvin Temperature of the CMB today
rhocrit=1.879e-26 # kg/m^-3 Critical density without h^2 factor
mproton = 1.67262e-27 # kg Proton mass
sigmaT_CGS = 6.65246e-25 # Thompson cross-section in cgs
sigmaT = 6.65246e-29 # Thompson cross-section in cgs

##### Units
parsec=3.0857e16 # Meter  Parsec in Meters
Mpc=1.0e6*parsec
MpcCGS=Mpc*100.0
Mpc2year=3.26382e6 # Year: 1Mpc in Years with c=1
hour=3600.0 # Second: 1 hour in seconds
deg2steradians=(math.pi/180)**2

Lsol=3.839e26  # Watts - Solar luminosity
LsolCGS=3.839e33  # erg/s

##### 21 cm physics
lambda21=0.2112 # Meter : 21cm wavelength in meters
nu21=1.420405e9 # Hz
