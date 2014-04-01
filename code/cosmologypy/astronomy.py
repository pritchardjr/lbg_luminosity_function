#
# Astronomy type stuff
#
#
#

import Cosmology
import cosmoconstants
import math
import numpy

defaultcosm=Cosmology.Cosmology()

####### Core magnitude definitions ###############################

def apparentFromAbsoluteMagnitude(M,dL):
    """ convert absolute M to relative magnitude m
    dL in Mpc
    """
    m=M+5.0*numpy.log10(dL/1.0e6)-1.0
    return m

def absoluteFromApparentMagnitude(m,dL):
    """ convert absolute M to relative magnitude m
    dL in Mpc
    """
    M=m-5.0*numpy.log10(dL/1.0e6)+1.0
    return m

def distanceModulus(m,M):
    """ Calculate distance modulus """
    mu=m-M
    return mu

def apparentMagnitudeFromLuminosityNu(Lnu,dL,z=0,Z=-48.60,band='AB'):
    """ Basic definition of apparent magnitude
    Lnu in erg/s/Hz
    dL in Mpc
    """
    Flux=Lnu*(1+z)/(4.0*math.pi*dL*dL)
    Flux/=cosmoconstants.MpcCGS**2
    m= -2.5*numpy.log10(Flux)-Z
    return m

def luminosityNuFromAbsoluteMagnitude():
    pass

####################################################################

def magnitudeFromLnu(Lnu,dL=-1,z=0,cosm=defaultcosm):
    """relative magnitude from Luminosity
    Units: dL in Mpc
           L in erg/s/Hz
    """
    if dL<0: dL=cosm.luminosityDistance(z);

    #calculate flux assuming luminosity in a freq. interval
    dLCGS=dL*cosmoconstants.MpcCGS
    flux=Lnu*(1.0+z)/(4.0*math.pi*dLCGS*dLCGS);
    m= -2.5*numpy.log10(flux)-48.60;

    return m

def magnitudeFromL(L,dL= -1,z=0,cosm=defaultcosm):
    """relative magnitude from Luminosity"""
    if dL<0: dL=cosm.luminosityDistance(z);

    #calculate flux assuming luminosity in a freq. interval
    dLCGS=dL*cosmoconstants.MpcCGS
    flux=L/(4.0*math.pi*dLCGS*dLCGS);
    m= -2.5*numpy.log10(flux)-48.60;

    return m

def absMagnitudeFromLnu(Lnu,dL= -1,z=0,cosm=defaultcosm):
    """absolute magnitude from Luminosity"""
    m=magnitudeFromLnu(Lnu,dL,z,cosm)
    M=m+2.5*numpy.log10(1+z)-5.0*numpy.log10(dL/parsec)+5.0
    return M

def luminosityFromMagnitude(m,dL= -1,z=0,cosm=defaultcosm):
    """ luminosity from relative magnitude """
    Lnu=pow(10.0,-0.4*(m+48.60))*4*math.pi/(1.0+z)
    return Lnu

def luminosityFromAbsoluteMagnitude(M,dL= -1,z=0,cosm=defaultcosm):
    """Luminosity from absolute magnitude"""
    m=M-2.5*numpy.log10(1+z)+5.0*numpy.log10(dL/parsec)-5.0
    Lnu=luminosityFromMagnitude(m,dL,z,cosm)
    return Lnu
