#
# Lyman break galaxy survey
# 
#
#

from Cosmology import *
from galaxies import Galaxy
import math
import numpy
import scipy.integrate

class LBG(Galaxy):
    """ Lyman break galaxy survey"""
    def __init__(self,area=4.0*math.pi,lowz=0.0,highz=1.0):
        Galaxy.__init__(self)
        self.initSurvey(area,lowz,highz)

    def luminosityFunctionParametersM(self,z):
        """ From Bouwens 2008 ApJ 686, 230 fit
        Returns (Mstar,phistar,alpha)
        Units: phistar Mpc^{-3}
        """
        z0=3.8
        
        M0=-21.02;
        sig_M0=0.09;
        M1=0.36;
        sig_M1=0.08;
   
        phi0=1.16;
        sig_phi0=0.20;
        phi1=0.024;
        sig_phi1=0.065;
   
        alpha0=-1.74;
        sig_alpha0=0.05;
        alpha1=0.04;
        sig_alpha1=0.05;

        Mstar=M0+M1*(z-z0);
        phistar=phi0*pow(10.0,phi1*(z-z0))*1.0e-3;
        alpha=alpha0+alpha1*(z-z0);

        #also get Lstar
        Lstar=Mstar

        return (Mstar,phistar,alpha,Lstar)
        
    def luminosityFunctionM(self,M,z,Mstar=-1,phistar=-1,alpha=-1):
        """ Luminosity function """

        if Mstar<0:
            Mstar,phistar,alpha,Lstar=self.luminosityFunctionParametersM(z)

        phi=self.schecterFunctionM(M,phistar,Mstar,alpha)

        return phi

    def starFormationRateFromLUV(self,LUV):
        """ dust-uncorrected SFR density
        Taken from Bouwens 2008 - assumes Salpeter
        Units SFR : Msol/yr
              LUV : erg/s/Hz at 1500A
        """
        norm=8.0e-27

        sfr=LUV/norm
        return sfr

    def surfaceDensity(self,m,z):
        """ surface density of sources N(m) """

        return Nm

    def luminosityDensity(self,Lmin,z):
        """ integrated luminosity density """
        pass

    def integratedStarFormationRate(self,Lmin,z):
        """ integrated star formation rate"""
        pass
