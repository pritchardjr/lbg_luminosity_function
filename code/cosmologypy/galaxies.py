#
# Look at galaxy properties
#
# Structured as a class for a galaxy survey
#

from Cosmology import *
import math
import numpy
import scipy.integrate

class Galaxy(Cosmology):
    """ Class for a galaxy survey """
    def __init__(self,area=4.0*math.pi,lowz=0.0,highz=1.0):
        Cosmology.__init__(self)
        self.initSurvey(area,lowz,highz)

    def initSurvey(self,area=4.0*math.pi,lowz=0.0,highz=1.0):
        """ Set the survey parameters"""
        self.area=area  #sky area in sterradians
        self.lowz=lowz
        self.highz=highz
        self.z=(highz+lowz)/2.0  # middle of survey
        volume=self.comovingVolumeToZ(highz)-self.comovingVolumeToZ(lowz)
        self.volume=volume*area/(4.0*math.pi)
        self.depth=self.comovingDistance(highz)-self.comovingDistance(lowz)
        self.width=self.comovingDistance(self.z)*math.sqrt(self.area)

    def infoSurvey(self):
        """ Basic survey information """
        surveydict={}
        surveydict['area_ster']=self.area
        surveydict['area_sqdeg']=self.area*math.pow(180.0/math.pi,2.0)
        surveydict['z']=self.z
        surveydict['lowz']=self.lowz
        surveydict['highz']=self.highz
        surveydict['volume']=self.volume  #comoving volume
        surveydict['depth']=self.depth    #comoving depth
        surveydict['width']=self.width    #comoving width - approximate

        for key in surveydict.keys():
            print key, "%g" % surveydict[key]

        return surveydict

    
    def schecterFunctionL(self,L,phistar,Lstar,alpha):
        """ Schecter function for galaxy LF
        Default to redshift ??? LBG properties
        """
        phi=phistar*numpy.pow(L/Lstar,alpha)*numpy.exp(-L/Lstar)

        return phi

    def schecterFunctionM(self,M,phistar,Mstar,alpha):
        """ Schecter function for galaxy LF
        Default to redshift ??? LBG properties
        """
        phi=phistar*math.log(10.0)/2.5;
        phi*=numpy.exp(-0.4*(M-Mstar)*(alpha+1.0)*math.log(10.0));
        phi*=numpy.exp(-numpy.exp(-0.4*(M-Mstar)*math.log(10.0)));

        return phi

    def luminosityFunctionParametersM(self,z):
        """ From Bouwens 2008 ApJ 686, 230 fit
        Returns (Mstar,phistar,alpha)
        Units: phistar Mpc^{-3}
        """
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

        return (Mstar,phistar,alpha,Lstar)
