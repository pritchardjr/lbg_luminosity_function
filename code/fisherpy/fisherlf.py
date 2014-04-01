#
# Code for fisher matrix for galaxy luminosity function
#
#
#
#

from math import *
import numpy as np
import cosmologypy.Cosmology as CCosm
from fisher import Fisher
import scipy


class FisherLF(Fisher):
    """ Fisher class for Galaxy Luminosity Function derived calculations """
    def __init__(self,paramdict={'oml':0.7,'hubble':0.7,'omm':0.3,'omk':0.0,'Ombhh':0.0225,'nscal':0.95,'Ascal':25.0}):
        Fisher.__init__(self,paramdict)
        self._fishertype="FisherCL"

    def calcFisherMatrix(self):
        #self._fisher=np.identity(self._nparam)
        raise Exception("calcFisherMatrix not defined")

################################################################

class GalaxySurvey:
    """ Class for galaxy survey and its properties
    angles in degrees

    NOTE: not entirely sure survey should own its own cosmology. May be better to separate survey fixed quantities from translation to sky
    """
    def __init__(self,name="GalaxySurvey",zmin=1.0,zmax=2.0,anglex=1.0,angley=1.0,nfield=1,nbin=1,maglim=-10.0):
        #cosm=CCosm.Cosmology()
        self.name=name
        self.zmin=zmin
        self.zmax=zmax
        self.anglex=anglex
        self.angley=angley
        self.maglim=maglim
        self.nfield=nfield
        self.nbin=nbin

    def __repr__(self):
        return "%s z=[%g,%g]; theta=[%g,%g]; nfield=%g; nbin=%g; maglim=%g" % (self.name,self.zmin,self.zmax,self.anglex,self.angley,self.nfield,self.nbin,self.maglim)

    def area(self):
        """ survey area in sq arcmin """
        return (self.anglex*60.0)*(self.angley*60.0)

    def areaSTER(self):
        """ survey area in steradians """
        return self.area()*(math.pi/180.0/60.0)**2

    def volume(self,cosm):
        """ survey volume in Mpc^3 """
        vol=cosm.comovingVolume(z1=self.zmin,z2=self.zmax)
        vol*=self.areaSTER()/(4.0*math.pi)*self.nfield;
        return vol

############################################################
class LuminosityFunction:
    """ Galaxy population from a Schecter luminosity function """
    def __init__(self,phi0=None,M0=None,alpha0=None,z=None,L0=None):

        if phi0 is None:
            if z is None:
                z=2.0
            phi0,M0,alpha0=self.getLFParam(z)

        #Can enter parameters either in terms of L or M
        if L0 is None:
            self.phi0=phi0
            self.M0=M0
            self.alpha=alpha0
            self.L0=0.0
            self.z=z
        else:
            print "NOT YET DEFINED HERE"

    def __repr__(self):
        return "%s(M0=%g phi0=%g alpha=%g z=%g)" % ("LF",self.M0,self.phi0,self.alpha,self.z)

    def getLFParam(self,z):
        """ Bouwyens+ 2008 ApJ 686, 230 fit to LF

        Units: phi Mpc^-3
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

        M=M0+M1*(z-z0);
        phi=phi0*pow(10.0,phi1*(z-z0))*1.0e-3;
        alpha=alpha0+alpha1*(z-z0);

        return (phi,M,alpha)

    def luminosityFunctionL(self,L):
        """
        Calculate galaxy luminosity function using the Schecter form 
        - Schechter (1976)
        three free parameters for the fit:
        phi0 : Mpc^{-3}
        alpha 
        L0   : ergs s^{-1}  or erg s^{-1} Hz^{-1} 
        units: phi: no. galaxies per unit comoving volume per unit luminosity
        """
        phi=self.phi0/self.L0;
        phi*=np.power(L/self.L0,self.alpha);
        phi*=np.exp(-L/self.L0);
        return phi

    def numberdensityL(self,Lmin):
        """ Number density of all galaxies with L>Lmin from integrating the Schecter function"""
        ngal=scipy.special.gammainc(1.0+self.alpha,Lmin/self.L0)
        ngal*=self.phi0
        return ngal

    def luminosityDensityL(self,Lmin):
        """ Total luminosity density above some minimum luminosity """
        lumdens=scipy.special.gammainc(2.0+self.alpha,Lmin/self.L0)
        lumdens*=self.phi0*self.L0
        return lumdens

    def luminosityFunctionM(self,M):
        """ Luminosity function dn/dM in absolute magnitude M """
        phi=self.phi0*np.log(10.0)/2.5
        phi*=np.exp(-0.4*(M-self.M0)*(self.alpha+1.0)*np.log(10.0))
        phi*=np.exp(-np.exp(-0.4*(M-self.M0)*np.log(10.0)))
        return phi

    def numberdensityM(self,Mmin):
        """ galaxies per Mpc^3 brighter than Mmin i.e. more negative M"""
        myinf= -np.inf
        myinf= -50.0
        ngal,error=scipy.integrate.quad(lambda x: self.luminosityFunctionM(x),myinf,Mmin)
        return ngal

    def absMagFromL(self,L,W=48.60):
        """
        Calculate absolute magnitude from specific luminosity Lnu
        Default value for the zero point is W=48.60 corresponding to AB
        magnitude system.

        Underlying form is redshift independent and applies to both
        L and Lnu provided appropriate zero point is used
        Units: Lnu erg s^-1 Hz^-1
        """
        M=-2.5*log10(Lnu/4.0/pi)+5-WAB
        return M

    def LfromAbsMag(self,M,W=48.60):
        """
        Calculate specific luminosity Lnu
        Default value for the zero point is W=48.60 corresponding to AB
        magnitude system
        Units: Lnu erg s^-1 Hz^-1
        """
        x=-0.4*(M-5+WAB)
        Lnu=4.0*pi*pow(10.0,x)
        return Lnu

    def magFromAbsMag(self,M,dL,z=0):
        """
        dL is luminosity distance in parsec
        redshift dependence only appropriate for Mnu from Lnu not M from L
        """
        m=M-5.0*+5.0*log10(dL)-2.5*log10(1+z)
        return m

    def absMagFromMag(self,m,dL,z=0):
        """
        dL is luminosity distance in parsec
        redshift dependence only appropriate for Mnu from Lnu not M from L
        """
        M=m+5.0*-5.0*log10(dL)+2.5*log10(1+z)
        return M   
        
