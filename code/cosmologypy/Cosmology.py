#
# This defines a Cosmology class
#
# All functions should be numpy compatible, so allowing vector arguments
#
#

import cosmoconstants
from math import *
import numpy
import scipy.interpolate
import scipy.integrate

MAXBINS=7  #constant for Press-Schecter calculations
from cosmoconstants import CRITDENMSOLMPC

class Cosmology:
    """ This is the cosmology class
    All functions should be numpy compatible, so allowing vector arguments
    """
    
    #def __init__(self,Omegam=0.3,Omegal=0.7,Omegak=0.0,Omegabhh=0.0225,hubble=0.7,nscal=0.95,Ascal=25.0,sigma8=0.9):
    def __init__(self,Omegam=0.3,Omegal=0.7,Omegak=0.0,Omegabhh=0.046*0.7*0.7,hubble=0.7,nscal=1.0,Ascal=25.0,sigma8=0.9):        
        """ Initialise the cosmological parameters: defaults specified"""
        #Establish core cosmological parameter set
        self.Omegam=Omegam
        self.Omegal=Omegal
        self.Omegak=Omegak
        self.hubble=hubble
        self.Omegabhh=Omegabhh
        self.nscal=nscal
        self.Ascal=Ascal*1.0e-10
        self.sigma8=sigma8

        self.w= -1

        #Basic derived quantities
        h=self.hubble
        self.h=h  #LAZY
        self.Omegab=self.Omegabhh/h/h
        self.Omegar=4.15e-5/h/h;
        self.Omegacmb=2.47e-5/h/h
        #self.Omegam=self.Omegam-self.Omegar  #rescale matter to allow for small radiation 
        self.Omegac=self.Omegam-self.Omegab
        self.Omegachh=self.Omegac*h*h
        self.Omegamhh=self.Omegam*h*h
        self.Omeganu=0.0
        # more specialised quantities

        self.Yp=0.247
        self.fHe=self.Yp/4/(1-self.Yp)
        self.TCMB=2.726   #in Kelvin

        #set up parameters for perturbation theory
        self.resetPowerSpectrum()

        #sigma(m) spline
        self.sigma_s=None
        self.sigma_sd=None
        
        return

    def __repr__(self):
        return "%s(Hubble=%g,Omm=%g,Oml=%g)" % ("Cosmology",self.hubble,self.Omegam,self.Omegal)

    def info(self):
        """ print basic information about the cosmology parameters"""
        cosmdict={}
        cosmdict['Omegam']=self.Omegam
        cosmdict['Omegal']=self.Omegal
        cosmdict['Omegak']=self.Omegak
        cosmdict['Omegabhh']=self.Omegabhh
        cosmdict['hubble']=self.hubble
        cosmdict['Omegar']=self.Omegar
        cosmdict['nscal']=self.nscal
        cosmdict['Ascal']=self.Ascal
        cosmdict['Omegachh']=self.Omegachh
        
        for key in cosmdict.keys():
            print key, cosmdict[key]
            
        return cosmdict

######## Summary - c.f. CosmoCalc output, which this mimics
    def summary(self,z):
        """ List common distances and quantities at redshift z """
        cosmdict={}
        cosmdict['z']=z
        cosmdict['a']=1.0/(1.0+z)
        cosmdict['dL']=self.luminosityDistance(z)
        cosmdict['dA']=self.angularDistance(z)
        cosmdict['dC']=self.comovingDistance(z)
        cosmdict['dP']=self.properDistance(z)
        cosmdict['E(z)']=self.EZ(z)
        cosmdict['mu']=self.distanceModulus(z)
        cosmdict['tlook']=self.lookbackTime(z)
        cosmdict['age']=self.ageUniverse(z)
        cosmdict['VolC']=self.comovingVolumeToZ(z)
        cosmdict['rhoCrit']=self.criticalDensity(z)

        for key in cosmdict.keys():
            print key, "%g" % cosmdict[key]
            
        return cosmdict

####### Memberfunctions
    def getH(self):
        return self.hubble

    def getOmega0(self):
        return self.Omegam

    def getOmega0hh(self):
        return self.Omegamhh
    
######## Basic distance functions
    def EZ(self,z):
        """ Calculate E(z)=H(z)/H0 """
        E2=self.Omegam*numpy.power(1.+z,3.)
        E2+=self.Omegal
        E2+=self.Omegar*numpy.power(1.+z,4.)
        E2+=self.Omegak*numpy.power(1.+z,2.)
        return numpy.sqrt(E2)

    def hubbleZ(self,z):
        """ Hubble parameter in 1/s """
        EZ=self.EZ(z)
        return EZ*cosmoconstants.H0*self.hubble

    def hubbleZMPC(self,z):
        """ Hubble parameter in 1/Mpc """
        return self.hubbleZ(z)*cosmoconstants.Mpc/cosmoconstants.c

    def comovingDistance(self,z):
        """ Calculate comoving distance in Mpc from z=0 to z=z
        This is c x conformalTime and gives the comoving distance
        """
        distance,error=scipy.integrate.quad(lambda x: 1.0/self.EZ(x),0,z)
        return distance*cosmoconstants.c/(cosmoconstants.H0*self.hubble)/cosmoconstants.Mpc

    def dtdz(self,z):
        """ d proper time/d redshift in seconds"""
        dtdz=1.0/self.hubbleZ(z)/(1.0+z)
        return dtdz

    def properDistance(self,z):
        """ Calculate proper distance in Mpc """
        distance,error=scipy.integrate.quad(lambda x: 1.0/self.EZ(x)/(1.0+x),0,z)
        return distance*cosmoconstants.c/(cosmoconstants.H0*self.hubble)/cosmoconstants.Mpc

    def angularDistance(self,z):
        """ proper angular diameter distance """
        return self.comovingDistance(z)/(1.0+z)

    def luminosityDistance(self,z):
        """ proper luminosity distance """
        return self.comovingDistance(z)*(1.0+z)

    def comovingVolumeElement(self,z):
        """ comoving volume element dVol/dz/dOmega in Mpc^3 """
        return math.pow(self.comovingDistance(z),2)/self.hubbleZMPC(z)

    def comovingVolumeToZ(self,z):
        """ total comoving volume from z=0 to z=z
        REDUNDANT
        """
        return self.comovingVolume(z2=z)

    def comovingVolume(self,z2,z1=0.0):
        """ total comoving volume from z1 to z2 in Mpc^3 """
        volume,error=scipy.integrate.quad(self.comovingVolumeElement,z1,z2)
        return volume*4.0*math.pi

    def distanceModulus(self,z):
        """ Calculate the distance modulus mu=5 log_10(dL/10pc)"""
        return 5.0*numpy.log10(self.luminosityDistance(z)*1.0e5)

    def lookbackTime(self,z):
        """ Look back time in years"""
        return self.properDistance(z)*cosmoconstants.Mpc2year

    def comovingHorizon(self,z):
        """ Calculate comoving horizon in Mpc i.e. from z=infty to z=z"""
        distance,error=scipy.integrate.quad(lambda x: 1.0/self.EZ(x),z,numpy.inf)
        return distance*cosmoconstants.cH0_Mpc/self.hubble

    def properHorizon(self,z):
        """ Calculate conformal horizon in Mpc i.e. from z=infty to z=z"""
        distance,error=scipy.integrate.quad(lambda x: 1.0/self.EZ(x)/(1.0+x),z,numpy.inf)
        return distance*cosmoconstants.cH0_Mpc/self.hubble

    def ageUniverse(self,z):
        """ age of the Universe at redshift z """
        return self.properHorizon(z)*cosmoconstants.Mpc2year

####### Redshift dependent abundances
    def omegamZ(self,z):
        """Omega_m(z) """
        return self.Omegam*(1+z)**3/self.EZ(z)**2

    def omegalZ(self,z):
        """Omega_lambda(z) """
        return self.Omegal/self.EZ(z)**2

    def getScale(self):
        return self.scale


####### Miscellaneous quantities
    def scaleEquality(self):
        """ scale factor at matter radiation equality"""
        return self.Omegar/self.Omegam

    def kEquality(self):
        """ return comoving wavenumber corresponding to horizon size at
        matter/radiation equality in 1/Mpc"""
        keq2=2.0*self.Omegam/(cosmoconstants.cH0_Mpc/self.hubble)**2/self.scaleEquality()
        return numpy.sqrt(keq2)
    
######## Abundances
    def criticalDensity(self,z):
        """ returns critical density in kg/m^3 """
        rho=cosmoconstants.rhocrit*self.hubble**2
        return rho*self.EZ(z)**2

    def nbaryon(self,z):
        """ number density of hydrogen atoms"""
        nb=self.Omegabhh*cosmoconstants.rhocrit/cosmoconstants.mproton*numpy.power(1.0+z,3.0)
        return nb

    def nhydrogen(self,z):
        """ number density of hydrogen atoms"""
        nh=(1-self.Yp)*self.nbaryon(z)
        return nh

####### Growth of structure
    def linearGrowth(self,z):
        """ Linear growth function D(z)"""
        growth=5.0*self.Omegam/2.0*self.EZ(z)
        f=lambda x: numpy.power(x*self.EZ((1.0/x)-1.0),-3.0)
        II,error=scipy.integrate.quad(f,0.0,1.0/(1.0+z))
        return growth*II
    
    def fGrowth(self,z,gamma=4.0/7.0):
        """ Approximation to growth function f=dlnD/dlna ~Omega_m^gamma """
        return numpy.power(self.omegamZ(z),gamma)


####### BBKS Power spectrum
    def transferBBKS(self,k):
        """ BBKS transfer function """
        x=k/self.kEquality()
        transfer=numpy.log(1.0+0.171*x)/(0.171*x)
        transfer*=(1 + 0.284*x+(1.18*x)**2 + (0.399*x)**3 + (0.490*x)**4)**-0.25
        return transfer
    
    def powerBBKS(self,k):
        """ BBKS power spectrum """
        deltaH=1.9e-5
        #deltaH=numpy.sqrt(self.Ascal)
        powerprim=2.0*pi**2*numpy.power(cosmoconstants.cH0_Mpc/self.hubble,3)*deltaH**2
        powerprim*=numpy.power(k*cosmoconstants.cH0_Mpc/self.hubble,self.nscal)
        power=powerprim*self.transferBBKS(k)**2
        return power

    def powerBBKSZ(self,k,z):
        """ BBKS power spectrum at redshift z
        NOTE: doesn't handle z as a vector, although that would be nice
        """
        return self.powerBBKS(k)*numpy.power(self.linearGrowth(z),2.0)

#######################################################################
    # Perturbation theory functions - originally from S. Furlanetto code
#######################################################################

    def DeltaC(self,zCurrent):
        """
        Returns fitted overdensity of collapsed halos *relative to
        critical*, according to chosen cosmology.  Taken from Bryan and
        Norman (1998); BL01 eq. 22.
        """
        d = self.omegamZ(zCurrent)-1
        if ((self.Omegal == 0.0) and (self.Omegam < 1.0)):
            return (18.0*pi*pi + 60.0*d - 32.0*d*d)
        if (abs(self.Omegal + self.Omegam-1.0) < 1.0e-3): 
            return (18.0*pi*pi + 82.0*d - 39.0*d*d)

        return 0.0
    

    def delCrit0(self,z):
        """
        Fit to the *linear* overdensity at which virialization takes place, 
        it is 1.69 with slight cosmology dependence.
        """
        t1 = 0.15*pow(12.0*pi,2.0/3.0);
        omZ = self.omegamZ(z);
        if (omZ > (1.0 - 1.0e-5)):
            return t1;
        else: 
            if (self.Omegal < 1.0e-5):
                t1 *= pow(omZ,0.0185);
                return t1;
            else:
                t1 *= pow(omZ,0.0055);
                return t1;

    def growthFac(self,z):
        """Returns inverse of linear growth from z to 0:
        * D = D1(z)/D1(z=0); from Eisenstein and Hu. """
        omZ = self.omegamZ(z)
        lamZ = self.omegalZ(z)
        D = ((1.0 + self.zEquality)/(1.0 + z)*5.0*omZ/2.0*pow(pow(omZ,4.0/7.0) - lamZ + (1.0+omZ/2.0)*(1.0+lamZ/70.0),-1.0))
        D /= ((1.0 + self.zEquality)*5.0/2.0*self.Omegam*pow(pow(self.Omegam,4.0/7.0) - self.Omegal + (1.0+self.Omegam/2.0)*(1.0+self.Omegal/70.0),-1.0))
        return D

    # **********************************************************************
    # ************************ Collapse Fraction ***************************
    # *********************************************************************/


    def fCollPSExact(self,z,mMin=None):
        """Collapse fraction above mMin for Press-Schechter mass function.
        Note mMin is an optional parameter; if not input it is taken to be the mass above halos with 10^4 K (the minimum cooling mass) """

        if mMin is None:
            mMin = self.coolMass(z);
        sigma = self.sigma0fM(mMin);
        dCritZ = self.delCrit0(z)/self.growthFac(z);
        ans = erfc(dCritZ/sqrt(2.0)/sigma);
        return ans;

#######################################################################
#####Power Spectrum Functions - originally from S. Furlanetto
#######################################################################
    def linearPowerSpectrum(self,k,z):
        """ Compute the linear power spectrum at redshift z.
        Uses the Eisenstein & Hu fit to the transfer function. 
        k = comoving wavenumber (Mpc^-1) """
        return self.powerSpectrum(k)*pow(self.growthFac(z),2.0)

    def powerSpectrum(self,k):
        """ Return power spectrum, assuming the EH transfer function
        k = wavenumber (comoving Mpc^-1) """
        temp = 2.0*numpy.power(pi*self.sNorm*self.TFMaster(k),2.0)
        temp *= numpy.power(k,self.nscal)
        return temp


    def resetPowerSpectrum(self):
        """ Calculates constants necessary for the calculation of the
        power spectrum and sigma """
  
        self.thetaCMB = self.TCMB/2.7
        #Set up transfer function parameters
        self.TFSetParameters()
        #Normalize sigma8 at present time
        scale = 8.0/self.hubble
        self.scale=scale
        acc = 1.0e-6
        sigma8norm = scipy.integrate.quad(self.sigmatop2,1.0e-6,0.001/scale,epsrel=acc)[0];
        sigma8norm += scipy.integrate.quad(self.sigmatop,log(0.001/scale),log(0.1/scale),epsrel=acc)[0];
        sigma8norm += scipy.integrate.quad(self.sigmatop,log(0.1/scale),log(1.0/scale),epsrel=acc)[0];
        sigma8norm += scipy.integrate.quad(self.sigmatop,log(1.0/scale),log(10.0/scale),epsrel=acc)[0]; 
        sigma8norm += scipy.integrate.quad(self.sigmatop,log(10.0/scale),log(100.0/scale),epsrel=acc)[0];
        sigma8norm = sqrt(sigma8norm);
        self.sNorm = self.sigma8/sigma8norm;


    def TFSetParameters(self):
        """
        Sets up transfer function parameters.  Ultimately from Hu and 
        Eisenstein (1998) online code 
        (http://www.sns.ias.edu/~whu/transfer/transfer.html)."""
        fNu = self.Omeganu/self.Omegam
        fBaryon = self.Omegab/self.Omegam

        thetaCMB=self.thetaCMB
        om0hh=self.Omegamhh
        ombhh=self.Omegabhh

        zEquality = 2.50e4*om0hh*pow(thetaCMB,-4.0) - 1.0;
        self.zEquality=zEquality
        kEquality = 0.0746*om0hh*pow(thetaCMB,-2.0);
        zDrag = 0.313*pow(om0hh,-0.419)*(1.0+0.607*pow(om0hh,0.674));
        zDrag = 1.0+zDrag*pow(ombhh,0.238*pow(om0hh,0.223));
        zDrag *= 1291.0*pow(om0hh,0.251)/(1.0 + 0.659*pow(om0hh,0.828));

        yD = (1.0 + zEquality)/(1.0+zDrag);
        RDrag = 31.5*ombhh*pow(thetaCMB,-4.0)*1000.0/(1.0+zDrag);
        REquality = 31.5*ombhh*pow(thetaCMB,-4.0)*1000.0/(1.0+zEquality);
        soundHorizon = 2.0/3.0/kEquality*sqrt(6.0/REquality)*log((sqrt(1.0+RDrag) + sqrt(RDrag+REquality))/(1.0+sqrt(REquality)));
        #Fit from Hu and Eisenstein
        soundHorizon = 44.5*log(9.83/om0hh)/sqrt(1.0+10.0*pow(ombhh,0.75));
        self.soundHorizon=soundHorizon

        fC = 1.0 - fNu - fBaryon;
        fCB = 1.0 - fNu;
        pC = -(5.0-sqrt(1.0+24.0*fC))/4.0;
        pCB = -(5.0-sqrt(1.0+24.0*fCB))/4.0;
        fNUB = fNu + fBaryon;

        alphaNu = (fC/fCB)*(2.0*(pC+pCB)+5.0)/(4.0*pCB+5.0);
        alphaNu *= (1.0-0.553*fNUB+0.126*pow(fNUB,3.0));
        alphaNu /= 1.0 - 0.193*sqrt(fNu) + 0.169*fNu;
        alphaNu *= pow(1.0+yD,pC-pCB);
        alphaNu *= (1.0+(pCB-pC)/2.0*(1.0+1.0/(4.0*pC+3.0)/(4.0*pCB+7.0))/(1.0+yD));
        self.alphaNu=alphaNu
        self.betaC = 1.0/(1.0 - 0.949*fNUB);      
        

    def TFMaster(self,k,iDeriv=False):
        """ Calculates and returns TF(k), and if iDeriv=1, calculates and stores
        dTF/dk.  Transfer function is from Hu and Eisenstein online code:
        http://www.sns.ias.edu/~whu/transfer/transfer.html
        Derivative calculation by Rennan Barkana """
        
        dTFdk=None
        t1=0.0
        dt1dqe=0.0

        q = k*self.thetaCMB*self.thetaCMB/self.Omegamhh
        gammaEff = (sqrt(self.alphaNu) + (1.0 - sqrt(self.alphaNu))/(1.0 + numpy.power(0.43*k*self.soundHorizon,4.0)))
        qEff = q/gammaEff

        tfMaster = numpy.log(exp(1.0)+1.84*self.betaC*sqrt(self.alphaNu)*qEff);
        if (iDeriv == 1):
            t1 = tfMaster
            dt1dqe = 1.84*self.betaC*sqrt(self.alphaNu)/exp(tfMaster)
            
        tfMaster = tfMaster/(tfMaster + qEff*qEff*(14.4+325.0/(1.0+60.5*numpy.power(qEff,1.11))))
        
        if iDeriv:
            t2 = tfMaster;
            v = t1/t2;
            dv = dt1dqe + 2.0*qEff*(14.4+325.0/(1.0+60.5*numpy.power(qEff,1.11))) - 325.0*60.5*1.11*numpy.power(qEff,2.11)/numpy.power(1.0+60.5*numpy.power(qEff,1.11),2.0);
            dgamedk = -(1.0-sqrt(self.alphaNu))*4.0*numpy.power(k,3.0)*pow(0.43*self.soundHorizon,4.0)/numpy.power(1.0+numpy.power(0.43*self.soundHorizon*k,4.0),2.0);
            dTFdk = self.thetaCMB*self.thetaCMB/self.Omegamhh*((v*dt1dqe-t1*dv)/v/v)*(1.0/gammaEff-k*dgamedk/gammaEff/gammaEff)

        if dTFdk is not None:
            return tfMaster,dTFdk
        else:
            return tfMaster

############################################################
    # Bias Calculations 
###########################################################
    def biasPS(self,mass,z):
        """
        Linear bias as a function of halo mass for Press-Shechter halos. 
        * See, e.g., Cooray & Sheth (2002), s. 3.4. """
        deltasc = self.delCrit0(z)/self.growthFac(z)
        deltasc0 = self.delCrit0(0.0)
        sig = self.sigma0fM(mass)
        bias = 1.0 + (deltasc*deltasc/sig/sig - 1.0)/deltasc0
        return bias

    def biasmPS(self,mass,z):
        """
        Linear bias as a function of halo mass for Press-Shechter halos. 
        * See, e.g., Cooray & Sheth (2002), s. 3.4.
        Uses splined values of sigma(m)
        """
        deltasc = self.delCrit0(z)/self.growthFac(z)
        deltasc0 = self.delCrit0(0.0)
        sig = self.sigm(mass)
        bias = 1.0 + (deltasc*deltasc/sig/sig - 1.0)/deltasc0
        return bias

    def biasmST(self,mass,z):
        """ Sheth-Torman Halo bias (linear)
        Uses splined values of sigma(m)
        """
        q=0.75
        p=0.3

        deltasc = self.delCrit0(z)/self.growthFac(z);
        deltasc0 = self.delCrit0(0.0);
        sig = self.sigm(mass)
        nu=deltasc*deltasc/sig/sig;
        bias = 1.0 + (q*nu - 1.0)/deltasc0;
        bias += 2.0*p/deltasc0/(1.0+numpy.power(q*nu,p));
        
        return bias

    def biasmTinker(self,mass,z):
        """ Tinker Halo bias (linear)
        Taken from Robertson (2010),
        but originates in Tinker+ (2010) as fit to numerical sim
        Uses splined values of sigma(m)
        """
        A=1.0
        a=0.1325
        B=0.183
        b=1.5
        C=0.265
        c=2.4

        deltasc = self.delCrit0(z)/self.growthFac(z);
        deltasc0 = self.delCrit0(0.0);
        sig = self.sigm(mass)
        nu=deltasc/sig;
        bias = 1.0
        bias -= A*pow(nu,a)/(pow(nu,a)+pow(deltasc0,a));
        bias += B*pow(nu,b)+C*pow(nu,c)
        
        return bias


############################################################
    # Press-Schechter functions ***********************
############################################################
    # Routines to calculate Press-Schecter mass function.  
    # Transfer function from Eisenstein and Hu (1998); does not include 
    # possibility of hot dark matter or neutrinos (despite the parameters!).
    # Translated from Fortran code provided by Rennan Barkana, transfer
    # function comes ultimately from Wayne Hu's website.


    def dndlM(self,z,tM):
        """ Calculates Press-Schechter mass function 
        tM = halo mass (Msun) """
        dCritZ = self.delCrit0(z)/self.growthFac(z);
        sigM,dsdM = self.sigma0fM(tM,True);
        dlsdlM = tM*dsdM/sigM;
        tdn = (sqrt(2.0/pi)*dCritZ*fabs(dlsdlM)*exp(-dCritZ*dCritZ/(2.0*sigM*sigM))/(tM*sigM));
        tdn *= CRITDENMSOLMPC*self.Omegamhh;
        return tdn;

    def dndlMSheth(self,z,tM):
        """    Calculates Sheth-Tormen (1999) mass function
        tM = halo mass (Msun)"""
        a=0.707
        A=0.322
        p=0.3

        dCritZ = self.delCrit0(z)/self.growthFac(z);
        sigM,dsdM = self.sigma0fM(tM,1);
        dlsdlM = tM*dsdM/sigM;
        nuc = dCritZ/sigM;
        tdn = A*sqrt(2.0*a/pi)*fabs(dlsdlM)/tM*nuc*(1.0+pow(nuc*nuc*a,-p));
        tdn *= exp(-a*nuc*nuc/2.0);
        tdn *= CRITDENMSOLMPC*self.Omegamhh;
        return tdn;


    def dndlMJenkins(self,z,tM):
        """ Calculates Jenkins et al. (2001) 
        *   tM = halo mass (Msun)
        """
        a=0.73
        A=0.353
        p=0.175

        dCritZ = self.delCrit0(z)/self.growthFac(z);
        sigM,dsdM = self.sigma0fM(tM,1);
        dlsdlM = tM*dsdM/sigM;
        nuc = dCritZ/sigM;
        tdn = A*sqrt(2.0*a/pi)*fabs(dlsdlM)/tM*nuc*(1.0+pow(nuc*nuc*a,-p));
        tdn *= exp(-a*nuc*nuc/2.0);
        tdn *= CRITDENMSOLMPC*self.Omegamhh;
        return tdn;


    def sigma0fM(self,tM,iDeriv=False):
        """Calculates and returns sigma(tM), and if iDeriv=1, calculates and stores
        * dsdM = d sigma(tM)/d tM.  Procedure comes from Rennan Barkana and is 
        * rather complex; could probably be improved.
        *   tM = halo mass (Msun)
        *   dsdM = place to store derivative (if iDeriv=0, ignored)
        """
        klimits=numpy.zeros(MAXBINS);
        dsdM=None

        #scale in comoving Mpc
        scale = pow(tM*3.0/(4.0*pi*CRITDENMSOLMPC*self.Omegamhh),1.0/3.0);
        self.scale=scale
        self.setSigmatop(0.0);
        self.setdsigmatop(0.0);
        acc = 2.0e-6;
        klimits[0] = min(1.0e-3/scale,0.2);
        klim = 2.0e2/scale;
        khi = klim;
        klhi = log(khi);
        if (klhi > klimits[0]):
            klim = min(klim,khi)
        nk = int(log(klim/klimits[0])/log(10.0)+0.5);
        if (nk < 5):
            nk = 5
        else:
            if (nk > 7):
                nk = 7
                
        kfac = log(klim/klimits[0])/(nk - 1.0);
        klimits[0] = log(klimits[0]);
        for i in range(1,nk):
            klimits[i] = klimits[i-1]+kfac
        klow = 1.0e-9;
        if (khi < klow):
            sigma0fm = 0.0;
        else: 
            if (klhi < klimits[0]):
                sigma0fm = scipy.integrate.quad(self.sigmatop2,klow,khi,epsrel=acc)[0]
                #cout << "In low sigma0fm =" << sigma0fm;
            else:
                sigma0fm = scipy.integrate.quad(self.sigmatop2,klow,exp(klimits[0]),epsrel=acc)[0]
                #cout << "In high, sigma0fm = " << sigma0fm << endl
                for j in range(1,nk):
                    sigma0fm += scipy.integrate.quad(self.sigmatop,klimits[j-1],klimits[j],epsrel=acc)[0]
                    #cout << "In high, sigma0fm = " << sigma0fm << endl;
        
        sigma0fm = self.sNorm*sqrt(sigma0fm);
        #cout << "Final sigm0fm = " << sigma0fm << endl;

        if iDeriv:
            if (khi < klow):
                dsigma0fm = 0.0;
            else:
                if (klhi < klimits[0]):
                    dsigma0fm = scipy.integrate.quad(self.dsigmatop2,klow,khi,epsrel=acc)[0];
                else:
                    dsigma0fm = scipy.integrate.quad(self.dsigmatop2,klow,exp(klimits[0]),epsrel=acc)[0];
                    for jj in range(1,nk):
                        dsigma0fm += scipy.integrate.quad(self.dsigmatop,klimits[jj-1],klimits[jj],epsrel=acc)[0];
	
            dsigma0fm += self.sigmatop(klhi)/scale;
            dsdM = 1.0/(pi*8.0*CRITDENMSOLMPC*self.Omegamhh*scale*scale);
            dsdM = self.sNorm*self.sNorm*dsdM*dsigma0fm/sigma0fm;

        if dsdM is not None:    
            return sigma0fm,dsdM
        else:
            return sigma0fm


    def setSigmatop(self,kl):
        """ Calculates sigmatophat, assuming k is entered as ln(k) """
        k = exp(kl);
        x = self.getScale()*k;
  
        sigtop = pow(k,3.0+self.nscal)*pow(self.TFMaster(k),2.0)*pow(3.0*(x*cos(x) - sin(x))/pow(x,3.0),2.0)
        #print x, k, sigtop
        return sigtop


    def sigmatop(self,kl):
        return self.setSigmatop(kl);

    def sigmatop2(self,k):
        """Calculates sigmatophat, assuming k is entered as linear. """
        return self.sigmatop(log(k))/k;


    def setdsigmatop(self,kl):
        """ Calculates derivative of sigmatophat, assuming k is entered as log(k).
        * Derivative calculation; used integration by parts to get a positive 
        * integrand """
        k = exp(kl);
        x = self.getScale()*k;

        t1,dTFdk = self.TFMaster(k,True);
        dsigmatop = -(pow(k,3.0+self.nscal)*t1*(2.0*dTFdk*k+t1*(3.0+self.nscal))*pow(3.0*(x*cos(x) - sin(x))/pow(x,3.0),2.0)/(self.getScale()));    
        return dsigmatop;


    def dsigmatop(self,kl):
        return self.setdsigmatop(kl);

    def dsigmatop2(self,k):
        """Calculates derivative of sigmatophat, assuming k is entered as linear."""
        return self.dsigmatop(log(k))/k;


    def sigm(self,m,flag=False):
        """Interpolation of sigma(m).  Returns sigma.  Slightly nonstandard, if:
        *     flag=None, initializes table
        *     flag=False, just calculates and returns sigma
        *     flag=True, calculates and returns sigma, also dsdm=derivative
        * Parameters are
        *     m = mass (Msun)
        *     dsdm = variable in which to store derivative; ignored if flag 
        *            does not equal 2
        * Uses spline interpolation from Numerical Recipes.
        """
        gridpts=300
        dsdm=None

        if (flag is None) or (self.sigma_s is None):
            # Initialize the spline fit
            print "initialising sigma spline"
            lgmass = numpy.zeros(gridpts);
            sigm = numpy.zeros(gridpts);
            sigsp = numpy.zeros(gridpts);
            sigderiv = numpy.zeros(gridpts);
            sigderivsp = numpy.zeros(gridpts);
            
            mMin = 1.0e3;
            mMax = 1.0e19;
            lgmMin = log(mMin);
            lgmMax = log(mMax);
            deltaM = (lgmMax-lgmMin)/float(gridpts);

            for i in range(gridpts):
                lgmass[i] = lgmMin + deltaM*i;
                mass = exp(lgmass[i]);
                sigm[i],deriv = self.sigma0fM(mass,True);
                sigderiv[i] = deriv;

            self.sigma_s=scipy.interpolate.UnivariateSpline(lgmass,sigm,s=0)
            self.sigma_sd=scipy.interpolate.UnivariateSpline(lgmass,sigderiv,s=0)
            print "sigma spline initialised"
            #return 0.0;
  
        lgmMin=self.sigma_s.get_knots()[0]
        lgmMax=self.sigma_s.get_knots()[-1]

        #This next sequence is not vectorised! :(
        #trouble lies in possibility of values outside range of spline
        #need to check bounds and handle exceptions separately
        #looks horrible and must be simpler way, especially not wild on
        #use of isinstance here...

        if isinstance(m,numpy.ndarray):
            #m is a numpy array so use vectorised version
            if (min(m) > exp(lgmMin) and max(m) < exp(lgmMax)):
                #good to use spline directly
                lgm = numpy.log(m);
                sig=self.sigma_s(lgm)
                if (flag == True):
                    dsdm=self.sigma_sd(lgm)
            else:
                #handle element by element to allow calc of elements outside spline
                raise Exception("not coded to handle vector with elements outside of spline")
                
        else:
            #m is just a number
            if (m > exp(lgmMin) and m < exp(lgmMax)):
                lgm = log(m);
                sig=self.sigma_s(lgm)
                if (flag == True):
                    dsdm=self.sigma_sd(lgm)
            else:
                sig,dsdm = self.sigma0fM(m,flag)

        if dsdm is None:
            return sig;
        else:
            return sig,dsdm

###########################################################
############## Halo Charactetics
#######################################################

    def jeansMass(self,z):
        """ Jeans mass; BL01 eq. 41 """
        return (6.2*pow(self.getOm0hh(),-0.5)*pow(self.getOmbhh(),-0.6)*pow(1.0+z,1.5));


    def filterMass(self,z):
        """ Filtered mass, assuming a continuous Jeans mass and adiabatic
        expansion past decoupling.  See Gnedin (2000) """
        zt = 137.0*pow(self.getOmbhh()/0.022,0.4) - 1.0;
        mass = 0.6*(1.0-2.0/3.0*sqrt((1.0+z)/(1.0+zt)));
        mass += log((1.0+zt)/(1.0+z))-2.0+2.0*sqrt((1.0+z)/(1.0+zt));
        mass = pow(3.0*mass,1.5);
        mass *= self.jeansMass(z);
        return mass;

    
    def minIonMass(self,z):
        """ Minimum mass for collapse assuming an ionized medium; condition is
        Tvir > 2x10^5 K.  Uses BL01 relations. """
        mass = 1.3e12/self.getH()/pow(1.0+z,1.5);
        mass *= sqrt(self.omegaZ(z)/self.getOmega0()/self.DeltaC(z));
        return mass;

    
    def coolMass(self,z):
        """ Analytic solution for when a halo hits 10^4 K, given z.  
        Uses BL01 relations. """
        m = 1.5e10/self.getH();
        m /= sqrt(self.getOmega0()/self.omegamZ(z)*self.DeltaC(z));
        m *= pow(1.0+z,-1.5);
        return m;


    def coolMassH2(self,z):
        """Analytic solution for when a halo hits 500 K, given z.  
        * Uses BL01 relations. """
        
        m = 1.69e8/self.getH();  #T=500K
        m /= sqrt(self.getOmega0()/self.omegamZ(z)*self.DeltaC(z));
        m *= pow(1.0+z,-1.5);
        return m;


    def rvir(self,mass,z):
        """ Returns virial radius in kpc; mass is in solar masses.  BL01 eq. 24 """
        ans = pow(mass*self.omegamZ(z)/self.getOmega0()/pow(self.getH(),2.0)/self.DeltaC(z),1.0/3.0);
        ans *= 9.5e-2/(1.0+z);
        return ans;


    def vcirc(self,mass,z):
        """ Returns circular velocity in km/s; mass is in solar masses. BL01 eq. 25 """
        ans = pow(mass*mass/self.omegamZ(z)*self.getOmega0()*pow(self.getH(),2.0)*self.DeltaC(z),1.0/6.0);
        ans *= 6.72e-3*sqrt(1.0+z);
        return ans;


    def tvir(self,mass,z,mu):
        """ Returns virial temperature in K; mass is in solar masses and mu
 * is mean molecular weight.  BL01 eq. 26 """
        ans = pow(mass*mass/self.omegamZ(z)*self.getOmega0()*pow(self.getH(),2.0)*self.DeltaC(z),1.0/3.0);
        ans *= 2.72e-3*mu*(1.0+z);
        return ans;


    def RComfromM(self,m):
        """ Convert a mass to comoving size """
        R = pow(3.0/4.0/pi*m/CRITDENMSOLMPC/self.getOmega0hh(),1.0/3.0);
        return R;


    def MfromRCom(self,R):
        """ Convert a comoving size to a mass """
        m = 4.0*pi/3.0*R*R*R*CRITDENMSOLMPC*self.getOmega0hh();
        return m;


