"""
Calculation of the uvdensity for an array

This is experiment focussed and can calculate the uv distribution
but also the power spectrum errors that one would expect
"""
from math import *
import numpy as np
import cosmologypy.Cosmology as Cosm
import cosmoconstants as CC
import scipy.integrate
import matplotlib.pyplot as plt

class Experiment:
    def __init__(self,Atot=6e4,Nant=1400,Dmax=1e3,B=8.0,z=7.5,tint=1000,cosm=Cosm.Cosmology(),eta=0.5,FoV=pi*8.6**2):
        self.Nant=Nant
        self.Dmax=Dmax
        
        self.Aeff=Atot/Nant
        self.tint=tint*CC.hour
        self.z=z
        self.B=B*1e6 #convert MHz to Hz
        self.eta=eta
        self.wavelength=CC.lambda21*(1.0+self.z)
        self.chi=cosm.comovingDistance(self.z)
        self.y=1.0/CC.nu21*(1+self.z)**2/cosm.hubbleZMPC(self.z)
        self.FoV=self.wavelength**2/self.Aeff
        self.FoV=FoV*CC.deg2steradians
        self.volume=self.surveyVolume()
        self.cosm=cosm

        filling=0.8
        powerlaw= 2.0
        Rmin=10.0
        self.configArrayPowerLaw(Rmin,Nant,filling,powerlaw)

    def expsummary(self):
        expdict={}
        expdict['x']=self.chi
        expdict['y']=self.y
        expdict['B']=self.B
        expdict['tint']=self.tint
        expdict['eta']=self.eta
        expdict['Aeff']=self.Aeff
        expdict['lambda']=self.wavelength
        expdict['z']=self.z
        expdict['Dmax']=self.Dmax
        expdict['Nant']=self.Nant
        if self.Rcore is not None: expdict['Rcore']=self.Rcore
        if self.Rnucleus is not None: expdict['Rnucleus']=self.Rnucleus
        for key in expdict.keys():
            print key, expdict[key]
        return

########## Stock experiment parameters
    def stockExperiments(self,case="ska"):

        #here taking values from Joudaki+ foor testing purposes
        if case=="ska":
            self.__init__(Atot=6e4,Nant=1400,Dmax=1e3,B=8.0,z=7.5,tint=1000,cosm=Cosm.Cosmology(),eta=0.5,FoV=pi*8.6**2)
        elif case=="lofar":
            self.__init__(Atot=32.0*590.0,Nant=32,Dmax=1e3,B=8.0,z=7.5,tint=1000,cosm=Cosm.Cosmology(),eta=0.5,FoV=2.0*pi*2.4**2)
        elif case=="mwa":
            self.__init__(Atot=500*13.0,Nant=500,Dmax=1e3,B=8.0,z=7.5,tint=1000,cosm=Cosm.Cosmology(),eta=0.5,FoV=pi*16.0**2)

########## UV coverage: Analytic results
    def uvdensityUniform(self,b,wavelength):
        """ Uniform UV coverage"""
        if b<self.Dmax:
            nbaseline=self.Nant**2/pi/self.Dmax**2
            return wavelength**2*nbaseline
        else:
            return 0

    def uvdensityFlat(self,b,wavelength):
        """ Analytic form for the baseline distribution given uniform antennae
        distribution of Nant out to Dmax"""

        a=self.Dmax/2.0
        Asec=(np.arccos(b/2.0/a))   #use identity that arcsec(x)=arccos(1/x)
        F= -b*np.sqrt(4.0*a*a-b*b)/a/a+4*Asec
        F/=(2.0*pi)
        F=np.nan_to_num(F) #handles baselines beyond 2*Dmax
        #print F

        #Alternative but equivalent trig form
        #J=b*np.sqrt(4.0*a*a-b*b)/a/a
        #J+=2*(np.arcsin(b/2.0/a)+np.arctan2(1,np.sqrt(4*(a/b)**2-1.0)))
        #J/=(2.0*pi)
        #J=np.nan_to_num(J)
        #print 1-J

        nbaseline=4.0*self.Nant**2/self.Dmax**2*F
        #nbaseline/=4.0*pi
        nbaseline/=2.0*pi
        return wavelength**2*nbaseline

################ UV coverage: Numerical calculation
    def configArrayPowerLaw(self,Rmin,Nant,filling,powerlaw= 2.0):
        """ calculate array configuration using the formula from
        Mao+
        """

         #Nucleus of 100% filled i.e. at maximal density
        rho0=1/self.Aeff
        Rnucleus=sqrt(filling*Nant/(rho0*pi))

        #Around nucleus have a core where rho~r^-powerlaw
        if powerlaw==2.0:
            Rcore=Rnucleus*exp((1-filling)/(2.0*filling))
        else:
            Rcore=Rnucleus*pow((2.0-powerlaw*(1.0-filling))/(2.0*filling),1.0/(2.0-powerlaw))

        #Mao includes an annulus around that, but ignore that here

        self.Rmin=Rmin
        self.Rnucleus=Rnucleus
        self.Rcore=Rcore
        self.Dmax=2.0*self.Rcore
        

    def uvdensityGeneral(self,b,wavelength):
        rmin=self.Rmin
        rmax=self.Dmax/2.0
        phimin=0.0
        phimax=2.0*pi

        #most simplistic version that works, but is slow
        nbaseline,error=scipy.integrate.dblquad(self.kernelCorrelation,rmin,rmax,lambda phi: phimin, lambda phi: phimax,args=(b,),)
        print nbaseline

        #slightly smarter version exploiting symmetry of integral
        #nbaseline,error=scipy.integrate.dblquad(self.kernelCorrelation,rmin,rmax,lambda phi: acos((rmax**2-b**2-r**2)/(2*r*b)), lambda phi: 2.0*pi,args=(b,),)
        #nbaseline*=2.0
        print nbaseline
        
        print b,nbaseline
        return wavelength**2*nbaseline
    

    def kernelCorrelation(self,phi,r,b):
        rb=sqrt(r*r+b*b+2.0*r*b*cos(phi))
        kernel=r*self.rho(r)*self.rho(rb)/2.0
        return kernel
    

    def rho(self,r):
        """ Filled nucleus + powerlaw core """
        #rho0=self.Nant/(pi*(self.Dmax/2.0)**2)        
        #Rcore=self.Dmax/2.0

        rho0=1.0/self.Aeff
        Rnucleus=self.Rnucleus

        if r<Rnucleus:
            return rho0
        elif r>self.Dmax/2.0:
            return 0
        else:
            return rho0*(r/Rnucleus)**-2

#################### Other things
    def surveyVolume(self):
        volume=self.chi**2*self.y*self.B*self.FoV        
        return volume
    
    def Tsky(self,wavelength):
        return 60.0*np.power(wavelength,2.55)

    def noisepower(self,uperp):
        """ Noise power spectrum as function of uperp
        would be nice to be able to feed this different uvdensity calculations
        as needed
        """
        nbaseline=self.uvdensityFlat(uperp,self.wavelength)
        #nbaseline=self.uvdensityUniform(uperp,self.wavelength)
        Tsys=self.Tsky(self.wavelength)
        Pnoise=(self.wavelength**2*Tsys/self.Aeff)**2
        Pnoise/=self.tint*nbaseline
        return Pnoise

    ######## K boundaries
    def maximumK(self):
        kmax=1.0
        return kmax

    def minimumK(self):
        kmin=0.002
        return kmin

    def fisherKStar(self):
        """ maximum perp k value in Mpc """
        kstar=2.0*pi/(self.wavelength/self.Dmax*self.chi)
        return kstar
    
    def thetaminK(self,k): 
        thetamin=np.arccos(np.minimum(self.y*self.B*k/(2.0*pi),1.0))
        return thetamin

    def thetamaxK(self,k):
        kstar=self.fisherKStar()
        thetamax=np.arcsin(np.minimum(kstar/k,1.0))
        return thetamax

    ####### u to k conversions
    def baselineFromU(self,uperp):
        return uperp*self.wavelength/(2*pi)
        
    def uperpFromK(self,kperp):
        return kperp*self.chi

    def uparaFromK(self,kpara):
        return kpara*self.y
    
    def kperpFromU(self,uperp):
        return uperp/self.chi

    def kparaFromU(self,upara):
        return upara/self.y

    def nuFromZ(self,z):
        """ Frequency corresponding to 21cm at redshift z """
        nu=CC.nu21/(1.0+z)
        return nu

    def zFromNu(self,nu):
        """ Redshift corresponding to 21cm observed at nu """
        z=CC.nu21/nu -1.0
        return z
        
    ######## Power spectrum errors
    def brightnessT(self,z,omegaHI=1.0e-3,cosm=None):
        """ Intensity mapping brightness tmperature Kelvin """
        if cosm is None:
            cosm=self.cosm
        
        tb=0.3*(omegaHI/1.0e-3)
        tb/=sqrt((cosm.Omegam+cosm.Omegal/(1.0+z)**3)/0.29)
        tb*=np.sqrt((1.0+z)/2.5)
        return tb/1.0e3

    
    def cosmicpower(self,uperp,upara):
        kpara=self.kparaFromU(upara)
        kperp=self.kperpFromU(uperp)    
        k=np.sqrt(kpara**2+kperp**2)
        pk=self.cosmicpowerK(k)
        pu=pk/(self.chi**2*self.y)
        return pu

    def cosmicpowerK(self,k):
        pk=self.cosm.linearPowerSpectrum(k,self.z)
        tb=0.02
        pk*=tb**2
        return pk

    def weightVariance(self,uperp,upara):
        """ delta P(u)^2 for Fisher matrix calc in uperp,upara cell
        missing a d^3u factor deliberately
        """
        pcosm=self.cosmicpower(uperp,upara)
        pcosm=0
        pnoise=self.noisepower(self.baselineFromU(uperp))
        ptot=pcosm+pnoise
        #print pcosm,pnoise,ptot
        Nc=2.0*pi*self.volume/(2.0*pi)**3
        weight=Nc/ptot
        #print ptot, weight
        return weight

    def weightVarianceK(self,theta,k):
        """ delta P(k)^2 for Fisher matrix calc in dk,dtheta cell"""
        kpara=k*cos(theta)
        kperp=k*sin(theta)
        upara=self.uparaFromK(kpara)
        uperp=self.uperpFromK(kperp)
        weight=self.weightVariance(uperp,upara)
        #convert from 1/P(u)^2 to 1/P(k)^2
        weight/=(self.chi**2*self.y)**2
        #print k,theta,upara,uperp,weight
        return weight*k*k*sin(theta)

    def powerErrorAveraged(self,k):
        """ calculate the error on the spherically averaged power spectrum
        carried out in uv-space
        """

        kstar=self.fisherKStar()
        #thetamin=np.arccos(np.minimum(self.y*self.B*k/(2.0*pi),1.0))
        #thetamax=np.arcsin(np.minimum(kstar/k,1.0))
        thetamin=self.thetaminK(k)
        thetamax=self.thetamaxK(k)

        #print k,thetamin,thetamax,self.y*self.B*k/(2.0*pi)
        
        deltaP,error=scipy.integrate.quad(self.weightVarianceK,thetamin,thetamax,args=(k))
        deltaP*=self.eta*k #multiply by d(log(k))
        print deltaP
        if deltaP>0:
            deltaP=1.0/sqrt(deltaP)
        else:
            deltaP=0.0
        return deltaP

    def testPlot(self):
        ek=np.linspace(-2,1)
        k=np.power(10,ek)
        pk=[self.powerErrorAveraged(kk) for kk in k]
        deltak=pk*k**3/2/pi**2*1e6 # in mK^2
        plt.loglog(k,deltak)
        pk=[self.cosmicpowerK(kk) for kk in k]
        deltakcosm=pk*k**3/2/pi**2*1e6 # in mK^2
        plt.loglog(k,deltakcosm)
        plt.show()
        return k,pk

    def testPlotUV(self):
        ek=np.linspace(0,4)
        u=np.power(10,ek)
        nbaseline=[self.uvdensityFlat(self.baselineFromU(kk),self.wavelength) for kk in u]
        plt.loglog(u,nbaseline)
        plt.show()

    def testPlotN(self):
        ek=np.linspace(0,log10(0.6*self.Dmax))
        b=np.power(10,ek)
        nbaselineF=[self.uvdensityFlat(bb,self.wavelength) for bb in b]      
        plt.loglog(b,nbaselineF)

        nbaseline=[self.uvdensityGeneral(bb,self.wavelength) for bb in b]      
        plt.loglog(b,nbaseline)

        plt.show()

        return nbaseline,nbaselineF


    def testPlotRhoAntennae(self):
        """ Plot rho distribution """
        low=log10(max(self.Rmin,1))  #plot from 1m or Rmin
        ek=np.linspace(low,log10(self.Dmax))
        r=np.power(10,ek)
        rho=[]
        for ri in r:
            rho.append(self.rho(ri))
        rho=np.array(rho)

        plt.loglog(r,rho)
        plt.show

        return (r,rho)
        

