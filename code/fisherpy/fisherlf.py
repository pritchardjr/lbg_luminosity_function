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
import luminosityfunction as CLF

################################################################
# Define class for Galaxy survey parameters
############################################################

class GalaxySurvey:
    """ Class for galaxy survey and its properties
    angles in degrees

    NOTE: not entirely sure survey should own its own cosmology. May be better to separate survey fixed quantities from translation to sky
    """
    def __init__(self,name="GalaxySurvey",zmin=1.0,zmax=2.0,anglex=1.0,angley=1.0,nfield=1,nbin=1,maglim=-13.0):
        #cosm=CCosm.Cosmology()
        self.name=name
        self.zmin=zmin
        self.zmax=zmax
        self.anglex=anglex
        self.angley=angley
        self.maglim=maglim
        self.nfield=nfield
        self.nbin=nbin

        self.z=(self.zmin+self.zmax)/2.0

    def __repr__(self):
        return "%s z=[%g,%g]; theta=[%g,%g]; nfield=%g; nbin=%g; maglim=%g" % (self.name,self.zmin,self.zmax,self.anglex,self.angley,self.nfield,self.nbin,self.maglim)

    def area(self):
        """ survey area in sq arcmin """
        return (self.anglex*60.0)*(self.angley*60.0)

    def areaSTER(self):
        """ survey area in steradians """
        return self.area()*(pi/180.0/60.0)**2

    def volume(self,cosm):
        """ survey volume in Mpc^3 """
        vol=cosm.comovingVolume(z1=self.zmin,z2=self.zmax)
        vol*=self.areaSTER()/(4.0*pi)*self.nfield;
        return vol


############################################################
#Create Fisher matrix for galaxy LF constraints following Robertson (2010)
############################################################
class FisherLF(Fisher):
    """ Fisher class for Galaxy Luminosity Function derived calculations """
    def __init__(self,paramdict={'oml':0.7,'hubble':0.7,'omm':0.3,'omk':0.0,'ombhh':0.0225,'nscal':0.95,'ascal':25.0,'sig8':0.9},survey=GalaxySurvey()):
        Fisher.__init__(self,paramdict)
        self._fishertype="FisherCL"
        self.survey=survey
        self.massfcn='TK'

        #calculate fiducial cosmology
        cosm=CCosm.Cosmology(paramdict=paramdict)
        self.cosm=cosm

        #set fiducial luminosity function
        self.LF=CLF.LuminosityFunction()

        #use fiducial cosmology to map survey parameters to comoving distances
        self.survey_r=cosm.comovingDistance(survey.z)
        self.survey_deltar=abs(cosm.comovingDistance(survey.zmax)-cosm.comovingDistance(survey.zmin))
        self.survey_vol=survey.volume(cosm)

        #calculate magnitude bins to use
        # take faintest magnitude from galaxy survey
        # take brightest magnitude to be arbitrarily high
        self.deltaMag=0.25
        self.magmax=survey.maglim #faintest objects in survey
        magmin= -29.0 #brightest possible
        self.Nbin= int((self.magmax-magmin)/self.deltaMag)
        self.magbins=self.magmax-np.arange(self.Nbin)*self.deltaMag
        self.magmin=self.magbins[-1]

        # Work around slow numerical integral that's cosmology dependent
        self.nasty=None

    def calcFisherMatrix(self):
        #self._fisher=np.identity(self._nparam)
        
        #First calculate the covariance matrix
        D=np.zeros([self.Nbin,self.Nbin])
        for i in range(self.Nbin):
            for j in range(self.Nbin):
                Cij=self.covarianceMatrix(i,j,self.cosm,self.LF)
                Pij=self.poissoncov(i,j,self.cosm,self.LF)
                D[i][j]=Cij+Pij

        invD=np.linalg.inv(D)

        #Now calculate the different terms of the Fisher mtrix
        
                

    def covarianceMatrix(self,i,j,cosm,LF):

        deltaMag=self.deltaMag
        Mi=self.magbins[i]
        Mj=self.magbins[j]
        ni=LF.numberdensityInMBin(Mi,deltaMag)
        nj=LF.numberdensityInMBin(Mj,deltaMag)
        mimin=LF.halomatch(Mi-deltaMag/2.0,cosm,self.massfcn)
        mimax=LF.halomatch(Mi+deltaMag/2.0,cosm,self.massfcn)
        mjmin=LF.halomatch(Mj-deltaMag/2.0,cosm,self.massfcn)
        mjmax=LF.halomatch(Mj+deltaMag/2.0,cosm,self.massfcn)

        print Mi,Mj,mimin,mimax,mjmin,mjmax

        bi=cosm.biasInBin(self.survey.z,mimin,mimax,self.massfcn)
        bj=cosm.biasInBin(self.survey.z,mjmin,mjmax,self.massfcn)
        
        #now caluclate the covariance matrix
        Sij=bi*bj*ni*nj/self.survey.nfield
        #Sij*=pow(cosm.growthFac(self.survey.z),2.0) - included in P(z,k)

        #now return 3D integral over P(k)
        if self.nasty==None:
            integral=self.nastyintegral(cosm,LF)
            self.nasty=integral
        else:
            integral=self.nasty
        print integral

        Sij*=integral[0]

        return Sij

    def poissoncov(self,i,j,cosm,LF):
        """ Poisson contribution to the covariance"""
        if i==j:
            Vi=self.survey_vol
            Pij=ni=LF.numberdensityInMBin(self.magbins[i],self.deltaMag)/Vi
        else:
            Pij=0.0
        return Pij

    def window(self,x,y,z):
        """ Survey window function in Fourier space

        Assumes Cartesian geometry, which would break down for a large
        area galaxy survey. Integrity applied for  angle<1 deg or delta_r/r<<1
        """
        def sinc(x):
            if x == 0:
                return 0
            else:
                return sin(x)/x
        WV=sinc(x)
        WV*=sinc(y)
        WV*=sinc(z)
        return WV

    def nastyintegral(self,cosm,LF):
        Rx=self.survey_r*self.survey.anglex/2.0
        Ry=self.survey_r*self.survey.angley/2.0
        Rz=self.survey_deltar/2.0

        #window function rapidlly kills kernal on scales > 1/R
        mymax=1.0e1
        xmin,xmax=-mymax/Rx,mymax/Rx
        ymin,ymax=-mymax/Ry,mymax/Ry
        zmin,zmax=-mymax/Rz,mymax/Rz
        print xmin,xmax,ymin,ymax,zmin,zmax

        result=scipy.integrate.tplquad(self.nastyKernal,xmin,xmax,lambda kx: ymin, lambda kx: ymax, lambda kx,ky: zmin, lambda kx,ky: zmax,args=(cosm,LF,Rx,Ry,Rz))

        return result

    def nastyKernal(self,kz,ky,kx,cosm,LF,Rx,Ry,Rz):
        """ Horrible 3D integral of window function over power spectrum"""

        kernal=pow(self.window(kx*Rx,ky*Ry,kz*Rz),2)
        k=sqrt(kx*kx+ky*ky+kz*kz)
        kernal*=cosm.linearPowerSpectrum(k,self.survey.z)
        kernal/=pow(2.0*pi,3)
        
        return kernal
        
############################################################
    """
    Functions to add

    Halo matching: nM from LF to halo mass

    Mass averaged bias
    
    """
