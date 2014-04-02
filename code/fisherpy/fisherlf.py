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

        self.z=(self.zmin+self.zmax)/2.0

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
#Create Fisher matrix for galaxy LF constraints following Robertson (2010)
############################################################
class FisherLF(Fisher):
    """ Fisher class for Galaxy Luminosity Function derived calculations """
    def __init__(self,paramdict={'oml':0.7,'hubble':0.7,'omm':0.3,'omk':0.0,'ombhh':0.0225,'nscal':0.95,'ascal':25.0},survey=GalaxySurvey()):
        Fisher.__init__(self,paramdict)
        self._fishertype="FisherCL"
        self.survey=survey

        #calculate fiducial cosmology
        cosm=CCosm.Cosmology()

        #use fiducial cosmology to map survey parameters to comoving distances
        self.survey_r=cosm.comovingDistance(survey.z)
        self.survey_deltar=abs(cosm.comovingDistance(survey.zmax)-cosm.comovingDistance(survey.zmin))

    def calcFisherMatrix(self):
        #self._fisher=np.identity(self._nparam)
        raise Exception("calcFisherMatrix not defined")

    def covarianceMatrix(i,j,cosm):
        
        Sij=bi*bj*ni*nj/self.survey.nfield
        Sij*=pow(cosm.growthFac(z),2.0)

        #now return 3D integral over P(k)
        

        return Sij

    def window(self,kx,ky,kz):
        sinc=lambda x: sin(x)/x
        WV=sinc(kx*self.survey_r*self.survey.anglex/2.0)
        WV*=sinc(ky*self.survey_r*self.survey.angley/2.0)
        WV*=sinc(kz*self.survey_deltar/2.0)

        
############################################################
    """
    Functions to add

    Halo matching: nM from LF to halo mass

    Mass averaged bias
    
    """
