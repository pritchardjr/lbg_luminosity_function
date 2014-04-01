#
# Hurried code to calculate an intensity mapping Fisher matrix
#
# This works in k-space
#
#

from Cosmology import *
import tbconstants
import cosmoconstants
import numpy
import math
import scipy.integrate

try:
    import cmbpyJP as cmbpy
except:
    print 'cmbpyJP not available'

class tbfisher(Cosmology):
    """ Fisher matrix built on top of fiducial cosmology """
    def __init__(self):
        Cosmology.__init__(self)
        self.setExperiment()

    def setExperiment(self,z=1.5,B=400.0,Aeff=10.0,tint=5000.0,Nant=1000,Dmax=100.0,eta=0.01):
        """ Establish experimental parameters """
        self.z=z
        self.B=B*1.0e6
        self.Aeff=Aeff
        self.tint=tint*cosmoconstants.hour
        self.Nant=Nant
        self.Dmax=Dmax
        self.eta=eta
        self.wavelength=tbconstants.lambda21*(1.0+self.z)

        #Derived parameters
        self.x=self.comovingDistance(self.z)
        self.y=self.bandwidthDistance(self.z,self.B)
        
        return

    def expsummary(self):
        expdict={}
        expdict['x']=self.x
        expdict['y']=self.y
        expdict['B']=self.B
        expdict['tint']=self.tint
        expdict['eta']=self.eta
        expdict['Aeff']=self.Aeff
        expdict['lambda']=self.wavelength
        expdict['z']=self.z
        expdict['Dmax']=self.Dmax
        for key in expdict.keys():
            print key, expdict[key]
        return

    ########## Brightness temperature
    def brightnessT(self,z,omegaHI=1.0e-3,cosm=Cosmology()):
        """ Intensity mapping brightness tmperature Kelvin """
        tb=0.3*(omegaHI/1.0e-3)
        tb/=math.sqrt((cosm.Omegam+cosm.Omegal/(1.0+z)**3)/0.29)
        tb*-numpy.sqrt((1.0+z)/2.5)
        return tb/1.0e3

    ########## Specific distances
    def nuFromZ(self,z):
        """ Frequency corresponding to 21cm at redshift z """
        nu=tbconstants.nu21/(1.0+z)
        return nu

    def zFromNu(self,nu):
        """ Redshift corresponding to 21cm observed at nu """
        z=tbconstants.nu21/nu -1.0
        return z
    
    def bandwidthDistance(self,z,B):
        """ comoving distance covered by bandwidth B [Hz] at z """
        nu0=self.nuFromZ(z)
        z1=self.zFromNu(nu0-B/2)
        z2=self.zFromNu(nu0+B/2)
        distance=self.comovingDistance(z1)-self.comovingDistance(z2)
        return distance

    ######### Experiment stuff

    def Tsys(self,z):
        """ System temperature in Kelvin """
        tsys=50.0
        return tsys
        
    def visibilityDensity(self,uperp):
        nvis=self.wavelength**2*self.Nant**2/self.Dmax**2/math.pi
        return nvis

    ######## K boundaries
    def maximumK(self):
        kmax=1.0
        return kmax

    def minimumK(self):
        kmin=0.002
        return kmin
    
    ######## Fisher matrix stuff

    def fisherE(self):
        """ McQuinn E parameter """
        EE=2.0*math.pi*math.sqrt(self.x**2*self.y)
        EE*=self.wavelength**3
        EE*=self.Tsys(self.z)**2
        EE/=self.Aeff**1.5*self.B*self.tint*math.sqrt(self.eta)
        return EE

    def fisherD(self):
        """ McQuinn D parameter """
        DD=(2.0*math.pi)**2*self.Aeff
        DD/=self.wavelength**2*self.x*self.x*self.y
        return math.sqrt(DD)

    def fisherKStar(self):
        """ maximum perp k value in Mpc """
        kstar=2.0*math.pi/(self.wavelength/self.Dmax*self.x)
        return kstar

    def powerErrorUniform(self,k):
        """ Power spectrum errors assuming uniform visibility coverage
        McQuinn (B3)
        """
        kstar=self.fisherKStar()
        if k<=kstar:
            error=min(self.y*k/(2.0*math.pi),1.0)**-0.5
        elif k>kstar:
            error=(1-math.sqrt(k**2-kstar**2)/k)**-0.5

        rho=self.visibilityDensity(0.0)
        error*=self.fisherE()/k**1.5/rho
                
        return error

    ########## That was a warm up. Now the angular integral to check errors

    def weightVariance(self,theta,k):
        """ really the weighting for the Fisher integral """
        Pk=0.0
        kperp=k*math.sin(theta)
        weight=self.fisherD()*Pk+self.fisherE()/self.visibilityDensity(kperp)
        return math.sin(theta)/weight**2

    def powerErrorAveraged(self,k):
        """ McQuinn (B1) error on angle averaged power spectrum """
        kstar=self.fisherKStar()
        thetamin=math.acos(min(self.y*k/(2.0*math.pi),1.0))
        thetamax=math.asin(min(kstar/k,1.0))
        
        II,error=scipy.integrate.quad(self.weightVariance,thetamin,thetamax,args=(k))
        return 1.0/math.sqrt(II*k**3)
        
    ######## Getting there: now for a fisher matrix element
    def weightFisher(self,theta,k,Pk,dP1,dP2):
        """ really the weighting for the Fisher integral """
        #Pk=0.0
        mu=math.cos(theta)
        f=self.fGrowth(self.z)
        angular=(1.0+f*mu*mu)**2
        kperp=k*math.sin(theta)
        weight=self.fisherD()*Pk*angular+self.fisherE()/self.visibilityDensity(kperp)
        #print k, self.fisherD()*Pk, self.fisherE()/self.visibilityDensity(kperp)
        return k**3*math.sin(theta)*dP1*dP2*angular**2/weight**2

    
    def fisherElement(self,kvec,Pk,dP1,dP2):
        """
        taking vector of power spectrum derivatives dP1 and dP2
        calculate the amplitude of the Fisher matrix element
        
        """
        kstar=self.fisherKStar()
        fisher12=0.0
        self.eta=abs(math.log(kvec[0])-math.log(kvec[1]))
        kmin=self.minimumK()
        kmax=self.maximumK()

        for (i,k) in enumerate(kvec):
            #print i,k, Pk[i]
            if k>kmin and k<kmax:
                thetamin=math.acos(min(self.y*k/(2.0*math.pi),1.0))
                thetamax=math.asin(min(kstar/k,1.0))
                II,error=scipy.integrate.quad(self.weightFisher,thetamin,thetamax,args=(k,Pk[i],dP1[i],dP2[i]))
                fisher12+=II

        return fisher12

######## Calculate power spectrum derivatives
    
    def powerDerivative(self,tagdict,cosm):
        """ using same tags as CAMB access calculate derivatives
        tagdict should have labels and stepsize tagdict[tag]=step

        My best accuracy on derivatives seems to be fractional error of 0.002
        Worst is 0.02 on omk. This is comparable with the accuracy of CAMB
        which is quoted as being at the level of 0.1-0.2%
        """

        resultdict={}

        for tag in tagdict.keys():
            step=tagdict[tag]/2.0

            #Establish base values
            mydict=cmbpy.cambio.defaultCosmDict()
            mydict['output_root']='test'
            mydict['get_transfer']='T'
            mydict['ombh2']=cosm.Omegabhh
            mydict['omch2']=cosm.Omegachh
            mydict['hubble']=cosm.hubble*100.0
            mydict['w']= cosm.w
            mydict['omk']=cosm.Omegak
            mydict['scalar_amp(1)']=cosm.Ascal
            mydict['scalar_spectral_index(1)']=cosm.nscal
            mydict['transfer_redshift(1)']=cosm.z
            mydict['transfer_k_per_logint'] = 5
            mydict['transfer_high_precision'] = 'T'

            #Now derivative
            plusdict=mydict.copy()
            minusdict=mydict.copy()
            plusdict['output_root']=tag+'_p'
            minusdict['output_root']=tag+'_m'
            plusdict[tag]=mydict[tag]+step
            minusdict[tag]=mydict[tag]-step
            cmbpy.cambio.runCAMBFromDict(plusdict)
            cmbpy.cambio.runCAMBFromDict(minusdict)

            #double plus
            pplusdict=mydict.copy()
            mminusdict=mydict.copy()
            pplusdict['output_root']=tag+'_pp'
            mminusdict['output_root']=tag+'_mm'
            pplusdict[tag]=mydict[tag]+step*2.0
            mminusdict[tag]=mydict[tag]-step*2.0
            cmbpy.cambio.runCAMBFromDict(pplusdict)
            cmbpy.cambio.runCAMBFromDict(mminusdict)

            #From the output power spectra calculate the derivative
            plusfile=plusdict['output_root']+'_matterpower.dat'
            minusfile=minusdict['output_root']+'_matterpower.dat'
            dataP=numpy.loadtxt(plusfile)
            dataM=numpy.loadtxt(minusfile)

            #double plus power spectra
            pplusfile=pplusdict['output_root']+'_matterpower.dat'
            mminusfile=mminusdict['output_root']+'_matterpower.dat'
            dataPP=numpy.loadtxt(pplusfile)
            dataMM=numpy.loadtxt(mminusfile)

            #For some reason can't get CAMB to always output same number of
            # steps. Step size is fixed but number varies by 1 or 2
            #
            # NOTE: doesn't yet account for CAMB output being in units of
            # [h/Mpc] and [(Mpc/h)^3]
            minlength=min(len(dataP),len(dataM),len(dataPP),len(dataMM))
            usedict=plusdict
            cosmP=Cosmology(Omegabhh=usedict['ombh2'],Omegak=usedict['omk'],hubble=usedict['hubble']/100.0,nscal=usedict['scalar_spectral_index(1)'],Ascal=usedict['scalar_amp(1)']*1.0e10,Omegam=((usedict['ombh2']+usedict['omch2'])/(usedict['hubble']/100.0)**2),Omegal=1.0-((usedict['ombh2']+usedict['omch2'])/(usedict['hubble']/100.0)**2)-usedict['omk'])
            usedict=minusdict
            cosmM=Cosmology(Omegabhh=usedict['ombh2'],Omegak=usedict['omk'],hubble=usedict['hubble']/100.0,nscal=usedict['scalar_spectral_index(1)'],Ascal=usedict['scalar_amp(1)']*1.0e10,Omegam=((usedict['ombh2']+usedict['omch2'])/(usedict['hubble']/100.0)**2),Omegal=1.0-((usedict['ombh2']+usedict['omch2'])/(usedict['hubble']/100.0)**2)-usedict['omk'])
            TbP=self.brightnessT(self.z,cosm=cosmP)
            TbM=self.brightnessT(self.z,cosm=cosmM)
            k=dataP[0:minlength,0]
            PkP=dataP[0:minlength,1]*TbP**2
            PkM=dataM[0:minlength,1]*TbM**2
            dP=(PkP-PkM)/(2.0*step)

            # double plus
            usedict=pplusdict
            cosmPP=Cosmology(Omegabhh=usedict['ombh2'],Omegak=usedict['omk'],hubble=usedict['hubble']/100.0,nscal=usedict['scalar_spectral_index(1)'],Ascal=usedict['scalar_amp(1)']*1.0e10,Omegam=((usedict['ombh2']+usedict['omch2'])/(usedict['hubble']/100.0)**2),Omegal=1.0-((usedict['ombh2']+usedict['omch2'])/(usedict['hubble']/100.0)**2)-usedict['omk'])
            usedict=mminusdict
            cosmMM=Cosmology(Omegabhh=usedict['ombh2'],Omegak=usedict['omk'],hubble=usedict['hubble']/100.0,nscal=usedict['scalar_spectral_index(1)'],Ascal=usedict['scalar_amp(1)']*1.0e10,Omegam=((usedict['ombh2']+usedict['omch2'])/(usedict['hubble']/100.0)**2),Omegal=1.0-((usedict['ombh2']+usedict['omch2'])/(usedict['hubble']/100.0)**2)-usedict['omk'])
            TbPP=self.brightnessT(self.z,cosm=cosmPP)
            TbMM=self.brightnessT(self.z,cosm=cosmMM)
            k=dataP[0:minlength,0]
            PkPP=dataPP[0:minlength,1]*TbPP**2
            PkMM=dataMM[0:minlength,1]*TbMM**2

            dPP=(PkP-PkM)*8.0/12.0/step-(PkPP-PkMM)/12.0/step

            resultdict[tag]=(mydict[tag],step,k,dPP,dP)
            
        return resultdict

    def calculateFisher(self,tagdict):
        """ calculate the Fisher matrix

        NEEDS fiducial value for Pk
        """
        N=len(tagdict.keys())
        fisherMatrix=np.zeros((N,N))
        for (i,tag1) in enumerate(sorted(tagdict.keys())):
            kvec=dP1=resultdict[tag1][2]
            for (j,tag2) in enumerate(sorted(tagdict.keys())):
                print i, j, tag1, tag2
                dP1=resultdict[tag1][3]
                dP2=resultdict[tag2][3]
                fisherMatrix[i,j]=tocm.fisherElement(kvec,Pk,dP1,dP2)

        return fisherMatrix

    def printConstraints(self,fisherMatrix):
        inverseFisher=numpy.linalg.inv(fisherMatrix)
        constraints=numpy.sqrt(inverseFisher.diagonal())

        return inverseFisher, constraints
