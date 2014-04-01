#
# EoR Fisher matrix calculation
#
#

import cosmologypy as cosmopy
import cosmologypy.Cosmology as CC
import cosmologypy.uvdensity as UVD
import tbconstants
import cosmoconstants
import numpy as np
import math
import scipy.integrate
from fisherpy.fisher import Fisher
import cmbpyJP as cmbpy
import os

#################################################################
# Fisher class definition
#################################################################
class FisherTB(Fisher):
    """ 21cm Brightness temperature Fisher matrix code
    Specify:
    parameters to vary: paramdict
    Experiment to use: Exp
    Fiducial cosmology: cosm
    """
    def __init__(self,paramdict={},Exp=UVD.Experiment()):
        """
        Need to decide on how to specify cosmology and parameters
        in best possible way
        """

        #currently super confused cosmology/parameter specification
        self.cosm=CC.Cosmology()

        if paramdict == {}:
            paramdict=self.defaultParamDict(self.cosm)
        
        #assign the radio experiment properties to be used
        self.Exp=Exp

        #Work out the fiducial cosmology from paramdict and standard values
        self.cosm=self.assignFiducialCosmology(paramdict)

        #path to precalculated power spectrum files from CAMB
        self.datadir=None

        Fisher.__init__(self,paramdict)

#################################################################
# Work out fiducial cosmology
##############################################################
    def assignFiducialCosmology(self,paramdict):
        #Lazy return of default -should update with paramdict values
        return CC.Cosmology() 

#################################################################
# Fisher matrix limits of integration
#################################################################
# These should be defined by the experiment and just wrapped to here

    def maximumK(self):
        return self.Exp.maximumK()

    def minimumK(self):
        return self.Exp.minimumK()

    def minimumThetaK(self,k):
        return self.Exp.thetaminK(k)

    def maximumThetaK(self,k):
        return self.Exp.thetamaxK(k)

#################################################################
# Fisher matrix calculation
#################################################################
    def calculateFisher(self,tagdict,z=0.0):
        """ calculate the Fisher matrix

        SHOULD TAKE PARAMETERS FROM Fisher.paramdict
        

        NEEDS fiducial value for Pk
        """

        #fiducial smology from initial declaration
        cosm=self.cosm

        #fiducial power spectrum
        k,Pk=self.powerFromCAMB(cosm)

        #first check to see if derivatives pre-calculated
        #if not calculate them anew
        #CURRENTLY CALCULATES FROM SCRATCH EACH TIME
        resultdict,tagdict=self.calculatePowerDerivatives(cosm,z,tagdict)

        N=len(tagdict.keys())
        fisherMatrix=np.zeros((N,N))
        for tag1 in tagdict.keys():
            kvec=resultdict[tag1][2]
            i=self._paramdict[tag1]['id']
            for tag2 in tagdict.keys():
                j=self._paramdict[tag2]['id']
                print i, j, tag1, tag2
                dP1=resultdict[tag1][3]
                dP2=resultdict[tag2][3]
                fisherMatrix[i,j]=self.fisherElement(kvec,Pk,dP1,dP2)

        self._fisher=fisherMatrix

        return fisherMatrix


    def fisherElement(self,kvec,Pk,dP1,dP2):
        """
        taking vector of power spectrum derivatives dP1 and dP2
        calculate the amplitude of the Fisher matrix element

        2D integral over first theta and then k
        
        """
        #kstar=self.fisherKStar()
        #thetamin=math.acos(min(self.y*k/(2.0*math.pi),1.0))
        #thetamax=math.asin(min(kstar/k,1.0))

        #remember step size
        self.eta=abs(math.log(kvec[0])-math.log(kvec[1]))

        #obtain limits of integration from class
        kmin=self.minimumK()
        kmax=self.maximumK()

        #carry out integration
        fisher12=0.0
        for (i,k) in enumerate(kvec):
            #print i,k, Pk[i]
            if k>kmin and k<kmax:
                II,error=scipy.integrate.quad(self.kernelFisher,self.minimumThetaK(k),self.maximumThetaK(k),args=(k,Pk[i],dP1[i],dP2[i]))
                fisher12+=II

        return fisher12


    def kernelFisher(self,theta,k,Pk,dP1,dP2):
        """ kernel for the Fisher integral """
        weight=self.Exp.weightVarianceK(theta,k)
        return dP1*dP2/weight**2

    def griddedKernel(self):

        #obtain limits of integration from class
        kperpmin=0.1
        kperpmax=1.0
        kparamin=0.1
        kparamax=1.0
        
        lkperpstep=0.01
        lkparastep=0.01

        lkperp=np.arange(np.log(kperpmin),np.log(kperpmax),lkperpstep)
        lkpara=np.arange(np.log(kparamin),np.log(kparamax),lkparastep)
        xx, yy = np.meshgrid(lkperp,lkpara)
        z=self.dummy(xx,yy)
        f = interpolate.interp2d(x, y, z, kind='linear')
        return xx,yy,f

    def dummy(self,kperp,kpara):
        k=kperp**2+kpara**2
        z=1.0
        mu=kpara/k
        f=self.cosm.fGrowth(z)
        angular=(1.0+f*mu**2)**2
        return self.cosm.linearPowerSpectrum(k,z)*angular
        return 1.0

######## Calculate power spectrum derivatives
    def defaultParamDict(self,cosm):
        paramdict={}
        paramdict['scalar_amp(1)']=cosm.Ascal
        paramdict['ombh2']=cosm.Omegabhh
        paramdict['omch2']=cosm.Omegachh
        paramdict['hubble']=cosm.hubble*100.0
        paramdict['w']= cosm.w
        paramdict['omk']=cosm.Omegak
        paramdict['scalar_spectral_index(1)']=cosm.nscal
        return paramdict
    
    def defaultTagDict(self,cosm):
        vary=0.05
        vary2=0.01
        tagdict={}
        tagdict['scalar_amp(1)']=cosm.Ascal*vary
        tagdict['ombh2']=cosm.Omegabhh*vary
        tagdict['omch2']=cosm.Omegachh*vary
        tagdict['hubble']=cosm.hubble*100.0*vary
        tagdict['w']= cosm.w*vary
        tagdict['omk']=vary2
        tagdict['scalar_spectral_index(1)']=cosm.nscal*vary
        return tagdict

    def calculatePowerDerivatives(self,cosm,z=0.0,tagdict=None):
        if tagdict is None:
            tagdict=self.defaultTagDict()
        resultdict=self.powerDerivative(tagdict,cosm,z)
        return resultdict,tagdict
    

    def powerFromCAMB(self,cosm,tag='fiducial',step=0.0,labeltag='',z=0.0):
        #Establish base values
        mydict=cmbpy.cambio.defaultCosmDict()
        mydict['output_root']='data/test'
        mydict['get_transfer']='T'
        mydict['ombh2']=cosm.Omegabhh
        mydict['omch2']=cosm.Omegachh
        mydict['hubble']=cosm.hubble*100.0
        mydict['w']= cosm.w
        mydict['omk']=cosm.Omegak
        mydict['scalar_amp(1)']=cosm.Ascal
        mydict['scalar_spectral_index(1)']=cosm.nscal
        mydict['transfer_redshift(1)']=z
        mydict['transfer_k_per_logint'] = 5
        mydict['transfer_high_precision'] = 'T'

        #crude hack
        mydict['fiducial']=0
        
        usedict=mydict.copy()
        usedict['output_root']=tag+labeltag
        usedict[tag]=mydict[tag]+step
        cmbpy.cambio.runCAMBFromDict(usedict)
        usefile=usedict['output_root']+'_matterpower.dat'
        data=np.loadtxt(usefile)
        k=data[:,0]
        usecosm=CC.Cosmology(Omegabhh=usedict['ombh2'],Omegak=usedict['omk'],hubble=usedict['hubble']/100.0,nscal=usedict['scalar_spectral_index(1)'],Ascal=usedict['scalar_amp(1)']*1.0e10,Omegam=((usedict['ombh2']+usedict['omch2'])/(usedict['hubble']/100.0)**2),Omegal=1.0-((usedict['ombh2']+usedict['omch2'])/(usedict['hubble']/100.0)**2)-usedict['omk'])
        Pk=self.powerSpectrumWrapper(k,z,usecosm,data[:,1])

        return k,Pk
    
    def powerDerivativeOld(self,tagdict,cosm,z=0.0):
        """ using same tags as CAMB access calculate derivatives
        tagdict should have labels and stepsize tagdict[tag]=step

        My best accuracy on derivatives seems to be fractional error of 0.002
        Worst is 0.02 on omk. This is comparable with the accuracy of CAMB
        which is quoted as being at the level of 0.1-0.2%
        """

        if not os.path.exists("./data"):
            os.mkdir("data")
            self.datadir=os.path.join(os.path.abspath("."),"/data/")

        resultdict={}

        for tag in tagdict.keys():
            step=tagdict[tag]/2.0

            #Establish base values
            mydict=cmbpy.cambio.defaultCosmDict()
            mydict['output_root']='data/test'
            mydict['get_transfer']='T'
            mydict['ombh2']=cosm.Omegabhh
            mydict['omch2']=cosm.Omegachh
            mydict['hubble']=cosm.hubble*100.0
            mydict['w']= cosm.w
            mydict['omk']=cosm.Omegak
            mydict['scalar_amp(1)']=cosm.Ascal
            mydict['scalar_spectral_index(1)']=cosm.nscal
            mydict['transfer_redshift(1)']=z
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
            dataP=np.loadtxt(plusfile)
            dataM=np.loadtxt(minusfile)

            #double plus power spectra
            pplusfile=pplusdict['output_root']+'_matterpower.dat'
            mminusfile=mminusdict['output_root']+'_matterpower.dat'
            dataPP=np.loadtxt(pplusfile)
            dataMM=np.loadtxt(mminusfile)

            #For some reason can't get CAMB to always output same number of
            # steps. Step size is fixed but number varies by 1 or 2
            #
            # NOTE: doesn't yet account for CAMB output being in units of
            # [h/Mpc] and [(Mpc/h)^3]
            minlength=min(len(dataP),len(dataM),len(dataPP),len(dataMM))
            usedict=plusdict
            cosmP=CC.Cosmology(Omegabhh=usedict['ombh2'],Omegak=usedict['omk'],hubble=usedict['hubble']/100.0,nscal=usedict['scalar_spectral_index(1)'],Ascal=usedict['scalar_amp(1)']*1.0e10,Omegam=((usedict['ombh2']+usedict['omch2'])/(usedict['hubble']/100.0)**2),Omegal=1.0-((usedict['ombh2']+usedict['omch2'])/(usedict['hubble']/100.0)**2)-usedict['omk'])
            usedict=minusdict
            cosmM=CC.Cosmology(Omegabhh=usedict['ombh2'],Omegak=usedict['omk'],hubble=usedict['hubble']/100.0,nscal=usedict['scalar_spectral_index(1)'],Ascal=usedict['scalar_amp(1)']*1.0e10,Omegam=((usedict['ombh2']+usedict['omch2'])/(usedict['hubble']/100.0)**2),Omegal=1.0-((usedict['ombh2']+usedict['omch2'])/(usedict['hubble']/100.0)**2)-usedict['omk'])

            #manipulate CAMB output
            k=dataP[0:minlength,0]

            khP=dataP[0:minlength,0]
            khM=dataM[0:minlength,0]
            PkhP=dataP[0:minlength,1]
            PkhM=dataM[0:minlength,1]

            kP=khP/cosmP.hubble
            PkP=PkhP*cosmP.hubble**3
            kM=khM/cosmM.hubble
            PkM=PkhM*cosmM.hubble**3            
            
            PfullP=self.powerSpectrumWrapper(k,z,cosmP,PkP)
            PfullM=self.powerSpectrumWrapper(k,z,cosmM,PkM)
            dP=(PfullP-PfullM)/(2.0*step)

            # double plus
            usedict=pplusdict
            cosmPP=CC.Cosmology(Omegabhh=usedict['ombh2'],Omegak=usedict['omk'],hubble=usedict['hubble']/100.0,nscal=usedict['scalar_spectral_index(1)'],Ascal=usedict['scalar_amp(1)']*1.0e10,Omegam=((usedict['ombh2']+usedict['omch2'])/(usedict['hubble']/100.0)**2),Omegal=1.0-((usedict['ombh2']+usedict['omch2'])/(usedict['hubble']/100.0)**2)-usedict['omk'])
            usedict=mminusdict
            cosmMM=CC.Cosmology(Omegabhh=usedict['ombh2'],Omegak=usedict['omk'],hubble=usedict['hubble']/100.0,nscal=usedict['scalar_spectral_index(1)'],Ascal=usedict['scalar_amp(1)']*1.0e10,Omegam=((usedict['ombh2']+usedict['omch2'])/(usedict['hubble']/100.0)**2),Omegal=1.0-((usedict['ombh2']+usedict['omch2'])/(usedict['hubble']/100.0)**2)-usedict['omk'])

            k=dataP[0:minlength,0]

            khPP=dataPP[0:minlength,0]
            khMM=dataMM[0:minlength,0]
            PkhPP=dataPP[0:minlength,1]
            PkhMM=dataMM[0:minlength,1]

            kPP=khPP/cosmPP.hubble
            PkPP=PkhPP*cosmPP.hubble**3
            kMM=khMM/cosmMM.hubble
            PkMM=PkhMM*cosmMM.hubble**3  

            PfullPP=self.powerSpectrumWrapper(k,z,cosmPP,dataPP[0:minlength,1])
            PfullMM=self.powerSpectrumWrapper(k,z,cosmMM,dataMM[0:minlength,1])

            dPP=(PfullP-PfullM)*8.0/12.0/step-(PfullPP-PfullMM)/12.0/step

            resultdict[tag]=(mydict[tag],step,k,dPP,dP)
            
        return resultdict

    def powerDerivative(self,tagdict,cosm,z=0.0):
        """ using same tags as CAMB access calculate derivatives
        tagdict should have labels and stepsize tagdict[tag]=step

        My best accuracy on derivatives seems to be fractional error of 0.002
        Worst is 0.02 on omk. This is comparable with the accuracy of CAMB
        which is quoted as being at the level of 0.1-0.2%
        """

        if not os.path.exists("./data"):
            os.mkdir("data")
            self.datadir=os.path.join(os.path.abspath("."),"/data/")

        resultdict={}

        #Establish base values
        mydict=cmbpy.cambio.defaultCosmDict()
        mydict['output_root']='data/test'
        mydict['get_transfer']='T'
        mydict['ombh2']=cosm.Omegabhh
        mydict['omch2']=cosm.Omegachh
        mydict['hubble']=cosm.hubble*100.0
        mydict['w']= cosm.w
        mydict['omk']=cosm.Omegak
        mydict['scalar_amp(1)']=cosm.Ascal
        mydict['scalar_spectral_index(1)']=cosm.nscal
        mydict['transfer_redshift(1)']=z
        mydict['transfer_k_per_logint'] = 5
        mydict['transfer_high_precision'] = 'T'

        for tag in tagdict.keys():
            step=tagdict[tag]/2.0
            
            #Run CAMB to get power spectrum
            kP,PkP,PfullP=self.pkFromCAMB(self,mydict,tag,step,label='_p')
            kM,PkM,PfullM=self.pkFromCAMB(self,mydict,tag,-step,label='_m')
            kPP,PkPP,PfullPP=self.pkFromCAMB(self,mydict,tag,2.0*step,label='_pp')
            kMM,PkMM,PfullMM=self.pkFromCAMB(self,mydict,tag,-2.0*step,label='_mm')

            #For some reason can't get CAMB to always output same number of
            # steps. Step size is fixed but number varies by 1 or 2
            #
            # NOTE: doesn't yet account for CAMB output being in units of
            # [h/Mpc] and [(Mpc/h)^3]
            minlength=min(len(PfullP),len(PfullM),len(PfullPP),len(PfullMM))

            #calculate derivatives from power spectrum values
            k=kP
            dP=(PfullP-PfullM)/(2.0*step)
            dPP=(PfullP-PfullM)*8.0/12.0/step-(PfullPP-PfullMM)/12.0/step

            resultdict[tag]=(mydict[tag],step,k,dPP,dP)
            
        return resultdict

    def pkFromCAMB(self,mydict,tag,step,label='_m'):
        """ Run CAMB for the specified cosmology and return both the basic matter
        power spectrum Pk and the brightness temperature one Pfull
        """

        #establish dictionary and call CAMB
        usedict=mydict.copy()
        usedict['output_root']=tag+label
        usedict[tag]=mydict[tag]+step
        cmbpy.cambio.runCAMBFromDict(usedict)

        #read in CAMB file and process
        usefile=usedict['output_root']+'_matterpower.dat'
        data=np.loadtxt(usefile)
        cosm=CC.Cosmology(Omegabhh=usedict['ombh2'],Omegak=usedict['omk'],hubble=usedict['hubble']/100.0,nscal=usedict['scalar_spectral_index(1)'],Ascal=usedict['scalar_amp(1)']*1.0e10,Omegam=((usedict['ombh2']+usedict['omch2'])/(usedict['hubble']/100.0)**2),Omegal=1.0-((usedict['ombh2']+usedict['omch2'])/(usedict['hubble']/100.0)**2)-usedict['omk'])

        #convert CAMB output from units of [Mpc/h] to [Mpc]
        kh=data[0:minlength,0]
        Pkh=data[0:minlength,1]
        k=kh/cosm.hubble
        Pk=Pkh*cosm.hubble**3

        #evaluate the brightness temperature for this matter power spectrum
        Pfull=self.powerSpectrumWrapper(k,z,cosm,Pk)
        return k,Pk,Pfull

    def powerSpectrumWrapper(self,k,z,cosm,pk):
        """ for power spectrum derivatives """
        Tb=self.Exp.brightnessT(z,cosm=cosm)
        return pk*Tb**2
