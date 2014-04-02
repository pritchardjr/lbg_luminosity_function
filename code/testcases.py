import numpy as np
import cosmologypy.Cosmology as CC
import fisherpy.luminosityfunction as CLF
import matplotlib.pyplot as plt
import scipy.optimize

##############################################################
# Test the halo matching
##############################################################

def halotest():
    cosm=CC.Cosmology()
    LF=CLF.LuminosityFunction()
    z=2
    ms=np.linspace(cosm.coolMass(z),1e15)
    ns=[]
    for mass in ms:
        ns.append(cosm.nCollObject(z,mass,massfcn='PS'))

    plt.figure(1)
    plt.loglog(ms,ns)

    plt.figure(2)
    ngals=[]
    mags=np.linspace(-25,-5)
    for M in mags:
        ngals.append(LF.numberdensityM(M))
        
    plt.yscale('log')
    plt.xscale('linear')
    plt.plot(mags,ngals)


    plt.show()

def halomatch():
    cosm=CC.Cosmology()
    LF=CLF.LuminosityFunction()
    z=2

    M= -15
    ngal=LF.numberdensityM(M)
    mmin=cosm.coolMass(z)/100.0
    mmax=1.0e20

    #throws error if the upper limit is high enough that dndlM is zero
    lim=cosm.dndlM(z,mmax,massfcn='PS')
    while(lim==0):
        mmax/=10.0
        lim=cosm.dndlM(z,mmax,massfcn='PS')
        
    
    mass=scipy.optimize.brentq(lambda x:cosm.nCollObject(z,x,massfcn='PS')-ngal,mmin,mmax)

    print mmax
    print mass, cosm.nCollObject(z,mass,'PS'),ngal
