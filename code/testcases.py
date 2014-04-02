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

def halomatch(M= -15):
    cosm=CC.Cosmology()
    LF=CLF.LuminosityFunction()
    z=2

    ngal=LF.numberdensityM(M)
    mass=cosm.halomatch(z,ngal,'PS')

    print mass, cosm.nCollObject(z,mass,'PS'),ngal
    return mass
