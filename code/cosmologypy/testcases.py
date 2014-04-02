import numpy as np
import Cosmology
import matplotlib.pyplot as plt

##############################################################
# Test the halo mass function code
##############################################################

def halotest():
    cosm=Cosmology.Cosmology()
    z=0
    ms=np.linspace(1e2,1e15)
    ns=[]
    for mass in ms:
        ns.append(cosm.dndlM(z,mass,massfcn='PS'))

    plt.figure(1)
    plt.loglog(ms,ns)

    plt.figure(2)
    cns=[]
    for mass in ms:
        cns.append(cosm.cumulativeHaloCount(z,mass,massfcn='PS'))

    plt.loglog(cns)
    plt.show()
