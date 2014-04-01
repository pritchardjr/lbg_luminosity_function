#
# Calculate the basic statistics commonly measured from a box
#
#
import numpy as np

def powerspectrum(kbox,nk,kcell=1):
    """  From a numpy 3D k-space box calculate the power spectrum
    optionally specify the size of a k-cell with kcell

    P(k)=(Sum_{cells of k} delta_k delta_k)/N_k

    NOTE: poorly coded and untested
    """

    kvec=np.zeros(nk)
    pk=np.zeros(nk)
    vark=np.zeros(nk)
    count=np.zeros(nk)
    kmin=0.0
    dim=len(kbox)
    kmax=np.sqrt(3.0)*dim
    kstep=(kmax-kmin)/float(nk)
    kvec=np.arange(nk)*kstep
    for i in range(dim):
        for j in range(dim):
            for k in range(dim):
                kk=np.sqrt(i*i+j*j+k*k)
                indx=int(kk/kstep)
                pk[indx]=pk[indx]+kbox[(i,j,k)]**2
                vark[indx]=vark[indx]+kbox[(i,j,k)]**4
                count[indx]=count[indx]+1

    for i in range(nk):
        pk[i]=pk[i]/float(count[i])

    return (kvec*kcell,pk)


def correlationFunction3D(box,nbin,Lcell=1.0):
    """ calculate 3D correlation function assuming homogeneity

    nbin= number of bins for correlation function
    
    Optionally specify Lcell to scale to physical size
    """

    #first determine bins for storing correlation function
    dim=len(box)
    rbins=np.zeros(nbin)
    xibins=np.zeros(nbin)

    rmin=1.0  #one cell spacing
    rmax=np.sqrt(3.0)*dim #diagonally opposite cells
    rstep=(rmax-rmin)/float(nbin)
    rbins=np.arange(nbin)*rstep
    rbins,rstep=np.linspace(rmin,rmax,nbin,retstep=True)

    #loop over all pairs of points to calculate correlation function
    #this is O(N^2) so expensive. Could do via FFT too as O(NlogN)

    #I'm lame so do this in O(3N^2)
    #first form list of all positions
    positions=[]
    for i in range(dim):
        for j in range(dim):
            for k in range(dim):
                print i,j,k
                positions.append(np.array((i,j,k)))

    #now have list of all positions so look at pairs
    for position1 in positions:
        for position2 in positions:
            r=np.linalg.norm(position2-position1)
            indx=int(r/rstep)
            xibins[indx]=xibins[indx]+1
    


    return rbins,xibins/float(dim**3)

def correlationFunctionByFFT3D(box,nbin,Lcell=1.0):
    """ use FFT """

    #fft to get k-space box
    kbox=np.fftn(box)
    #from k-space box get power spectrum
    pk=1
    #ifft to get correlation function
    rbox=np.ifft(pk)

    return rbox
