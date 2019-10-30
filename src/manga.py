#!/usr/bin/env python
import pystan
import numpy
import matplotlib.pyplot as plt
import pickle

if __name__ == "__main__":
    plot = False
    chains=1

    lmin=3800
    lmax=10000
    lmin=5000
    lmin=7000
    dlam = 5.

    # Model the underlying spectrum with a spectrum with the following wavelengths

    lam_m = numpy.arange(lmin-5*dlam,lmax+5.5*dlam,dlam)

    # Read in data

    filenames = ['manga-7991-12701-38-26.txt','manga-7991-12701-38-38.txt']

    filedata = []
    for filename in filenames:
        filedata.append(numpy.loadtxt('../data/'+filename))

    # truncate the wavelength range for the data
    w0 = numpy.logical_and(filedata[0][:,0]>= lmin,filedata[0][:,0] < lmax)
    w1 = numpy.logical_and(filedata[1][:,0]>= lmin,filedata[1][:,0] < lmax)

    nlam0= w0.sum()
    nlam1 = w1.sum()

    # the K matrix usses square exponontial kernal with rho = 1 Angstrom
    # This is constant for spectrum 0
    K=numpy.subtract.outer(filedata[0][w0,0],lam_m)
    K = numpy.exp(-K**2/2/dlam**2)

    data = {'lam0': filedata[0][w0,0], 'flux0': filedata[0][w0,1], 'lam1': filedata[1][w1,0], 'flux1': filedata[1][w1,1], 'K_0': K,
        'lam_m':lam_m, 'nlam_m': len(lam_m),'dlam':dlam, 'nlam0':nlam0,'nlam1':nlam1}

    ones = numpy.interp(lam_m,filedata[0][:,0],filedata[0][:,1])
    ones[ones <=0] = 1e-6
    init = [{
        'flux': ones, 'a_1': 1, 'n_1':1,'sigma':0.05
    } for i in range(chains)]
    if plot:
        for i in ['0','1']:
            plt.plot(data['lam'+i],data['flux'+i])
        plt.show()

    sm = pystan.StanModel(file='manga_lerp.stan')

    fit = sm.sampling(data=data, iter=300,init=init,n_jobs=1,chains=chains)
    pickle.dump( fit, open( "fit_lerp.pkl", "wb" ) )
    print(fit)
