#!/usr/bin/env python
import pystan
import numpy
import matplotlib.pyplot as plt
import pickle

if __name__ == "__main__":

    chains=4
    iteration=300
    n_jobs=1

    #range of (1+v) that are possible
    arange =1e-3
    amin=10**-arange
    amax=10**arange

    # data range observer frame
    lmin=3800
    lmax=10000


    # # Model the underlying spectrum with a spectrum with the following wavelengths
    # # Probably would be better for the model spectrum to be natively in log(lambda)
    # dlam = 2.
    # lam_m = numpy.arange(lmin*amin-dlam,lmax*amax+dlam*1.5,dlam)

    #native 1-pix resolution is 
    native_R = 4342.444731562085
    R = 2000.
    lam_m = numpy.arange(numpy.log(lmin)-arange, numpy.log(lmax)+arange+1/R/2.,1/R)
    lam_m = numpy.exp(lam_m)


    # plot = False
    # chains=1

    # lmin=3800
    # lmax=10000
    # lmin=5000
    # lmin=7000
    # dlam = 5.

    # # Model the underlying spectrum with a spectrum with the following wavelengths

    # lam_m = numpy.arange(lmin-5*dlam,lmax+5.5*dlam,dlam)

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

    order =3
    p = numpy.polyfit(numpy.log(filedata[0][w0,0]),filedata[0][w0,1], order)
    res1 = filedata[0][:,1] - numpy.polyval(p,numpy.log(filedata[0][:,0]))
    w = numpy.logical_or(w0, numpy.abs(res1) <0.5)
    p = numpy.polyfit(numpy.log(filedata[0][w,0]),filedata[0][w,1], order) 
    res1 = filedata[0][w0,1] - numpy.polyval(p,numpy.log(filedata[0][w0,0]))

    p = numpy.polyfit(numpy.log(filedata[1][w1,0]),filedata[1][w1,1], order)
    res2 = filedata[1][:,1] - numpy.polyval(p,numpy.log(filedata[1][:,0]))
    w = numpy.logical_or(w0, numpy.abs(res2) <0.5)
    p = numpy.polyfit(numpy.log(filedata[1][w,0]),filedata[1][w,1], order) 
    res2 = filedata[1][w1,1] - numpy.polyval(p,numpy.log(filedata[1][w1,0]))


    # the K matrix usses square exponontial kernal with rho = 1 Angstrom
    # This is constant for spectrum 0
    K=numpy.subtract.outer(numpy.log(filedata[0][w0,0]),lam_m)
    K = numpy.exp(-K**2/2/(1/R)**2)

    data = {'lam0': numpy.log(filedata[0][w0,0],) 'flux0': res1, 'lam1': numpy.log(filedata[1][w1,0]), 'flux1': res2, 'K_0': K,
        'lam_m':lam_m, 'nlam_m': len(lam_m),'R2':R**2, 'nlam0':nlam0,'nlam1':nlam1}

    ones = numpy.interp(numpy.log(lam_m),numpy.log(filedata[0][:,0],res1))+ numpy.interp(numpy.log(lam_m),numpy.log(filedata[1][:,0],res2))
    ones=ones/2
    # ones[ones <=0] = 1e-6
    init = [{
        'flux': ones, 'a_1': 1, 'n_1':1,'sigma':0.03
    } for i in range(chains)]
    if plot:
        for i in ['0','1']:
            plt.plot(data['lam'+i],data['flux'+i])
        plt.show()

    sm = pystan.StanModel(file='manga.stan')

    fit = sm.sampling(data=data, iter=300,init=iOOnit,n_jobs=1,chains=chains)
    pickle.dump( fit, open( "fit.pkl", "wb" ) )
    print(fit)
