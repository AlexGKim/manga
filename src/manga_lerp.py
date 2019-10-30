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


    # Model the underlying spectrum with a spectrum with the following wavelengths
    # Probably would be better for the model spectrum to be natively in log(lambda)
    dlam = 5
    lam_m = numpy.arange(lmin*amin-dlam,lmax*amax+dlam*1.5,dlam)

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

    # Model spectrum in log(lambda), STAN runs faster with sum rather than multiplication


    data = {'lam0': numpy.log(filedata[0][w0,0]), 'flux0': filedata[0][w0,1], 'lam1': numpy.log(filedata[1][w1,0]), 'flux1': filedata[1][w1,1],
        'lam_m':numpy.log(lam_m), 'nlam_m': len(lam_m),'nlam0':nlam0,'nlam1':nlam1,'arange':arange}

    # As an initial condition take the average spectrum
    ones = numpy.interp(numpy.log(lam_m),numpy.log(filedata[0][:,0]),filedata[0][:,1]) + numpy.interp(numpy.log(lam_m),numpy.log(filedata[1][:,0]),filedata[1][:,1]) 
    ones = ones/2.
    ones[ones <=0] = 1e-2

    init = [{
        'flux': ones, 'a_scale':0.0, 'n':[0.5,0.5], 'sigma':0.03
    } for i in range(chains)]

    # init = [{
    #     'flux': ones, 'n':[0.5,0.5], 'sigma':0.03
    # } for i in range(chains)]

    sm = pystan.StanModel(file='manga_lerp.stan')

    fit = sm.sampling(data=data,init=init,chains=chains,n_jobs=n_jobs,iter=iteration)

    with open('fit_lerp.pkl', 'wb') as handle:
        pickle.dump(fit.extract(), handle)
    print(fit)

    # plt.plot(data['lam0'],data['flux0'],label='0')
    # plt.plot(data['lam1'],data['flux1'],label='1')
    # plt.plot(lam_m,ones,label='avg')
    # plt.legend()
    # plt.show()
    # wef
