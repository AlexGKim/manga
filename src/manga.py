#!/usr/bin/env python
import pystan
import numpy
import pickle

'''
Given a set of central wavelengths of a Manga spectrum, return a mask of elements that have telluric features.

That this could be made much more efficient but since it isn't run often I don't bother.
'''

def telluric_mask(lam):
    data = numpy.loadtxt('../data/telluric.txt')
    logicals=[]
    # I don't know exactly what the edges are of the MANGA spectral bins so I make a guess
    edges = (lam+numpy.roll(lam,1))/2
    edges[0]=lam[0]-(lam[1]-lam[0])/2
    edges = numpy.append(edges, lam[-1]+(lam[-1]-lam[-2])/2)

    mask = []
    # loop over data bins
    for i in range(len(lam)):
        temp=[]
        # Loop over telluric bins
        for lmin,lmax in zip(data[:,0],data[:,1]):
            # Do the spectral and telluric bins overlap
            temp.append((edges[i] <= lmin  and lmin<=edges[i+1]) or (lmin<=edges[i+1] and edges[i+1]<=lmax))
        mask.append(numpy.any(temp))

    return numpy.logical_not(mask)

if __name__ == "__main__":

    # 
    chains=4
    iteration=2000
    warmup=1000
    n_jobs=4

    # range of ln(1+v) that are possible
    # This helps STAN out a bit
    arange =1e-3
    amin=10**-arange
    amax=10**arange

    # data range observer frame
    # The data will be cut by this
    # The model range will be this expanded by the possible velocity range
    lmin=3800
    lmax=10000


    # Model the underlying spectrum with a spectrum with the following wavelengths
    # Note that the native resolution of MANGA is constant R, so it makes sense to
    # work with a constant resolution grid

    #native 1-pix resolution is 
    native_R = 4342.444731562085
    R = 2000    # This is somewhat arbitrary, could be increase a little
    lam_m = numpy.arange(numpy.log(lmin)-arange-5/R, numpy.log(lmax)+arange+5.5/R,1/R)
    lam_m = numpy.exp(lam_m)

    # Read in data
    filenames = ['manga-7991-12701-38-26.txt','manga-7991-12701-38-38.txt']
    filedata = []
    for filename in filenames:
        filedata.append(numpy.loadtxt('../data/'+filename))

    # truncate the wavelength range for the data
    w0 = numpy.logical_and.reduce((filedata[0][:,0]>= lmin,filedata[0][:,0] < lmax,filedata[0][:,3]==1,telluric_mask(filedata[0][:,0])))
    w1 = numpy.logical_and.reduce((filedata[1][:,0]>= lmin,filedata[1][:,0] < lmax,filedata[1][:,3]==1,telluric_mask(filedata[1][:,0])))

    nlam0= w0.sum()
    nlam1 = w1.sum()

    # subtract out a third order polynomial to get rid of smooth continuum
    # I did try low-pass filter with FFT but these leave periodic features
    # we want to avoid

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

    # Do the STAN HMC
    # Note that the wavelengths are in log space.
    data = {'lam0': numpy.log(filedata[0][w0,0]), 'flux0': res1, 'flux0_var':1/filedata[0][w0,2], 'lam1': numpy.log(filedata[1][w1,0]), 'flux1': res2, 'flux1_var':1/filedata[1][w1,2],
        'lam_m':numpy.log(lam_m), 'nlam_m': len(lam_m),'nlam0':nlam0,'nlam1':nlam1,'arange':arange, 'R': R}

    # As an initial condition take the average spectrum
    ones = numpy.interp(numpy.log(lam_m),numpy.log(filedata[0][w0,0]),res1) + numpy.interp(numpy.log(lam_m),numpy.log(filedata[1][w1,0]),res2) 
    ones = ones/2.

    init = [{
        'flux': ones, 'a_scale':0.0, 'n':[0.5,0.5], 'sigma2':0.03**2
    } for i in range(chains)]

    sm = pystan.StanModel(file='manga.stan')
    fit = sm.sampling(data=data,init=init,chains=chains,n_jobs=n_jobs,iter=iteration,warmup = warmup)

    # Save results
    with open('fit.pkl', 'wb') as handle:
        pickle.dump(fit.extract(), handle)
    print(fit)