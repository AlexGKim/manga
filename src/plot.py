#!/usr/bin/env python
import pystan
import numpy
import matplotlib.pyplot as plt
import pickle
from chainconsumer import ChainConsumer

with open('fit_K.pkl', 'rb') as handle:
    ans=pickle.load(handle)

dum = numpy.percentile(numpy.exp(ans['scale_a'])-1,(50,50-34,50+34))
dum=dum*3e5
print(r'${:8.0f}_{{-{:8.0f}}}^{{+{:8.0f}}}$'.format(dum[0],dum[0]-dum[1],dum[2]-dum[0]))


# plt.plot(ans['flux'].mean(axis=0)[10:-10])
# plt.show()



# plt.plot(ans['scale_a'])
# plt.show()

# plt.plot(ans['sigma2'])
# plt.show()

c = ChainConsumer()
# print(ans['scale_a'].shape)
# print(ans['norm'].shape)
# sdfd
c.add_chain([ans['scale_a'],ans['norm'][:,0],ans['norm'][:,1],numpy.sqrt(ans['sigma2'])], parameters=["$a_2-a_1$", "$n_1$", "$n_2$", r"$\sigma$"])
fig = c.plotter.plot()
fig.tight_layout()
fig.savefig('corner.png')


# plt.hist((numpy.exp(ans['scale_a'])-1)*3e5)
# plt.xlabel(r'$(\lambda_2/\lambda_1 -1)*3e5$')
# plt.savefig('scale.png')

# plt.plot(ans['flux'][:,1000])
# plt.show()

# #range of (1+v) that are possible
# arange =1e-3
# amin=10**-arange
# amax=10**arange

# # data range observer frame
# lmin=3800
# lmax=10000


# Read in data

# filenames = ['manga-7991-12701-38-26.txt','manga-7991-12701-38-38.txt']

# filedata = []
# for filename in filenames:
#     filedata.append(numpy.loadtxt('../data/'+filename))

# # truncate the wavelength range for the data
# w0 = numpy.logical_and(filedata[0][:,0]>= lmin,filedata[0][:,0] < lmax)
# w1 = numpy.logical_and(filedata[1][:,0]>= lmin,filedata[1][:,0] < lmax)

# nlam0= w0.sum()
# nlam1 = w1.sum()

# #native 1-pix resolution is 
# native_R = 4342.444731562085
# R = 2000.
# lam_m = numpy.arange(numpy.log(lmin)-arange, numpy.log(lmax)+arange+1/R/2.,1/R)
# lam_m = numpy.exp(lam_m)


# plt.clf()
# plt.plot(filedata[0][:,0],filedata[0][:,1],label='0')
# plt.plot(filedata[1][:,0],filedata[1][:,1],label='1')
# plt.xlim((lmin,lmax))
# plt.legend()
# plt.savefig('spectra.png')

# plt.clf()
# plt.plot(filedata[0][w0,0],res1,label='0 - background',alpha=0.7)
# plt.plot(filedata[1][w1,0],res2,label='1- background',alpha=0.7)
# plt.plot(lam_m,spec,label='fit',alpha=0.7)
# plt.xlim((lmin,lmax))
# # plt.ylim((0,0.5))
# plt.legend()
# plt.savefig('spectra_sub.png')