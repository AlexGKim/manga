#!/usr/bin/env python
import pystan
import numpy
import matplotlib.pyplot as plt
import pickle

with open('fit_lerp.pkl', 'rb') as handle:
    ans=pickle.load(handle)

# plt.plot(ans['scale_a'])
# plt.show()

# plt.plot(ans['sigma'])
# plt.show()



dum = numpy.percentile(numpy.exp(ans['scale_a'])-1,(50,50-34,50+34))
dum=dum*3e5
plt.hist((numpy.exp(ans['scale_a'])-1)*3e5)
plt.xlabel(r'$(\lambda_2/\lambda_1 -1)*3e5$')
plt.savefig('scale.png')
print(r'${:8.0f}_{{-{:8.0f}}}^{{+{:8.0f}}}$'.format(dum[0],dum[0]-dum[1],dum[2]-dum[0]))

# plt.plot(ans['flux'][:,1000])
# plt.show()

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

spec = ans['flux'].mean(axis=0)

plt.plot(filedata[0][:,0],filedata[0][:,1],label='0')
plt.plot(filedata[1][:,0],filedata[1][:,1],label='1')
plt.plot(lam_m,spec,label='fit')
plt.xlim((lmin,lmax))
plt.ylim((0,5))
plt.legend()
plt.savefig('spectra.png')

plt.plot((ans['n_1']))
plt.show()

plt.plot((ans['norm']))
plt.show()

