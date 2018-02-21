#!/usr/bin/env python

from obspy.core import read, UTCDateTime
import glob
import operator
import numpy as np
from matplotlib.mlab import csd
from obspy.signal.invsim import paz_to_freq_resp
from obspy.signal.spectral_estimation import get_nlnm, get_nhnm
import matplotlib.pyplot as plt

import matplotlib as mpl
mpl.rc('font',family='serif')
mpl.rc('font',serif='Times') 
mpl.rc('text', usetex=True)
mpl.rc('font',size=18)

# Compute Cross power
def cp(tr1, tr2, lenfft, lenol, delta):
    sr = 1./delta
    cpval, fre = csd(tr1.data, tr2.data,
                     NFFT=lenfft, Fs=sr, noverlap=lenol,
                     scale_by_freq=True)
    fre = fre[1:]
    cpval = cpval[1:]
    return cpval, fre

# Compute the response
def computeresp(resp, delta, lenfft):
    respval = paz_to_freq_resp(resp['poles'],
                               resp['zeros'],
                               resp['sensitivity']*resp['gain'],
                               t_samp=delta,
                               nfft=lenfft, freq=False)
    respval = np.absolute(respval*np.conjugate(respval))
    respval = respval[1:]
    return respval


# Compute self noise and PSD in dB
def compNoise(idx, length, overlap, delta, instresp):
    cpFix = lambda i1, i2: cp(st[i1], st[i2], length,
                              overlap, delta)

    # We could do these as permutations but that gets confusing
    # Instead we will just hard code the indices
    pp, f = cpFix(idx, idx)
    if idx == 0:
        noisetemp = pp - \
                    cpFix(1, 0)[0]*cpFix(0, 2)[0]/cpFix(1, 2)[0]
    elif idx == 1:
        noisetemp = pp - \
                    cpFix(2, 1)[0]*cpFix(1, 0)[0]/cpFix(2, 0)[0]
    elif idx == 2:
        noisetemp = pp - \
                    cpFix(1, 2)[0]*cpFix(2, 0)[0]/cpFix(1, 0)[0]
    else:
        print 'Bad index, crash landing.'
        sys.exit()
    # Convert to acceleration
    noisetemp *= (2.*np.pi*f)**2
    pp *= (2.*np.pi*f)**2
    # Remove the response
    noisetemp *= 1./instresp
    pp *= 1./instresp
    # Convert to dB
    noisetemp = 10.*np.log10(np.absolute(noisetemp))
    pp = 10.*np.log10(np.absolute(pp))
    return noisetemp, pp, f


debug = True
stime = UTCDateTime('2018-052T01:00:00.0')
etime = stime + 60.*60.

#Here is my estimated response
paz= {'poles': [(149.67189745697971+24.351503981794483j), (149.67189745697971-24.351503981794483j), -8.1574647449752771, -517.98353468592916, 
    (-18.688964026369025-96.597388009454505j)], 'sensitivity': 495550.0, 
    'zeros': [(-15324.27744543273+9331.9159286808681j), (-15324.27744543273-9331.9159286808681j), 
    -653.22466280016897, 3.5520283683062983, 0.13662296274942703], 'gain': 1.0}

files = glob.glob('/tr1/telemetry_days/AM*/' + str(stime.year) + '/*' + str(stime.julday).zfill(3) + '/*EH*')


st = reduce(operator.add, map(read,files))
# You need to be careful as the clock gitter will make the trim function flaky
st.trim(starttime=stime, endtime=etime)
if debug:
    print(st)

# Now we have our data why not compute some spectra?
length =2**12
overlap=2**8
delta = st[0].stats.delta
instresp = computeresp(paz, delta, length)

n1,p1,f = compNoise(0, length, overlap, delta, instresp)
n2,p2,f = compNoise(1, length, overlap, delta, instresp)
n3,p3,f = compNoise(2, length, overlap, delta, instresp)

NLNMper,NLNMpower = get_nlnm()
NHNMper,NHNMpower = get_nhnm()


fig =plt.figure(1,figsize=(12,12))
plt.semilogx(f,p1,label=st[0].id+ ' Power')
plt.semilogx(f,p2,label=st[1].id+ ' Power')
plt.semilogx(f,p3, label=st[2].id + ' Power')
plt.semilogx(f,n1,label=st[0].id + ' Noise')
plt.semilogx(f,n2,label=st[1].id + ' Noise')
plt.semilogx(f,n3, label=st[2].id + ' Noise')
plt.semilogx(1./NLNMper, NLNMpower, linewidth=2., color='k')
plt.semilogx(1./NHNMper, NHNMpower, linewidth=2., color='k',label='NLNM/NHNM')
plt.xlim((200.,.1))
plt.ylim((-200., -100.))
plt.ylabel('Power (dB rel. 1 $(m/s^2)^2/Hz)$')
plt.xlabel('Frequency (Hz)')
plt.title('Self-Noise Estimates for ' + str(stime.year) + ' ' + str(stime.julday).zfill(3) + str(stime.hour).zfill(2) + ':' + str(stime.minute).zfill(2) + ' Duration 1 Hour')
plt.legend()
plt.show()



