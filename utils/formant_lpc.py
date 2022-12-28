# -*- coding: utf-8 -*-

import numpy as np
from .lpc import levinson_1d, lpc_ref
from scipy.signal import lfilter, fftconvolve
from scipy.signal.windows import hamming, kaiser

def intensity(audio,sr, winlen=0.025, winstep=0.01):
    nWin = int(winlen*sr)
    nStep = int(winstep*sr)
    # audio = np.concatenate((audio,np.zeros((nWin,))))
    nPts = len(audio)
    nFrames = int((nPts-nWin)/nStep + 1)
    I = np.zeros((nFrames,))
    
    # G = np.zeros(nFrames,)
    k = 0
    audio = lfilter([1., -.975], 1, audio) 
    for j in range(0,nPts-nWin,nStep):    
        Ia = fftconvolve(audio[j:j+nWin]**2,kaiser(nWin, beta=20),mode='same')
        I[k] = 20*np.log10(np.mean(Ia)/2e-5)
        k = k+1
    return I

def format_lpc(audio,sr, nFormReq=4, maxFreq = 4000, winlen=0.01, winstep=0.01):
    if (0.5*sr > maxFreq):
        nForm = int(0.5*sr/1000)
    if (0.5*sr < maxFreq):
        maxFreq = 0.5*sr
    order = 2*nForm+1
    nWin = int(winlen*sr)
    nStep = int(winstep*sr)
    # audio = np.concatenate((audio,np.zeros((nWin,))))
    nPts = len(audio)
    nFrames = int((nPts-nWin)/nStep + 1)
    F = np.zeros((nFormReq,nFrames))
    B = np.zeros((nFormReq,nFrames))
    # G = np.zeros(nFrames,)
    k = 0
    audio = lfilter([1., -.975], 1, audio) 
    for j in range(0,nPts-nWin,nStep):
        wAudio = audio[j:j+nWin]*hamming(nWin)
        # a, _, _ = levinson_1d(wAudio,order)
        a = lpc_ref(wAudio,order)
        rs = np.roots(a)
        ff = np.angle(rs)*sr*0.5/np.pi
        fs = np.sort(ff)
        fi = np.argsort(ff)
        bb = -np.log(np.abs(rs))*sr/np.pi
        bs = bb[fi]
        idx = np.where((fs>0)*(fs<maxFreq))[0]
        if (len(idx) > nFormReq):
            idx = idx[:nFormReq]
        if (len(idx) < nFormReq):
            fsel = np.concatenate((fs[idx],np.zeros((nFormReq-len(idx),))))
            bsel = np.concatenate((bs[idx],np.zeros((nFormReq-len(idx),))))
        else:
            fsel = fs[idx]
            bsel = bs[idx]
        F[:,k] = fsel
        B[:,k] = bsel
        k = k+1
    # if (nFormReq < nForm):
    #     F = F[:,:nFormReq]
    #     B = B[:,:nFormReq]
    # if (nFormReq > nForm):
    #     F = np.concatenate((F,np.zeros((nFormReq-nForm,nFrames))),axis=0)
    #     B = np.concatenate((B,np.zeros((nFormReq-nForm,nFrames))),axis=0)
    return F, B

