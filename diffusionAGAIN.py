#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May  3 21:27:57 2021

@author: gabeseymour
"""

import numpy as np
import matplotlib.pyplot as plt
import diffFxns as dF
import warnings
warnings.filterwarnings('ignore')
#%% ---------- Define constants and coordinates ----------

X = 1 # In ums
Nx = 1000 # Number of x steps. Somewhat arbitrary. Choosing so spacing is about 1 µm
D = .01 # In µm/s (for 1D)
T = .15 # Diffusion time in seconds
pltSteps = 15 # Number of plots 
tSteps = 100000 # Number of time steps between plots. Total timesteps are
    # pltSteps*tSteps. Total time is T*pltSteps

# ----- Define coordinates in real and fourier space -----
dx = X/Nx # size of each step in ums
xCoords = np.arange(-X/2, X/2+dx, dx)

# Center of xCoords is xCoords[int(len(xCoords)/2-.5)]
fxCoords = np.arange(0, X/dx+dx) # 
fxCoords = fxCoords/fxCoords[-1] # normalize
fxCoords = (fxCoords - fxCoords[-1]/2)*len(fxCoords) # Shift everything and scale

# ----- Define the species/intensity details -----
species = np.zeros_like(xCoords)

# --- For an impulse centered at f = 0 ---
speciesImpulse = species.copy()
speciesImpulse[int(len(species)/2-.5)] = 1

# --- For a window centered in the vial ---
window = 2*.1
speciesWindow = species.copy() + 1
windowMask = np.abs(xCoords) - window/2 < 0
speciesWindow = speciesWindow*windowMask
speciesWindow = np.roll(speciesWindow, -1)
species = speciesWindow

# --- Two windows centered around middle ---
species = np.zeros_like(xCoords)
window = 2*.05
shift = .175
speciesWindow = species.copy() + 1
twoWindowMask0 = np.abs(xCoords) - (window/2 + shift) < 0
twoWindowMask1 = np.abs(xCoords) > shift
twoWindowMask = np.logical_and(twoWindowMask0, twoWindowMask1)
speciesWindow = speciesWindow*twoWindowMask
#speciesWindow = np.roll(speciesWindow, -1)
species = speciesWindow

#%% ---------- Test fourier transforms ----------

fourierTestFlag = 0
if fourierTestFlag:
    # ----- Test cosine -> impulse transforms
    f1 = 20
    f2 = 100

    test1 = np.zeros_like(xCoords)
    test1 = np.cos(2*np.pi*f1*xCoords)
    test2 = np.zeros_like(xCoords)
    test2 = np.cos(2*np.pi*f2*xCoords)
    
    plt.figure()
    plt.plot(xCoords, test1, label = 'f = ' + str(f1))
    plt.plot(xCoords, test2, label = 'f = ' + str(f2))
    plt.title("Cosine transform tests - real space")
    plt.axis((-.1, .1, -1.1, 1.1))
    plt.legend()
    
    ftest1 = np.zeros_like(test1)
    ftest1 = np.fft.fftshift(np.fft.fft(test1))
    ftest2 = np.zeros_like(test2)
    ftest2 = np.fft.fftshift(np.fft.fft(test2))
    
    plt.figure()
    plt.plot(fxCoords, ftest1, label = 'f = ' + str(f1))
    plt.plot(fxCoords, ftest2, label = 'f = ' + str(f2))
    plt.title("Cosine transform tests - Fourier space")
    x1, x2, y1, y2 = plt.axis()
    plt.axis((0, 120, y1, y2))
    plt.legend()
    
    # ----- Test window -> sinc transform
    # Ft of a window of width w becomes a sinc where first null is at 2*pi/w
    test3 = np.zeros_like(xCoords)
    test4 = np.zeros_like(xCoords)
    w3 = 2*.25 # in µm. Width of the window
    w4 = 2*.45 # in µm. Width of the window
    test3 = test3 + 1
    test4 = test4 + 1
    windowMask3 = np.abs(xCoords)- w3/2 < 0
    windowMask4 = np.abs(xCoords)- w4/2 < 0
    test3 = test3*windowMask3
    test4 = test4*windowMask4
    
    plt.figure()
    plt.plot(xCoords, test3, label = 'w = ' + str(w3))
    plt.plot(xCoords, test4, label = 'w = ' + str(w4))
    plt.title("Window transform tests - real space")
    x1, x2, y1, y2 = plt.axis()
    plt.legend()
    
    ftest3 = np.zeros_like(test3)
    ftest3 = np.fft.fftshift(np.fft.fft(test3))
    ftest4 = np.zeros_like(test4)
    ftest4 = np.fft.fftshift(np.fft.fft(test4))
    
    plt.figure()
    plt.plot(fxCoords, ftest3, label = 'w = ' + str(w3))
    plt.plot(fxCoords, ftest4, label = 'w = ' + str(w4))
    plt.title("Window transform tests - Fourier space")
    x1, x2, y1, y2 = plt.axis()
    plt.axis((-10, 10, y1, y2))
    plt.legend()
    
#%% ---------- Test diffFxns.py ---------
dF.prnt("test")
I = species.copy()

plt.figure(21)
plt.plot(xCoords, species/max(species), label = 'Normed species at t = 0')
plt.title("Normed species at t = 0")
Tstep = T/(pltSteps)
for i in range(1, pltSteps + 1):
    print("i = ", i)
    species = dF.timeSim1D(species, I, D, Tstep*(i-1), Tstep*i, tSteps, xCoords, fxCoords)
    plt.figure()
    plt.plot(xCoords, species, label = 'Normed species at t = ' + str(i*Tstep))
    plt.title("Normalized species distribution for " + str(pltSteps) + " iterations of " + str(tSteps) + " steps")
    plt.legend()
    
#%% ---------- Testing offcenter distributions



# # --- for a window centered in the vial ---
# species = np.zeros_like(xCoords)
# window = 2*.2
# specieswindow = species.copy() + 1
# windowmask = np.abs(xCoords) - window/2 < 0
# specieswindow = specieswindow*windowmask
# specieswindow = np.roll(specieswindow, -100)
# species = specieswindow
# i = species.copy()
# plt.figure(11)
# plt.plot(xCoords, species, label = 'prediff')
# species = dF.diffusion1D(species, D, T, fxCoords, xCoords)
# plt.figure(11)
# plt.plot(xCoords, species, label = 'posdiff')
# plt.title("species distribution")
# plt.legend()
# species = dF.reaction1D(species, I, .1, xCoords)
# plt.figure(11)
# plt.plot(xCoords, species, label = 'postrxn')
# plt.title("species distribution")
# plt.legend()
# species = dF.diffusion1D(species, D, T*2, fxCoords, xCoords)
# plt.figure(11)
# plt.plot(xCoords, species, label = 'posdiff')
# plt.title("species distribution")
# plt.legend()
# species = dF.reaction1D(species, I, .1, xCoords)
# plt.figure(11)
# plt.plot(xCoords, species, label = 'postrxn')
# plt.title("species distribution")
# plt.legend()
# species = dF.diffusion1D(species, D, T*3, fxCoords, xCoords)
# plt.figure(11)
# plt.plot(xCoords, species, label = 'posdiff')
# plt.title("species distribution")
# plt.legend()
# species = dF.reaction1D(species, I, .1, xCoords)
# plt.figure(11)
# plt.plot(xCoords, species, label = 'postrxn')
# plt.title("species distribution")
# plt.legend()