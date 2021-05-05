#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May  3 21:27:57 2021

@author: gabeseymour
"""
import numpy as np
import matplotlib.pyplot as plt
import sys, os 
dir_path = os.path.dirname(os.path.realpath(__file__))
sys.path.append(dir_path)
import diffFxns as dF
import warnings
#from mpl_toolkits.mplot3d import Axes3D
#from matplotlib import cm
#from matplotlib.ticker import LinearLocator

warnings.filterwarnings('ignore')


"""
This file is just to test the functions in diffFxns.py
These are all "code - to - code" tests, but qualitatively it is easy to see the 
way that diffusion progresses.
"""
#%% ---------- Define constants and coordinates ----------
X = 1 # In ums
Nx = 10000 # Number of x steps. Somewhat arbitrary.
D = .01 # In µm/s (for 1D)
T = .01 # Diffusion time in seconds
pltSteps = 15 # Number of plots 
tSteps = 10000 # Number of time steps between plots. Total timesteps are
    # pltSteps*tSteps. Total time is T*pltSteps

# ----- Define coordinates in real and fourier space -----
dx = X/Nx # size of each step in ums
#xCoords = np.arange(-X/2, X/2, dx)
xCoords = np.arange(0, X/dx+dx)
xCoords = xCoords/xCoords[-1]
xCoords = xCoords - xCoords[-1]/2
xCoords = xCoords*X
"""
 ~~~~~ This section doesn't matter right now becuase I couldn't get the 
 ~~~~~ transfer function diffusion approach to behave properly. However,
 ~~~~~ I'm not deleting it because I will be expanding on this code in the future
 ~~~~~ and I don't want to rewrite anything
 
# Center of xCoords is xCoords[int(len(xCoords)/2-.5)]
fxCoords = np.arange(0, X/dx+dx) # 
fxCoords = fxCoords/fxCoords[-1] # normalize
fxCoords = (fxCoords - fxCoords[-1]/2)*len(fxCoords) # Shift everything and scale
fyCoords = np.arange(0, Y/dy+dy) # 
fyCoords = fyCoords/fyCoords[-1] # normalize
fyCoords = (fyCoords - fyCoords[-1]/2)*len(fyCoords) # Shift everything and scale
"""

#%% Create species distributions for testing
# ----- Define the species/intensity details -----
species1D = np.zeros_like(xCoords)



# --- For an impulse centered at f = 0 ---
speciesImpulse1D = species1D.copy()
speciesImpulse1D[int(len(species1D)/2-.5)] = 1
speciesImpulseShifted1D = np.roll(speciesImpulse1D, -1000)


# --- For a window centered in the vial ---
# Just creating a mask of width window1D for a top hat function 
window1D = 2*.1
speciesWindow1D = species1D.copy() + 1
windowMask1D = np.abs(xCoords) - window1D/2 < 0
speciesWindow1D = speciesWindow1D*windowMask1D
speciesWindowShifted1D = np.roll(speciesWindow1D, -1000) # shifting to make sure
    # my code doesn't have any symmetry dependencies

# --- Two windows centered around middle ---

window1D = 2*.05 # Window width
shift1D = .175 # Amount to shift by
speciesDoubleWindow1D = species1D.copy() + 1
# Just made two masks cause it was easier to debug than one
twoWindowMask01D = np.abs(xCoords) - (window1D/2 + shift1D) < 0
twoWindowMask11D = np.abs(xCoords) > shift1D
twoWindowMask1D = np.logical_and(twoWindowMask01D, twoWindowMask11D)
speciesTwoWindow1D = speciesDoubleWindow1D*twoWindowMask1D
speciesTwoWindowShifted1D = np.roll(speciesTwoWindow1D, -500)
#species1D = speciesWindow1D

plt.figure(1)
fig, axs = plt.subplots(1, 2, sharex='col', sharey='row')
axs[0].plot(xCoords, speciesImpulse1D, label = "Impulse on center")
axs[0].set_title('impulse on center')
axs[1].plot(xCoords, speciesImpulseShifted1D, label = "Shifted impulse")
axs[1].set_title('shifted impulse')

plt.figure(2)
fig, axs = plt.subplots(1, 2, sharex='col', sharey='row')
axs[0].plot(xCoords, speciesWindow1D, label = "Window on center")
axs[0].set_title('Window on center')
axs[1].plot(xCoords, speciesWindowShifted1D, label = "Shifted window")
axs[1].set_title('Shifted window')

plt.figure(3)
fig, axs = plt.subplots(1, 2, sharex='col', sharey='row')
axs[0].plot(xCoords, speciesTwoWindow1D, label = "Double window on center")
axs[0].set_title('Double window on center')
axs[1].plot(xCoords, speciesTwoWindowShifted1D, label = "Shifted window")
axs[1].set_title('Shifted double window')
#plt.plot(xCoords, speciesWindow1D, label="Top hat function on center")
#plt.plot(xCoords, speciesWindowShifted1D, label="Shifted top hat function")

#%% ---------- Test Surface Plots ----------
"""
Started using the following when moving to 2D, but with the transfer function
approach not working, 2D has to wait
"""
# specs2D = np.zeros_like(speciesImpulse2D)


# fSquared = (fxCoords2D**2 + fyCoords2D**2)
# specs2D = np.exp(-fSquared/150**2)
# specs2D = specs2D.reshape(fxCoords2D.shape)
# species2D = species2D.reshape(fxCoords2D.shape)
# plt.figure()
# #fig, ax = plt.subplots(subplot_kw={"projection": "3d"})
# ax = plt.axes(projection='3d')
# #surf = ax.plot_surface(xCoords, yCoords, specs2D, cmap=cm.coolwarm, linewidth=0, antialiased=False)
# ax.plot_surface(fxCoordsTest, fyCoordsTest, species2D, cmap='binary')

# ax.zaxis.set_major_locator(LinearLocator(10))
# ax.zaxis.set_major_formatter('{x:.02f}')
# # Add a color bar which maps values to colors.
# #fig.colorbar(surf, shrink=0.5, aspect=5)
# plt.show()


#%% ---------- Test fourier transforms ----------
"""
Used the following to verify FFTs and IFFs were working properly. Not as 
important now that the transfer function approach is failing
"""
fourierTestFlag = 0
if fourierTestFlag:
    Y = 1 # In ums
    Ny = 10000 # Number of y steps. Somewhat arbitrary. 
    dy = Y/Ny
    
    fxCoords = np.arange(0, X/dx+dx) # 
    fxCoords = fxCoords/fxCoords[-1] # normalize
    fxCoords = (fxCoords - fxCoords[-1]/2)*len(fxCoords) # Shift everything and scale
    fyCoords = np.arange(0, Y/dy+dy) # 
    fyCoords = fyCoords/fyCoords[-1] # normalize
    fyCoords = (fyCoords - fyCoords[-1]/2)*len(fyCoords) # Shift everything and scale
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
    
#%% ---------- Test diffusion ----------

T = .01 # Total sim time

# ----- First, impulse diffusion
plt.figure(4)
speciesImpulse1D_diff = np.zeros_like(speciesImpulse1D)
fig, axs = plt.subplots(1, 2, sharex='col', sharey='row')
axs[0].plot(xCoords, speciesImpulse1D, label = "Impulse on center")
axs[0].set_title('impulse on center')
speciesImpulse1D_diff, _ = dF.diffusion1DPDE(speciesImpulse1D, D, T, xCoords)
axs[0].plot(xCoords, speciesImpulse1D_diff/max(speciesImpulse1D_diff), label='T = ' + str(T))
x1, x2, y1, y2 = axs[0].axis()
axs[0].axis((x1, x2, 0, y2+0.1))
axs[0].legend()
axs[0].legend(loc=2, prop={'size': 6})

speciesImpulseShifted1D_diff = np.zeros_like(speciesImpulse1D)
speciesImpulseShifted1D_diff, _ = dF.diffusion1DPDE(speciesImpulseShifted1D, D, T, xCoords)
axs[1].plot(xCoords, speciesImpulseShifted1D/max(speciesImpulseShifted1D), label = "Shifted impulse")
axs[1].set_title('shifted impulse')
axs[1].plot(xCoords, speciesImpulseShifted1D_diff/max(speciesImpulseShifted1D_diff), label = "Shifted impulse")
x1, x2, y1, y2 = axs[1].axis()
axs[1].axis((x1, x2, 0, y2+0.1))
axs[1].legend()
axs[1].legend(loc=2, prop={'size': 6})


plt.figure(5)
speciesWindow1D_diff = np.zeros_like(speciesWindow1D)
speciesWindowShifted1D_diff = np.zeros_like(speciesWindowShifted1D)
speciesWindow1D_diff, _ = dF.diffusion1DPDE(speciesWindow1D, D, T, xCoords)
speciesWindowShifted1D_diff, _ = dF.diffusion1DPDE(speciesWindowShifted1D, D, T, xCoords)

fig, axs = plt.subplots(1, 2, sharex='col', sharey='row')
axs[0].plot(xCoords, speciesWindow1D, label = "Window on center")

axs[0].plot(xCoords, speciesWindow1D_diff, label = "Dual window on center")

axs[0].set_title('Window on center')
axs[0].legend()
axs[0].legend(loc=2, prop={'size': 6})

axs[1].plot(xCoords, speciesWindowShifted1D, label = "Shifted window")
axs[1].plot(xCoords, speciesWindowShifted1D_diff, label = "Shifted window")
axs[1].set_title('Shifted window')
axs[1].legend()
axs[1].legend(loc=2, prop={'size': 6})

plt.figure(6)
speciesTwoWindow1D_diff = np.zeros_like(speciesTwoWindow1D)
speciesTwoWindowShifted1D_diff = np.zeros_like(speciesTwoWindowShifted1D)
speciesTwoWindow1D_diff, _ = dF.diffusion1DPDE(speciesTwoWindow1D, D, T, xCoords)
speciesTwoWindowShifted1D_diff, _ = dF.diffusion1DPDE(speciesTwoWindowShifted1D, D, T, xCoords)
fig, axs = plt.subplots(1, 2, sharex='col', sharey='row')
axs[0].plot(xCoords, speciesTwoWindow1D, label = "Double window on center")
axs[0].plot(xCoords, speciesTwoWindow1D_diff, label = "Double window on center")
axs[0].set_title('Double window on center')
axs[0].legend()
axs[0].legend(loc=2, prop={'size': 6})
axs[1].plot(xCoords, speciesTwoWindowShifted1D, label = "Shifted window")
axs[1].plot(xCoords, speciesTwoWindowShifted1D_diff, label = "Shifted window")
axs[1].set_title('Shifted double window')
axs[1].legend()
axs[1].legend(loc=2, prop={'size': 6})
#plt.plot(xCoords, speciesWindow1D, label="Top hat function on center")
#plt.plot(xCoords, speciesWindowShifted1D, label="Shifted top hat function")


    
#%% ---------- Testing offcenter distributions
"""
Implemented in a much better way with test diffusion
"""


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
# species = dF.diffusion1DTransferFxn(species, D, T, fxCoords, xCoords)
# plt.figure(11)
# plt.plot(xCoords, species, label = 'posdiff')
# plt.title("species distribution")
# plt.legend()
# species = dF.reaction1D(species, I, .1, xCoords)
# plt.figure(11)
# plt.plot(xCoords, species, label = 'postrxn')
# plt.title("species distribution")
# plt.legend()
# species = dF.diffusion1DTransferFxn(species, D, T*2, fxCoords, xCoords)
# plt.figure(11)
# plt.plot(xCoords, species, label = 'posdiff')
# plt.title("species distribution")
# plt.legend()
# species = dF.reaction1D(species, I, .1, xCoords)
# plt.figure(11)
# plt.plot(xCoords, species, label = 'postrxn')
# plt.title("species distribution")
# plt.legend()
# species = dF.diffusion1DTransferFxn(species, D, T*3, fxCoords, xCoords)
# plt.figure(11)
# plt.plot(xCoords, species, label = 'posdiff')
# plt.title("species distribution")
# plt.legend()
# species = dF.reaction1D(species, I, .1, xCoords)
# plt.figure(11)
# plt.plot(xCoords, species, label = 'postrxn')
# plt.title("species distribution")
# plt.legend()

#%% ---------- Testing rxn
"""
Only testing one thing becuase with no diffusion, as t->inf, the radical 
distribution approaches k_I/k_t * I. These constants are experimentally found
and I set them arbitrarily. I set the ratio to be fairly small. Just going to show
rxn for an arbitrary amount of time
"""
T = .01
plt.figure(7)
speciesTwoWindow1D_rxn = np.zeros_like(speciesTwoWindow1D)
I = speciesTwoWindow1D
speciesTwoWindow1D_rxn = dF.reaction1D(speciesTwoWindow1D, I, T, xCoords)
plt.plot(xCoords, speciesTwoWindow1D, label = 'Double window on center')
plt.plot(xCoords, speciesTwoWindow1D_rxn, label = 'Double window rxn')
plt.title("Double window on center")
plt.label()

#%% ---------- Testing timeSim1d
"""
Have 10 mins to submit. Only going to show for centered double window
"""
T = 0.01
plt.figure(8)
plt.plot(xCoords, speciesTwoWindow1D, label = 'Species at t = 0')
#plt.title("Normed species at t = 0")
pltSteps = 5
Tstep = T/(pltSteps)
speciesTwoWindow1DFull = speciesTwoWindow1D.copy()

for i in range(1, pltSteps + 1):
    #print("i = ", i)
    speciesTwoWindow1DFull = dF.timeSim1D(speciesTwoWindow1DFull, I, D, Tstep*(i-1), Tstep*i, tSteps, xCoords, xCoords)
    #plt.figure()
    #plt.plot(xCoords, np.abs(species1D), label = 'Species at t = ' + str(i*Tstep) + "\nmax = " + str(max(species1D)))
    plt.plot(xCoords, speciesTwoWindow1DFull, label = 'Species at t = ' + str(i*Tstep))
    #plt.plot(xCoords, (np.abs(species1D)- min(species1D))/(max(species1D) - min(species1D)), label = 'Rescaled species with range ' + str(max(species1D) - min(species1D)))
    plt.title("Species distribution for " + str(pltSteps) + " iterations of " + str(tSteps) + " steps", y = 1.08)
    x1, x2, y1, y2 = plt.axis()
    plt.legend(loc=2, prop={'size': 6})
    #plt.axis((x1, x2, 0, max(species1D) + .01))
    #plt.legend()

