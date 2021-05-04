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
#%% ---------- Define constants and coordinates ----------
twoDflag = 0
X = 1 # In ums
Nx = 10000 # Number of x steps. Somewhat arbitrary. Choosing so spacing is about 1 µm
Y = 10 # In ums
Ny = 10000 # Number of y steps. Somewhat arbitrary. Choosing so spacing is about 1 µm
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
dy = Y/Ny # size of each step in ums
yCoords = np.arange(-Y/2, Y/2, dy)
yCoords = np.arange(0, Y/dy+dy)
yCoords = yCoords - yCoords[-1]/2
yCoords = yCoords*Y
# Center of xCoords is xCoords[int(len(xCoords)/2-.5)]
fxCoords = np.arange(0, X/dx+dx) # 
fxCoords = fxCoords/fxCoords[-1] # normalize
fxCoords = (fxCoords - fxCoords[-1]/2)*len(fxCoords) # Shift everything and scale
fyCoords = np.arange(0, Y/dy+dy) # 
fyCoords = fyCoords/fyCoords[-1] # normalize
fyCoords = (fyCoords - fyCoords[-1]/2)*len(fyCoords) # Shift everything and scale


# ----- Define the species/intensity details -----
species1D = np.zeros_like(xCoords)
if twoDflag:
    species2D = np.zeros((len(xCoords), len(yCoords)))


# --- For an impulse centered at f = 0 ---
speciesImpulse1D = species1D.copy()
speciesImpulse1D[int(len(species1D)/2-.5)] = 1

if twoDflag:
    speciesImpulse2D = species2D.copy()
    xlen, ylen = speciesImpulse2D.shape
    speciesImpulse2D[int(xlen/2-.5)][int(ylen/2-.5)] = 1



# --- For a window centered in the vial ---
window1D = 2*.1
speciesWindow1D = species1D.copy() + 1
windowMask1D = np.abs(xCoords) - window1D/2 < 0
speciesWindow1D = speciesWindow1D*windowMask1D
#speciesWindow1D = np.roll(speciesWindow1D, -1)
species1D = speciesWindow1D

if twoDflag:
    fxCoords2D, fyCoords2D = np.meshgrid(fxCoords, fyCoords)
    fSquared = (fxCoords2D**2 + fyCoords2D**2)
    window2D = 2*.1
    speciesWindow2D = species2D.copy() + 1
    windowMask2D = np.abs(fSquared) - window2D/2 < 0
    speciesWindow2D = np.multiply(speciesWindow2D, windowMask2D)
    species2D = speciesWindow2D

# --- Two windows centered around middle ---

species1D = np.zeros_like(xCoords)
window1D = 2*.05
shift1D = .175
speciesWindow1D = species1D.copy() + 1
twoWindowMask01D = np.abs(xCoords) - (window1D/2 + shift1D) < 0
twoWindowMask11D = np.abs(xCoords) > shift1D
twoWindowMask1D = np.logical_and(twoWindowMask01D, twoWindowMask11D)
speciesWindow1D = speciesWindow1D*twoWindowMask1D
#speciesWindow1D = np.roll(speciesWindow1D, -1)
species1D = speciesWindow1D

plt.plot(xCoords, species1D)
#%% ---------- Test Surface Plots ----------
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
I = species1D.copy()

plt.figure(21)
#species1D = np.roll(species1D, -500)
plt.plot(xCoords, species1D/max(species1D), label = 'Normed species at t = 0')
plt.title("Normed species at t = 0")
Tstep = T/(pltSteps)

for i in range(1, pltSteps + 1):
    print("i = ", i)
    species1D = dF.timeSim1D(species1D, I, D, Tstep*(i-1), Tstep*i, tSteps, xCoords, fxCoords)
    #plt.figure()
    #plt.plot(xCoords, np.abs(species1D), label = 'Species at t = ' + str(i*Tstep) + "\nmax = " + str(max(species1D)))
    plt.plot(xCoords, np.abs(species1D), label = 'Species at t = ' + str(i*Tstep))
    #plt.plot(xCoords, (np.abs(species1D)- min(species1D))/(max(species1D) - min(species1D)), label = 'Rescaled species with range ' + str(max(species1D) - min(species1D)))
    plt.title("Species distribution for " + str(pltSteps) + " iterations of " + str(tSteps) + " steps", y = 1.08)
    x1, x2, y1, y2 = plt.axis()
    plt.legend(loc=2, prop={'size': 6})
    #plt.axis((x1, x2, 0, max(species1D) + .01))
    #plt.legend()
    
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

#%% ---------- Testing diffusion
species1D = np.zeros_like(xCoords)

# --- Two windows centered around middle ---

species1D = np.zeros_like(xCoords)
window1D = 2*.05
shift1D = .175
speciesWindow1D = species1D.copy() + 1
twoWindowMask01D = np.abs(xCoords) - (window1D/2 + shift1D) < 0
twoWindowMask11D = np.abs(xCoords) > shift1D
twoWindowMask1D = np.logical_and(twoWindowMask01D, twoWindowMask11D)
speciesWindow1D = speciesWindow1D*twoWindowMask1D
species1D = speciesWindow1D
plt.figure()
plt.plot(xCoords, species1D, label = 'Species at t = 0')
T = .01
Nstep1 = 1
Nstep2 = 100

dT1 = T/Nstep1
dT2 = T/Nstep2
for i in range(1, Nstep1+1):
    species1D, _ = dF.diffusion1DPDE(species1D, D, dT1*i, xCoords)
    if i:
        print("i = ", i)
        plt.plot(xCoords, species1D, label='i = ' + str(i))
        
plt.title("Nstep1 = " + str(Nstep1))
x1, x2, y1, y2 = plt.axis()
plt.axis((.1, .3, 0, 1.1))
plt.legend()


species1D = np.zeros_like(xCoords)
window1D = 2*.05
shift1D = .175
speciesWindow1D = species1D.copy() + 1
twoWindowMask01D = np.abs(xCoords) - (window1D/2 + shift1D) < 0
twoWindowMask11D = np.abs(xCoords) > shift1D
twoWindowMask1D = np.logical_and(twoWindowMask01D, twoWindowMask11D)
speciesWindow1D = speciesWindow1D*twoWindowMask1D
species1D = speciesWindow1D
plt.figure()
plt.plot(xCoords, species1D, label = 'Species at t = 0')

for i in range(1, Nstep2+1):
    species1D, _ = dF.diffusion1DPDE(species1D, D, dT2*i, xCoords)
    if i % 10 == 0:
        print("i = ", i)
        plt.plot(xCoords, species1D, label='i = ' + str(i))
        
plt.title("Nstep2 = " + str(Nstep2))
x1, x2, y1, y2 = plt.axis()
plt.axis((.1, .3, .4, 1.1))
plt.legend()
