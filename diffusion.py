#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr  8 20:23:23 2021

@author: gabeseymour
"""

import numpy as np
import matplotlib.pyplot as plt

#fun = 1
"""
Adding this for fun. Make fun = 0 to turn off
"""
#if fun:
#    import mplcyberpunk # Kind of a cool thing *shrug*
#    plt.style.use("cyberpunk")

def FFT1(real):
    """
    

    Parameters
    ----------
    real : Real space 1D function

    Returns
    -------
    shifted FFT of real: fourier space representation of the passed function. 
        Shifted as to center the Fourier spectrum

    """
    return np.fft.ifftshift(np.fft.fft(real))

def IFFT1(fourier):
    """
    

    Parameters
    ----------
    fourier : 1D Fourier spectrum centered in the middle of the array

    Returns
    -------
    inverse shifted real space function: inverse shifted and IFFT'd function

    """
    return np.abs(np.fft.ifftshift(np.fft.ifft(fourier)))

def defineCoords(dimensions, steps):
    ### NOT CURRENTLY USED
    """
    Creates real space numpy array and the corresponding Fourier 
    numpy array based on physical (dimensions) X, Y, and Z dimensions, (steps) dx, dy, and dz spacial
    steps. 
    """
    
    print("1D")
    xCoords = np.arange(-X/2, X/2+dx, dx) # 1D 
    fxCoords = np.arange(0, X/dx+dx)
    fxCoords = fxCoords - fxCoords[-1]/2 # Shift everything over so the center of array is at f = 0
    realCoords = xCoords
    fourierCoords = fxCoords
    realSpace = np.zeros(Nx+1)
    xRealSpace = np.zeros(Nx)
    #realSpace = np.zeros(Nx) #1D
    
    fourierSpace = np.zeros_like(realSpace, complex)

    return realCoords, fourierCoords, realSpace, fourierSpace

    
def testCoords1D(realSpace, fourierSpace, X, dx):
    """
    Tests a set of 1D fourier transforms to make sure our coordinates
    are defined correctly
    """
    
    def cosines(freqs, Nyquist, real, fourier, x):
        for i in range(0, len(freqs)):
            real = real + np.cos(2*np.pi * freqs[i] * x)
            
        fourier = np.fft.fftshift(np.fft.fft(real))
        return real, fourier
    
    def gauss(sigma, mu, real, fourier, x): # Still need to fix
        real = 1/(sigma * np.sqrt(2 * np.pi)) * np.exp( - (x - mu)**2 / (2 * sigma**2) )
        fourier = np.fft.fftshift(real)
        return real, fourier
    
    def window(w, real, fourier, x):
        real = real*0 + 1
        mask = np.abs(x)-w/2 < 0
        real = real*mask
        #fourier = np.fft.fftshift(np.fft.fft(real))
        fourier = FFT1(real)
        #f1 = np.fft.fft(real)
        #ifft = np.fft.ifft(np.fft.fftshift(fourier))
        #ifft = np.fft.ifft(np.fft.ifftshift(fourier))
        ifft = IFFT1(fourier)
        return real, fourier, ifft
    
    N = len(realSpace)
    xVals = np.arange(-X/2, X/2 + dx, dx)
    f = np.linspace(0, X, N)/X
    const = 0*realSpace + 1
    impulse = np.fft.fft(const)
    #plt.plot(f, impulse)
    cosines, impulses = cosines([1], X/2, const*0, f, xVals)
    #plt.plot(xVals, cosines)
   # plt.plot(f*X, np.abs(impulses))
    gauss, gaussFT = gauss(1, 0, const*0, f, xVals)
    window, sinc, iwindow = window(X/4, const, f, xVals)
    plt.plot(f - f[-1]/2, np.abs(sinc))
    #plt.plot(xVals, iwindow)
   # plt.plot(xVals, np.abs(gauss))
    #plt.plot(f, np.abs(gaussFT))
    
def findMiddle(inList):
    """
    I'm lazy. This is easy. Taken from 
    https://stackoverflow.com/questions/38130895/find-middle-of-a-list/38131003

    Parameters
    ----------
    inList : List of length N

    Returns
    -------
    middleCoords : indices of middle
    
    middleVals : Values of the middle of the list

    """
    middle = float(len(inList))/2
    if middle % 2 != 0:
        return int(middle - .5), inList[int(middle - .5)]
    else:
        return (int(middle), int(middle-1)), (inList[int(middle)], inList[int(middle-1)])
        
def diffusion1(species, diffusivity, t_D, f_coords, x_coords, plotflg):
    """
    

    Parameters
    ----------
    species : Array of floats
        Each element's value corresponds with the number of that species in 
        that real-space coordinate
    diffusivity : Float
        D of the material. In µm^2/s. For my holographic material, the 
        diffusivity is 1 µm^2/s. This is fairly low. For testing, use 10x this
        which is arbitrarily picked
    t_D : float
        Diffusion time
    f_coords : Array of floats
        The frequency coordinate of each element
        
    Returns
    -------
    species : Array of floats
        The species distribution after diffusion

    """
    # plt.figure(21)
    # plt.plot(x_coords, species)
    # x1, x2, y1, y2 = plt.axis()
    # plt.axis((-.1, .1, 0, y2))
    # plt.title('species pre diffusion')
    
    
    
    #print("len(f_coords) = ", len(f_coords))
    #print("len(species) = ", len(species))
    #print("len(f_coords) = ", len(f_coords))
    f_c = 1/np.sqrt(diffusivity*t_D) # Characteristic spatial frequency
    #print("dT = " + str(t_D) + ", f_c = " + str(f_c))
    #print("f_c = ", f_c)
    H = np.exp(-np.square(np.abs(f_coords*len(f_coords))/f_c)) # Transfer function
    if plotflg:
        plt.figure(95)
        plt.plot(f_coords, np.fft.fftshift(H), label = 't_D = ' + str(t_D))
        plt.title("H plots")
        plt.legend(loc=2, prop={'size': 6})
        #plt.legend()
        
    Hmax = max(range(len(H)), key=H.__getitem__)
    """
    Hmax = max(range(len(H)), key=H.__getitem__)
    print("Hmax index= ", Hmax)
    Hmaxval = H[Hmax]
    print("Hmax value= ", Hmaxval)
    print("f_coords[Hmaxindex] = ", f_coords[Hmax])
    H1 = H[:Hmax + 1]
    H2 = np.flip(H[Hmax:int(len(H))])
    Htest = H1 - H2
    print("max Htest = ", max(Htest))
    """
   #  plt.figure(69)
   #  plt.plot(f_coords, H)
   #  x1, x2, y1, y2 = plt.axis()
   # # plt.axis((-.01, .01, .999, 1.0))
   #  plt.title('H')
    #plt.figure(1)
    #print("f_coords = ", np.fft.fftshift(f_coords))
    #print("H = ", H)
    #plt.plot(x_coords, species)
    f_species = FFT1(species)
    #f_species = np.fft.fft(species)
    
    # plt.figure(70)
    # plt.plot(f_coords, np.abs(f_species), label = 'max = ' + str(max(f_species)) + '\nf_species[Hmax_index] = ' + str(f_species[Hmax]))
    # x1, x2, y1, y2 = plt.axis()
    # #plt.axis((-.1, .1, 0.8, y2))
    # plt.title('f_species before norm')
    # plt.legend()
    
    #f_species = f_species/f_species[int(len(f_species)/2-.5)] # Gotta normalize
        # So species at f = 0 is 1 
    f_species1 = f_species.copy()
    #plt.figure(2)
    f_species = np.multiply(f_species, H)
    f_species = f_species/f_species[int(len(f_species)/2-.5)] # Gotta normalize
    f_species2 = f_species.copy()
    
    # plt.figure(71)
    # plt.plot(f_coords, f_species, label = 'max = ' + str(max(f_species)) + '\nf_species[Hmax_index] = ' + str(f_species[Hmax]))
    # x1, x2, y1, y2 = plt.axis()
    # #plt.axis((-0.5, .5, 0, y2))
    # plt.title('f_species * H')
    # plt.legend()
    
    #fSpeciesMax = max(range(len(f_species)), key=f_species.__getitem__)
    #print("f_species max index= ", fSpeciesMax)
    #f_speciesval = f_species[fSpeciesMax]
    #print("f_species value= ", f_speciesval)
    #print("f_coords[fSpeciesMaxindex] = ", f_coords[fSpeciesMax])
    #fS1 = f_species[:fSpeciesMax + 1]
    #fS2 = np.flip(f_species[fSpeciesMax:int(len(f_species))])
    #fStest = fS1 - fS2
    #print("max fStest = ", max(fStest))
    #print("f_species*H[middle - 1] = ", str(f_species[fSpeciesMax - 1]))
    #print("f_species*H[middle + 1] = ", str(f_species[fSpeciesMax + 1]))
    
    #plt.plot(f_coords, np.abs(f_species))
    species = np.fft.ifftshift(IFFT1(f_species))
    #species = np.fft.ifftshift(np.fft.ifft(np.abs(f_species)))
    #species = np.fft.ifftshift(np.fft.ifft(f_species))
    ### !!! ADDING A BANDAID FOR NOW
    species = np.roll(species, -1)
    #speciesMax = max(range(len(species)), key=species.__getitem__)
    #print("species max index= ", speciesMax)
    #speciesval = species[speciesMax]
    #print("species value= ", speciesval)
    #print("x_coords[speciesMaxindex] = ", x_coords[speciesMax])
    #S1 = species[:speciesMax + 1]
    #S2 = np.flip(species[speciesMax:int(len(species))])
    #Stest = S1 - S2
    #print("max Stest = ", max(Stest))
    #print("IFFT{species*H[middle - 1]} = ", str(species[speciesMax - 1]))
    #print("IFFT{species*H[middle + 1]} = ", str(species[speciesMax + 1]))
    
    
    
    # plt.figure(72)
    # plt.plot(x_coords, np.abs(species), label = 'max = ' + str(max(species)) + '\nspecies[Hmax_index] = ' + str(species[Hmax]))
    # x1, x2, y1, y2 = plt.axis()
    # plt.axis((-.1, .1, 0, y2))
    # plt.title('IFFT{species*H}')
    # plt.legend()
    #print("IFFT{species*H}[middle - 1] = ", str(np.abs(species[Hmax - 1])))
    #print("IFFT{species*H}[middle + 1] = ", str(np.abs(species[Hmax + 1])))
    
    # plt.figure(22)
    # plt.plot(x_coords, species)
    # x1, x2, y1, y2 = plt.axis()
    # plt.axis((-.1, .1, 0, 1.00001))
    # plt.title('species post diffusion')
    """
    if t_D == 10/10:
        plt.figure(2)
        middle, middleVals = findMiddle(f_species)
        
        plt.plot(f_coords, np.abs(f_species1), label = "f_species")
        plt.plot(f_coords, H, label = "H")
        plt.plot(f_coords, np.abs(f_species2), label = "f_species*H")
        x1, x2, y1, y2 = plt.axis()
        plt.axis((-.1, .1, 0, y2))
        plt.legend()
    """    
    """
    plt.figure(2)
    plt.plot(f_coords, H, label = 'dt = ' + str(t_D))   
    plt.axis((-.05, .05, 0, 1))
    plt.legend()
    #print("plotting f vs F{species}")
    #plt.figure()
    #plt.plot(f_coords, np.abs(f_species), label='fourier space')
    
    #plt.figure()
    #print("plotting x vs. species")
    #plt.plot(x_coords, species, label='real space')
    plt.figure(3)
    plt.plot(x_coords, species)
    """
    return species
    
def rxn1(species, I, dTau, xCoords, i):
    """
    This function finds the species distribution after one time step.
    This is from a solution to a pde which takes into account the 
    generation and reaction rates of materials. To find these rates, 
    you have to run a few experiments so for this project I'm just 
    worrying about the shape of the distribution after a time step for now.
    
    r(tau + dTau) = c*exp{-R_m/R_r(tau + dTau)} + I where 
    c = [r(tau) - I]*exp{R_m/R_r*tau}
    Therefore, r(tau + dTau) = [r(tau) - I]*exp{R_m/R_r*tau}*exp{-R_m/R_r(tau + dTau)} + I
    = [r(tau) - I]*exp{-R_m/R_r*tau} + I
    
    Parameters
    ----------
    species : 1D array of floats
        Describes the real-space distribution of the species in question
    intensity : 1D array of floats
        Describes the spacial intensity distribution which is constant over 
            a time step. 
    dTau : float
        One time step

    Returns
    -------
    species : 1D array of floats
        Describes the real-space distribution of the species in question

    """
    # plt.figure(23)
    # plt.plot(xCoords, species)
    # x1, x2, y1, y2 = plt.axis()
    # plt.axis((-.1, .1, 0, y2))
    # plt.title('species pre reaction')
    
    
    if i == 0:
        plt.figure(111+i)
        plt.plot(xCoords, I, label = 'I')
        plt.plot(xCoords, species, label = 'Species')
        
    R_m = 1
    R_r = .1
    #species_dTau = (species - I)*np.exp(-R_m/R_r*dTau) + I
    species_dTau = species*np.exp(-R_m/R_r * dTau) + I*(1 - np.exp(-R_m/R_r*dTau))
    if i == 0:
        plt.plot(xCoords, species_dTau, label = 'Species dTau')
        plt.legend()
        x1, x2, y1, y2 = plt.axis()
        plt.axis((-.1, .1, 0, y2))
        plt.title("Reaction step")
    
    # plt.figure(24)
    # plt.plot(xCoords, species)
    # x1, x2, y1, y2 = plt.axis()
    # plt.axis((-.1, .1, 0, y2))
    # plt.title('species post reaction')
    return species_dTau

#rlCoords, frCoords, rlSpace, frSpace = defineCoords([10, 0, 0], [.1, 0, 0])

# ---------- Contants ----------
X = 1 # In ums
Nx = 1000 # Number of x steps. Somewhat arbitrary. Choosing so spacing is about 1 µm
D = 1 # In µm/s (for 1D)
t_D = .001 # Diffusion time in seconds

# ---------- Define coordinates in real and fourier space ----------
dx = X/Nx # size of each step in ums
xCoords = np.arange(-X/2, X/2+dx, dx)
#print("len(xCoords) = ", str(len(xCoords)))
# Center of xCoords is xCoords[int(len(xCoords)/2-.5)]
fxCoords = np.arange(0, X/dx+dx) # 
fxCoords = fxCoords/fxCoords[-1] # normalize
fxCoords = fxCoords - fxCoords[-1]/2
#print("fxCoords[1] = ", fxCoords[1])
#print("fxCoords[-2] = ", fxCoords[-2])
#print("fxCoords[int(len(fxCoords)/2-.5)] = ", fxCoords[int(len(fxCoords)/2-.5)])
#print("dfx = ", fxCoords[1] - fxCoords[0])
# plt.figure(0)
# plt.plot(fxCoords, fxCoords)
# plt.plot(fxCoords[int(len(fxCoords)/2-.5)], fxCoords[int(len(fxCoords)/2-.5)], '.', color='red')
# plt.plot(fxCoords[int(len(fxCoords)/2-.5) + 1], fxCoords[int(len(fxCoords)/2-.5) + 1], '.', color='green', label = 'fxCoords[int(len(fxCoords)/2-.5) + 1] = ' + str(fxCoords[int(len(fxCoords)/2-.5) + 1]))
# plt.plot(fxCoords[int(len(fxCoords)/2-.5) + 2], fxCoords[int(len(fxCoords)/2-.5) + 2], '.', color='green', label = 'fxCoords[int(len(fxCoords)/2-.5) + 2] = ' + str(fxCoords[int(len(fxCoords)/2-.5) + 2]))

# plt.plot(fxCoords[int(len(fxCoords)/2-.5) - 1], fxCoords[int(len(fxCoords)/2-.5) - 1], '.', color='black', label = 'fxCoords[int(len(fxCoords)/2-.5) - 1] = ' + str(fxCoords[int(len(fxCoords)/2-.5) - 1]))
# plt.axis((-.01, .01, -.01, .01))
# plt.legend()

# plt.title("fxCoords")
# ---------- Create species array ----------

# ----- For impulse -----
"""
species = np.zeros_like(xCoords)
species[int(len(species)/2-.5)] = 1
plt.plot(xCoords, species)
#species = diffusion1(species, D, t_D, fxCoords, xCoords)

for i in range(1, 11):
    idt = t_D/(11-i)
    species = diffusion1(np.abs(species), D, idt, fxCoords, xCoords)
    #species = species/species[int(len(species)/2-.5)] # Gotta normalize
    plt.plot(xCoords, np.abs(species), label = 'tau = ' + str(idt)) # Plot the species profile. Should be an impulse in the middle for now
    plt.legend()

"""
# ----- For window -----
"""
species = np.zeros_like(xCoords)
w = 1 # in µm. Width of the window
species = species + 1
windowMask = np.abs(xCoords)- w/2 < 0
species = species*windowMask
plt.figure(1)
plt.plot(xCoords, species)
for i in range(1, 11):
    
    idt = t_D/(11-i)
    species = diffusion1(np.abs(species), D, idt, fxCoords, xCoords)
    species = species/species[int(len(species)/2-.5)] # Gotta normalize
    plt.figure(1)
    plt.plot(xCoords, np.abs(species*(11-i)*.1), label = 'tau = ' + str(idt)) # Plot the species profile. Should be an impulse in the middle for now
    
plt.legend(loc=2, prop={'size': 7})
"""

# ----- For testing reaction with an intensity impulse distribution
species = np.zeros_like(xCoords)
w = .1 # in µm. Width of the window
species = species + 1
windowMask = np.abs(xCoords)- w/2 < 0
species = species*windowMask
#species[int(len(species)/2-.5)] = 1
intensity = species.copy()
#species = species
plt.figure(1)
plt.plot(xCoords, np.abs(species))
plt.title('Initial species')
#plt.axis((-.1, .1, 0, 1))

# plt.figure(2)
# plt.plot(xCoords, intensity)
# plt.title('Intensity distribution')
# plt.axis((-.1, .1, 0, 1))
steps = 100001
#specsDiff = np.zeros((len(species), steps))
#specsRxn = np.zeros((len(species), steps))
#specsDiff[:, 0] = species.copy()
#specsRxn[:, 0] = species.copy()

pltStep = int(steps/10)
dt = t_D/steps
plotflag = 0

for i in range(1, steps + 1):
    plotflag = 0
    dTau = t_D/steps
    dTau = dt
    #idt = t_D/(steps - i)
    idt = dt*i
    if i % 10000 == 0 or i >= (steps - pltStep) or i == 1:
        print("i = ", i, "\t ", steps-i, " steps to go \t t = ", idt)
    
    num = i + 20
    #print("num = ", num)
    if i % pltStep == 0 or i == 1 or i == steps:
        plotflag = 1
       
        
    species = np.abs(diffusion1(species, D, idt, fxCoords, xCoords, plotflag))
    #specsDiff[:, i] = species.copy()
    if i % pltStep == 0 or i == 1 or i == steps:
        plt.figure()
        plt.plot(xCoords, species, label = 'diffusion ' + str(i))
        plt.legend()
        
    species = np.abs(rxn1(np.abs(species), intensity, dTau, xCoords ,i))

    #specsRxn[:, i] = species.copy()
    #if i % pltStep == 0 or i == 1 or i == steps:
        #plt.plot(xCoords, np.abs(species), label = 'reaction')
        
    if i % pltStep == 0 or i == 1 or i == steps:
        plt.title('Steps = ' + str(i) + ', t = ' + str(idt))
        x1, x2, y1, y2 = plt.axis()
        #plt.axis((-0.1, 0.1, 0, y2))
        plt.legend()

    

#species = diffusion1(species, D, t_D, fxCoords, xCoords)


#species = diffusion1(species, 10000, t_D, fxCoords, xCoords)
# for i in range(1, 11):
#     idt = t_D/(11-i)
#     species = diffusion1(species, D, idt, fxCoords, xCoords)
#     species = species/species[int(len(species)/2-.5)] # Gotta normalize
#     plt.plot(xCoords, np.abs(species)) # Plot the species profile. Should be an impulse in the middle for now



"""
shift = np.fft.fftshift(species)

plt.figure()
plt.plot(xCoords, shift)
plt.title("shift")


Fshift = np.fft.fft(shift)

plt.figure()
Fshift = Fshift
plt.plot(fxCoords, Fshift)
plt.title("Fshift")

plt.figure()
Fspecies = np.fft.fft(species)
plt.plot(fxCoords, Fspecies)
plt.title("Fspecies")

plt.figure()
IFshift =  np.fft.ifftshift(Fspecies)
plt.plot(fxCoords, IFshift)
plt.title("IFshift")

plt.figure()
IFFT_IFshift = np.abs(np.fft.ifft(IFshift))
plt.plot(xCoords, IFFT_IFshift)
plt.title("IFFT of IFshift")

plt.figure()
IFshift_IIFT_IFshift = np.fft.ifftshift(IFFT_IFshift)
plt.plot(xCoords, IFshift_IIFT_IFshift)
plt.title("IFshift of the IFFT of IFshift")

f1 = np.fft.fft(species)
f2 = np.fft.ifftshift(f1)
f3 = np.fft.ifft(f2)
f4 = np.abs(f3)
f5 = np.fft.ifftshift(f4)

plt.figure()
plt.plot(xCoords, species)
plt.title("species")

plt.figure()
plt.plot(fxCoords, f1)
plt.title("fft of species")

plt.figure()
plt.plot(fxCoords, f2)
plt.title("IFFT shift of fft of species")

plt.figure()
plt.plot(fxCoords, f3)
plt.title("ifft of IFFT shift of fft of species")

plt.figure()
plt.plot(xCoords, f4)
plt.title("abs of IFFT shift of fft of species")

plt.figure()
plt.plot(xCoords, f5)
plt.title("reconstructed")
"""

#%%

test = np.zeros_like(xCoords)
test = np.cos(2*np.pi*20*xCoords)
plt.plot(xCoords, test)

fTest = np.fft.fftshift(np.fft.fft(test))
plt.plot(fxCoords*len(fTest), fTest)
x1, x2, y1, y2 = plt.axis()
plt.axis((-50, 50, 0, y2))