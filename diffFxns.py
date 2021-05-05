#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May  3 22:02:02 2021

@author: gabeseymour
"""

import numpy as np
import matplotlib.pyplot as plt

#%% diffusion1DTransferFxn(species, D, T, fCoords, xCoords)

def diffusion1DTransferFxn(species, D, T, fCoords, xCoords):
    #print("diffusion1D")
    # ----- Testing control panel -----
    pltH = 0
    pltFspecies = 0
    pltSpecies = 0
    """
    This function takes the Fourier transform of the species array, multiplies
    it element-wise by the transfer function of Fick's second law, 
    H = exp{ -( |f| / f_c )^2 } where the characteristic frequency is 
    f_c = 1 / sqrt(D * T) where D is the diffusivity of the material and T
    is the time from the start of diffusion till present 

    Parameters
    ----------
    species : Array of floats
        A 1D array describing the distribution of species
    T : Float
        The time from the start of diffusion till present
    fCoords : Array of floats
        A 1D array where each element's value corresponds to the corresponding
            frequency coordinate of that element. If len(fCoords) is odd, then 
            fCoords[int(len(fCoords) - 0.5)] == 0
    xCoords : Array of floats
        A 1D array where each element's value corresponds to the corresponding
            real-space coordinate of that element. If len(xCoords) is odd, then 
            xCoords[int(len(xCoords) - 0.5)] == 0

    Returns
    -------
    species : Array of floats
        A 1D array describing the distribution of species

    """
    
    f_c = 1/np.sqrt(D*T) # Define characteristic frequency
    #H = np.exp(-np.square(np.abs(fCoords/(f_c*max(fCoords))))) # Define transfer fxn
    H = np.exp(-np.square(np.abs(fCoords/(f_c)))) # Define transfer fxn
    
    fSpecies = np.fft.fftshift(np.fft.fft(species)) # FFT of species array
    #fSpecies = fSpecies/fSpecies[int(len(fSpecies)/2 - 0.5)]
        # YOU MUST USE FFTSHIFT! This centers fSpecies so the f=0 component
        # is in the middle of the 
    Hshift = np.fft.fftshift(H)
    fSpeciesH = np.flip(np.multiply(fSpecies, H)) # Adding the flip fixes 
        # a species array flip when the array isn't symmetric
    # --- Need to normalize so the f = 0 component == 1. Doing shifting because 
        # --- I'm confused otherwise. Norming conserves mass ---
    if fCoords[int(len(fSpeciesH)/2 - 0.5)] != 0:
        print("AHHHHHHH SOMETHING IS WRONG")
        print("\t fCoords[int(len(fSpeciesH)/2 - 0.5)] = ", fCoords[int(len(fSpeciesH)/2 - 0.5)])
    
    fSpeciesH = fSpeciesH/fSpeciesH[int(len(fSpeciesH)/2 - 0.5)]

    # --- 
    

    speciesDiff = np.abs(np.fft.fft(fSpeciesH))
    
    if pltH:
        plt.figure()
        plt.plot(fCoords, H, label = 'H')
        plt.plot(fCoords, Hshift, label = 'H shifted')
        plt.title("H for T = " + str(T))
        #plt.axis((-50, 50, 0, 1))
        plt.legend()
    
    if pltFspecies:
        plt.figure()
        plt.plot(fCoords, fSpecies, label = 'fSpecies')
        plt.plot(fCoords, fSpeciesH, label = 'fSpecies*H')
        plt.title("fSpecies for T = " + str(T))
        plt.legend()
        
    if pltSpecies:
        plt.figure()
        plt.plot(xCoords, species, label = 'pre diffusion')
        plt.plot(xCoords, speciesDiff/max(speciesDiff), label = 'post diffusion (normalized)\nMax = ' + str(max(speciesDiff)))
        plt.title("Species for T = " + str(T))
        plt.legend()
        
    return speciesDiff

#%% diffusion1DPDE(species, D, T, xCoords)

def diffusion1DPDE(species, D, T, xCoords):
    """
    Well, I'm not having luck implementing the transfer function diffusion fxn. 
    So, now I'm trying the PDE solution like we discussed in class. Fick's second law
    is du/dt = D * d^2*u/dx^2 + F where u is the species concentration, D
    is the diffusivity, and F is a forcing function. The forcing function
    takes into account termination and generation of species. IF WE CAN DECOUPLE
    the diffusion and reaction steps (which I believe we can for small enough dt), 
    then du/dt = D * d^2*u/dx^2 is the diffusion term. This has the same form 
    as the heat equation. This has time step soln
    u(x, t + dt) = u(x, t) + eta * [u(x + dx, t) + u(x - dx, t) - 2 * u(x,t)]
    where eta = D*dt/(dx^2).
    
    NOTE: I think the for loop method (set runSuperSlow = 1 to see this) is more
    accurate than the roll method IF the boundary values of the species array are 
    non-zero. However, in VAM, we are projecting into a cylindrical vial. 
    Since we don't have perfect index-matching/ray deviation correction to
    correct for the cylindrical vial, we don't print in the outer regions of 
    the material. So, I feel justified in using the roll method! I also ran 
    an initial test where I used ∆x = 10^-4, ∆t = 10^-10 for an eta/D value of
    10^-2, and 20,000 steps. The roll method took less than 4 seconds while 
    the for loop method took about 2 mins and 20 seconds. This ~35x speedup 
    qualitatively yielded the same results for both methods! I also decreased
    ∆t by an order of magnitude and increased Nsteps by an order of magnitude
    and qualitatively got the same results again! 
    eta = 10^-1

    Parameters
    ----------
    species : Array of floats
        A 1D array describing the distribution of species
    D : float
        Diffusivity of the material in µm^/s (for 3D case)
    T : float
        Total diffusion timestep
    xCoords : Array of floats
        A 1D array where each element's value corresponds to the corresponding
            real-space coordinate of that element. If len(xCoords) is odd, then 
            xCoords[int(len(xCoords) - 0.5)] == 0

    Returns
    -------
    species : Array of floats
        A 1D array describing the distribution of species

    """
    runSuperSlow = 0
    dx = np.abs(xCoords[0] - xCoords[1])
    eta = 0.1
    dt = eta * (dx**2) / D
    Nsteps = int(T/dt)
    #eta = D*dt/(dx**2)
    #print("Nsteps = ", Nsteps)
    if runSuperSlow:
        speciesArr = species.copy()
        valsAdd = np.zeros_like(speciesArr)

        for j in range(1, len(speciesArr) - 1):
            valsAdd[j] = speciesArr[j + 1] + speciesArr[j - 1] - 2*speciesArr[j]
            #  if xCoords[j] <= .1:
            #print("valsAdd[x = ", xCoords[j], "] = ", valsAdd[j])
        
        speciesArr = speciesArr + eta*valsAdd
        species = speciesArr/max(speciesArr)
    else:
        for i in range(1, Nsteps+1):
            # ========== Add in some padding on the end ==========
            padded = np.insert(species, 0, 0)
            padded = np.insert(padded, len(padded), 0)
            #padded = p
            paddedMin1 = np.roll(padded, -1)
            paddedPlus1 = np.roll(padded, 1)
            specs = padded + eta*(paddedPlus1 + paddedMin1 - 2*padded)
            species = specs[1:len(specs)-1]
   # print("eta = ", eta)
    
    return species, dt

""" Test nonsense
# X = 1
# Nx = 10000
# D = 10
# T = 1E-6
# dx = X/Nx
# xCoords = np.arange(0, X/dx+dx)
# xCoords = xCoords/xCoords[-1]
# xCoords = xCoords - xCoords[-1]/2
# xCoords = xCoords*X
# species1D = np.zeros_like(xCoords)
# window1D = 2*.1
# speciesWindow1D = species1D.copy() + 1
# windowMask1D = np.abs(xCoords) - window1D/2 < 0
# speciesWindow1D = speciesWindow1D*windowMask1D
# #speciesWindow1D = np.roll(speciesWindow1D, -1)
# species1D = speciesWindow1D


# plt.plot(xCoords, species1D, label = 'initial species')

# species1D, deltaT = diffusion1DPDE(species1D, D, T, xCoords)
# plt.plot(xCoords, species1D, label = 'Diffused species')

# # #species1D = diffusion1DPDE(species1D, 10, np.abs(xCoords[0] - xCoords[1])/10000, xCoords)
# # steps = 20000
# # deltat = np.abs(xCoords[0] - xCoords[1])/100000
# # T = steps*deltat
# # printStep = steps/5

# # for i in range(1, steps):
# #     species1D = diffusion1DPDE(species1D, 10, deltat, xCoords)
# #     #if i%100 == 0:
# #         #print("i = ", i)
# #     if i % printStep == 0:
# #      #   print("------ \ti = ", i, "\t-----")
# #         plt.plot(xCoords, species1D, label = 'species ' + str(i))
        
# plt.title("T = " + str(T) + "\n∆t = " + str(deltaT), y = 1.1)
# plt.axis((0, .2, 0, 1.1))
# plt.legend() 
"""
#%% reaction1D(species, I, dt, xCoords):

def reaction1D(species, I, dt, xCoords):
    pltSpeciesFlg = 0
    """
    This function finds the species distribution after one time step with 
    unimolecular termination occuring.
    This is from a solution to a pde which takes into account the 
    generation and reaction rates of materials, [r'(t)] = k_I * I - k_t * [r(t)].
    To find these rates, you have to run a few experiments so for this project I'm just 
    worrying about the shape of the distribution after a time step for now.
    
    Sol'n: [r(t)] = k_I/k_t * I + exp{-k_t * t}*C. Solving for C gives
    C = [r(t)] - k_I/k_t * I. Therefore, 
    [r(t + dt)] = k_I/k_t * I * (1 - exp{-k_t * dt}) + [r(t)] * exp{-k_t * dt}

    Parameters
    ----------
    species : Array of floats
        A 1D array describing the distribution of species
    I : Array of floats
        A 1D array describing the intensity distribution in the material
    dt : float
        Time step size
    xCoords : Array of floats
        A 1D array where each element's value corresponds to the corresponding
            real-space coordinate of that element. If len(xCoords) is odd, then 
            xCoords[int(len(xCoords) - 0.5)] == 0

    Returns
    -------
    species : Array of floats
        A 1D array describing the distribution of species

    """
    # ---------- Define constants ----------
    # Arbitrarily chosen constants. I know that as t -> inf, r(inf) -> k_I/k_T*I
    # So I'm setting k_I/k_t < 1 to see a noticable reaction change
    k_I = .1
    k_t = 100 # NOTE: The time constant is tau = 1/k_t
    # ----------
    
    speciesRxn = k_I/k_t * I * (1 - np.exp(-k_t * dt)) + species * np.exp(-k_t * dt)
    
    if pltSpeciesFlg:
        plt.figure()
        plt.plot(xCoords, species, label = 'species pre rxn')
        plt.plot(xCoords, speciesRxn, label = 'species post rxn')
        plt.title("Unimolecular termination and generation")
        plt.legend()
        x1, x2, y1, y2 = plt.axis()
        plt.axis((x1, x2, 0, 1))
    
    return speciesRxn

#%% timeSim1D(species, I, D, T0, T, steps, xCoords, fCoords)
def timeSim1D(species, I, D, T0, T, steps, xCoords, fCoords):
    print("in timeStim1D")
    """
    This function simulates the diffusion, termination, and generation of 
    radicals in a material of diffusivity, D, over a total time period, T, for 
    a specified amount of steps given an initial species distribution and a 
    constant intensity distribution.

    Parameters
    ----------
    species : Array of floats
        A 1D array describing the distribution of species
    I : Array of floats
        A 1D array describing the intensity distribution in the material
    D : float
        Diffusivity of the material in µm^/s (for 3D case)
    T0 : float
        Initial runtime of the simulation
    T : float
        The total runtime of the simulation
    steps : Int
        The number of steps the simulation will run though. dt = T/steps
    xCoords : Array of floats
        A 1D array where each element's value corresponds to the corresponding
            real-space coordinate of that element. If len(xCoords) is odd, then 
            xCoords[int(len(xCoords) - 0.5)] == 0
    fCoords : Array of floats
        A 1D array where each element's value corresponds to the corresponding
            frequency coordinate of that element. If len(fCoords) is odd, then 
            fCoords[int(len(fCoords) - 0.5)] == 0

    Returns
    -------
    species : Array of floats
        A 1D array describing the distribution of species

    """
    dt = (T - T0)/steps
    for i in range(1, steps + 1):
        #if i*1000 % steps == 0:
         #   print("\tstep = ", i)
            
        """
        Each step of the for loop first simulates the diffusion, and then 
        simulates the generation and termination times
        """
        dT = dt * i + T0# Total sim time for ith step
        # species = diffusion1DTransferFxn(species, D, dT, fCoords, xCoords)
        species, _ = diffusion1DPDE(species, D, dt, xCoords) 
        #diffusion1DPDE(species, D, T, xCoords)
        species = reaction1D(species, I, dt, xCoords)
    
    return species




