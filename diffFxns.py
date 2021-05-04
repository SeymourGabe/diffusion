#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May  3 22:02:02 2021

@author: gabeseymour
"""

import numpy as np
import matplotlib.pyplot as plt

def prnt(strng):
    print(strng)
    

def diffusion1D(species, D, T, fCoords, xCoords):
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
    H = np.exp(-np.square(np.abs(fCoords/f_c))) # Define transfer fxn
    
    fSpecies = np.fft.fftshift(np.fft.fft(species)) # FFT of species array
        # YOU MUST USE FFTSHIFT! This centers fSpecies so the f=0 component
        # is in the middle of the 
    Hshift = np.fft.fftshift(H)
    fSpeciesH = np.flip(np.multiply(fSpecies, H)) # Adding the flip fixes 
        # a species array flip when the array isn't symmetric
    # --- Need to normalize so the f = 0 component == 1. Doing shifting because 
        # --- I'm confused otherwise. Norming conserves mass ---
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
    k_t = 1 # NOTE: The time constant is tau = 1/k_t
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

def timeSim1D(species, I, D, T0, T, steps, xCoords, fCoords):
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
    dt = T/steps
    for i in range(1, steps + 1):
        if i*10 % steps == 0:
            print("\tstep = ", i)
        """
        Each step of the for loop first simulates the diffusion, and then 
        simulates the generation and termination times
        """
        dT = dt * i # Total sim time for ith step
        species = diffusion1D(species, D, dT, fCoords, xCoords)
        species = reaction1D(species, I, dt, xCoords)
    
    return species




