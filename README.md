# Diffusion and Reaction of a Species in a Diffuse Material
This repo consists of a python module developed in Phys 5070 at CU Boulder during the Spring 2021 semester. This module simulates the distribution of a species given arbitrary, spatially invariant diffusivities and arbitrary, spatially varying initial species concentrations for a specified time.

## Background
I'm in my first year of a Ph.D. in electrical engineering in [Professor Robert McLeod's research group](https://www.colorado.edu/faculty/mcleod/) at CU Boulder where I work on Volumetric Additive Manufacturing (VAM). VAM is a new 3D printing technique in which images are sequentially projected in a rotating vial of photocurable resin. In one voxel of resin, if the accumulated optical dose in that voxel crosses a threshold known as the gelation threshold, t<sub>g</sub>, then the resin in that voxel will gel, effectively turning from a liquid to a solid. 
### Resin
The most basic photocurable resin consists of a writing monomer and a photoreactive radical generator. When the photoreactive radical generator creates a radical, it is then free to either react with a writing monomer, or it can terminate. If it reacts with a writing monomer, the writing monomer then has a radical on it and is free to react with other monomers or polymer chains. There radical-monomer groups then react and create linear growth polymers. These polymers can then crosslink and gel. 
### Fick's second law
Fick's second law gives a relationship between the change in concentration of a species with respect to time and the change in the spatial concentration of the same species (dø/dt = D*d(dø/dx)/dx where D is the diffusivity of our material and ø is the concentration of our species). It's important to note that I am specifically interested in modeling the concentration of radicals in my resin as this concentration directly maps to where I do and do not have gelation/printing. I can add a forcing term to Fick's - a term that describes the spatial distribution of the light dose applied to the resin which, with some scaling factor, maps to a radical generation term, and another term that accounts for radical termination. Then, with a small enough time step, it should be a fine enough approximation to numerically solve the PDE describing JUST diffusion (not termination), and then numerically solve the PDE describing termination. In other words, I can figure out the way the radicals diffuse in my material for some small time step, and then I can figure out which radicals are terminated.

### Assumptions:
- I'm assuming only unimolecular termination is occuring 
- Diffusivity is spatially invariant
- Diffusion and termination are separable. This assumption is a bit iffier, but as there are no VAM (volumetric Additive Manufacturing) groups in the world modeling diffusion, a module in our code doing so, even with this assumption, will be a big step up in the world of VAM. With a small enough time scale, this approximation shouldn't be an issue. 
## Theory
A simple way to conceptualize the printing process is that in one local spot of the resin, if there are enough polymer radical-monomer reactions, there will be gelation (i.e. printing) in that local spot. How we typical talk about this is if our dose in a local spot crosses a GELATION THRESHOLD (caps because this is important), there is gelation in that spot. This language implies an applied dose->radical generation from photoreactive species reacting to applied dose->radical-monomer rxn ->radical on monomer-monomer reaction (forming a polymer)->radical on polymer-monomer reaction (which repeats) process. These polymers then cross-link and solidify. This causes the local density to increase which increases the diffusivity in this region. In other words, applying dose to a local spot causes more stuff to get shoved into that local spot. More stuff means it's harder for species to push through that local spot (i.e. diffuse through that spot). Since most of the "important" diffusion occurs before gelation, it is a fine enough approximation to say that diffusivity is constant spatially and temporally. There are groups that have spent a few PhDs to model the change in diffusivity as cross-link density changes, and no one gives any validity to their findings because there are too many variables unaccounted for and lots of handwaving went into finding the results they found. 

## Motivation
One of our collaborators prints into materials much less viscous than ours (their diffusivity is higher). When they try to print something with sharp edges in this material, they get blobby edges. To fix this, they go into their CAD software and cut away some of the sharp region to "deblob." They then print the new model and see if they're close to the desired geometry. As you can imagine, this is a slow, iterative process. In upcoming work my group is doing, we will need some sort of approach to account for diffusion in order to achieve higher fidelity prints. Adding a module to model the diffusion of radicals is debatably the best way to do this. This can also be added to an algorithm  (recently submitted for publication by someone in my lab) which iteratively changes the model to achieve better prints. If done right, this project may lead to a minor publication (likely after changing the code to fit in our current algorithm after the class is over). 

## Functions
Surprisingly, there are just four functions used in my code. At this point, I have only implemented this in 1D. Unfortunately, my initial approach to finding the diffusion of a time step didn't pan out, though it is significantly less computationally expensive for large arrays than what I've implemented.
### diffusion1DTransferFxn(species, D, T, fCoords, xCoords):
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
