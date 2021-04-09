#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr  8 20:23:23 2021

@author: gabeseymour
"""

import numpy as np
from scipy.special import gamma
import matplotlib.pyplot as plt

def defineCoords(M, N, O):
    """
    Creates an MxNxO real space numpy array and the corresponding Fourier 
    numpy array 
    """
    realSpace = np.zeros([M, N, O])
    return realSpace

print(defineCoords(2, 3, 4))