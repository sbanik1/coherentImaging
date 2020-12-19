#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This code sets the directories for running the test codes
Created on Sat Dec 19 12:51:22 2020
@author: Swarnav Banik
sbanik1@umd.edu
"""
# %% Import all ###############################################################
import sys
import matplotlib.pyplot as plt

# %% Add necessary paths ######################################################
sys.path.insert(1, '/Users/swarnav/Google Drive/Work/Projects/Imaging/src')
# %% Define the output directory ##############################################
saveDir = '/Users/swarnav/Google Drive/Work/Projects/Imaging/test/out'
# %% Set some default values ##################################################
plt.rcParams.update({'font.size': 14})