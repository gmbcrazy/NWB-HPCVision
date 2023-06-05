# -*- coding: utf-8 -*-
"""
Created on Wed Oct  5 17:06:25 2022

@author: lzhang481
"""
##pip install nwbwidgets
##pip install ghostscript
#pip install ghostipy
#pip install numpy --upgrade
#pip install scipy
#import numpy
#numpy.__version__


from pynwb import NWBHDF5IO
from nwbwidgets import nwb2widget

#import ghostscript

import ghostipy
import matplotlib.pyplot as plt


file_path='E:/000061/sub-Rat08/sub-Rat08_ses-Rat08-20130709_ecephys+image.nwb';
io = NWBHDF5IO(file_path, mode='r')
nwb = io.read()


timeI=list(range(50000,52000))
Chan=5
tempData=nwb.processing['ecephys'].data_interfaces['LFP'].electrical_series['LFP'].data[timeI,Chan];
plt.plot(tempData)


