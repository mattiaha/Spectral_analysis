# -*- coding: utf-8 -*-
"""
Created on Tue May  9 11:03:06 2023

@author: matti
"""
from gammapy.astro.darkmatter import (DarkMatterAnnihilationSpectralModel)
import matplotlib.pyplot as plt
import numpy as np
"""
scale = np.load("s_vt23.npy")
mass = np.load("m_vt3.npy")
err_= np.load("errors3.npy")"""
therm = 3e-26
import plotly.graph_objs as go
scale = np.load("s_vt.npy")

#exclude last 11 entries as they were undefined for y
y = np.load("s_vt23.npy")[:-11]*therm
x = np.load("m_vt3.npy")[:-11]
err= np.load("errorst3.npy")[:-11]*therm

y_upper = y+err
y_lower = y-err




plt.figure(0)

plt.yscale("log")
plt.ylabel("<$\sigma v$ > $(cm^3s^{-1}$")
plt.xscale("log")
plt.xlabel("Mass [TeV]")

plt.fill_between(x, y_lower, y_upper)
plt.hlines(y=therm,xmin=min(x),xmax = max(x), color='r', linestyle='-')



plt.show()

