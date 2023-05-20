# -*- coding: utf-8 -*-
"""
Created on Tue Feb  7 16:26:11 2023

@author: matti
"""
import os
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt

import os

import astropy
import astropy.units as u
from astropy.io import fits


from astropy.coordinates import SkyCoord

import gammapy
from gammapy.maps import Map, MapAxis, WcsGeom
from gammapy.irf import load_cta_irfs
from gammapy.datasets import MapDataset, MapDatasetEventSampler
from gammapy.makers import MapDatasetMaker, SafeMaskMaker
from gammapy.data import Observation, observatory_locations, DataStore
from gammapy.modeling import Fit
from gammapy.modeling.models import PowerLawNormSpectralModel, Models,  FoVBackgroundModel,GaussianSpatialModel, SkyModel, TemplateSpatialModel
from gammapy.astro.darkmatter import (
    DarkMatterAnnihilationSpectralModel)
from astropy.coordinates import Angle, SkyCoord
from regions import CircleSkyRegion

path =  Path(os.environ["CTADATA"])/"data/baseline/gc"

irf_path = Path(os.environ["CALDB"])/"South_z20_50h"



paths = list(path.rglob("*.fits"))

#Here we can inspect the observations.
data_store = DataStore.from_events_files(paths,irf_path)
table = data_store.obs_table



pos_obs = SkyCoord(table['GLON_PNT'], table['GLAT_PNT'], frame='galactic', unit='deg')
pos_target = SkyCoord(0, 0, frame='galactic', unit='deg')
offset = pos_target.separation(pos_obs).deg

mask = (offset > 0.5) & (offset <1)

table = table[mask]


table.show_in_browser(jsviewer=True)
length = len(table)
nr_o = 50
obs_ids = []
id_int = np.linspace(0,length-1, nr_o,  dtype = int)

for i in range(nr_o):
    a = id_int[i]
    obs_ids.append(table[a][0])
np.save("obs_ids_dc.npy",obs_ids)

