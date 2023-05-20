# -*- coding: utf-8 -*-
"""
Created on Mon Jan  9 14:26:42 2023

@author: matti
"""

import numpy as np
import matplotlib.pyplot as plt

import os

from pathlib import Path
import astropy
import astropy.units as u
from astropy.io import fits
from astropy.coordinates import SkyCoord

import gammapy
from gammapy.maps import Map, MapAxis, WcsGeom, WcsNDMap
from gammapy.irf import load_cta_irfs
from gammapy.datasets import MapDataset, MapDatasetEventSampler
from gammapy.makers import MapDatasetMaker, SafeMaskMaker
from gammapy.data import Observation, observatory_locations, DataStore
from gammapy.modeling import Fit
from gammapy.modeling.models import PowerLawNormSpectralModel, Models,  FoVBackgroundModel,GaussianSpatialModel, SkyModel, TemplateSpatialModel
from gammapy.astro.darkmatter import (
    DarkMatterAnnihilationSpectralModel,JFactory,
    PrimaryFlux,
    profiles,
)
from astropy.coordinates import Angle, SkyCoord
from regions import CircleSkyRegion
from astropy.time import Time
#Load the models we will use to define the galaxy center
sky_models = Models.read("gc_models.yaml")
#Load the irfs of our choice
irfs = load_cta_irfs("$CALDB/South_z20_50h/irf_file.fits")

y_pos = np.linspace(-1.0,1.0,5)
x_pos = np.array([359.0,359.5,0.0,0.5,1.0])

livetime = 0.5 * u.hr
location = observatory_locations["cta_south"]
#Decide where we want to create our observation



energy_axis = MapAxis.from_energy_bounds("0.1 TeV", "20 TeV", nbin=20, per_decade=True)
energy_axis_true = MapAxis.from_energy_bounds(
    "0.03 TeV", "30 TeV", nbin=30, per_decade=True, name="energy_true"
)
migra_axis = MapAxis.from_bounds(0.5, 2, nbin=150, node_type="edges", name="migra")
bkg_model = FoVBackgroundModel(dataset_name="Obs-1")


profile = profiles.EinastoProfile(r_s=20 * u.kpc , alpha=0.17,
rho_s=0.081* u.Unit("GeV / cm3"))
profiles.DMProfile.DISTANCE_GC = 8.5*u.kpc
profiles.DMProfile.LOCAL_DENSITY = 0.39 * u.Unit("GeV / cm3")

# This changes the normalisation parameter: rho_s
profile.scale_to_local_density()

position = SkyCoord(0.0, 0.0, frame="galactic", unit="deg")
geom = WcsGeom.create(binsz=0.1, skydir=position, width=6,frame="galactic")

jfactory = JFactory(
    geom=geom, profile=profile, distance=profiles.DMProfile.DISTANCE_GC
)
jfact = jfactory.compute_jfactor()
jfact_map = WcsNDMap(geom=geom, data=jfact.value, unit=jfact.unit)
spatial_model = TemplateSpatialModel(jfact_map)
channel = "W"
massDM = 10.0*u.Unit("TeV")
jfactor = jfact.sum()
j_Einasto = 7.1*10**22*u.Unit("GeV2/cm5")

DarkMatterAnnihilationSpectralModel.THERMAL_RELIC_CROSS_SECTION
modelDM = DarkMatterAnnihilationSpectralModel(mass=massDM,
channel=channel, jfactor=j_Einasto)
dm_model = SkyModel(
    spatial_model=spatial_model, spectral_model=modelDM,
    name="model-simu"
)
models = Models([sky_models[0],sky_models[1],sky_models[2],sky_models[3],sky_models[4],
                 sky_models[5],sky_models[6], dm_model,  bkg_model])
models.name = "Models"
models.datasets_names = "models"
obs_id = 1000
maker = MapDatasetMaker(selection=["exposure", "background", "psf", "edisp"])

maker_safe_mask = SafeMaskMaker(
    methods=["aeff-max"], aeff_percent = 10
)
numb = 0

for a in range(5):
    x = x_pos[a]
    for b in range(5):
        y = y_pos[b]
        pointing = SkyCoord(x, y, frame="galactic", unit="deg")
        #We want to avoid simulations of the galactic center for our spectral analysis
        if x==0 and y == 0:
            continue
        #We want to avoid centering on the huge gamma ray source around this position as it will be masked
        elif x==359 and y == -1.0:
            continue
        else:
            
            for j in range(10):
                observation = Observation.create(
                    obs_id=obs_id,
                    pointing=pointing,
                    livetime=livetime,
                    irfs=irfs,
                    location=location,
                    deadtime_fraction=0.02,
                    reference_time =Time("2020-01-01 00:00:00"))
                obs_id += 1
                
                geom = WcsGeom.create(skydir = pointing,             #The coordinates to which the central position of the skymap corresponds
                                binsz=0.02, 
                                width=(4,4),                #Width (in degrees) of the skymap
                                frame ="galactic", 
                                axes=[energy_axis])
                
                empty = MapDataset.create(
                    geom,
                    energy_axis_true=energy_axis_true,
                    migra_axis=migra_axis,
                    name ='Obs-1',
                    reference_time='2020-01-01 00:00:00'
                )
                
                dataset = maker.run(empty, observation)
                dataset = maker_safe_mask.run(dataset,observation)
                
                dataset.models = models
            
                sampler = MapDatasetEventSampler(random_state=j)
                events = sampler.run(dataset, observation)
                
                filename = "C:/Users/matti/Master/On_Off_ww/events" + str(10000+numb) +".fits"
                numb += 1
                primary_hdu = fits.PrimaryHDU()
                hdu_evt = fits.BinTableHDU(events.table)
                hdu_gti = fits.BinTableHDU(dataset.gti.table, name="GTI")
                hdu_all = fits.HDUList([primary_hdu, hdu_evt, hdu_gti])
                hdu_all.writeto(filename, overwrite=True)
            