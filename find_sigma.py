

import numpy as np
import matplotlib.pyplot as plt

import os

import astropy
import astropy.units as u
from astropy.coordinates import Angle,SkyCoord
from regions import CircleSkyRegion, RectangleSkyRegion
from IPython.display import display
from gammapy.data import DataStore, Observation,observatory_locations, Observations
from gammapy.datasets import (
    Datasets,
    FluxPointsDataset,
    SpectrumDataset,
    MapDataset,
    SpectrumDatasetOnOff,
    
)
from gammapy.utils.regions import make_orthogonal_rectangle_sky_regions
from gammapy.estimators import FluxPointsEstimator, FluxPoints, FluxProfileEstimator
from gammapy.makers import (
    ReflectedRegionsBackgroundMaker,
    SafeMaskMaker,
    SpectrumDatasetMaker,MapDatasetMaker
)
from gammapy.irf import load_cta_irfs
from gammapy.astro.darkmatter import (
    DarkMatterAnnihilationSpectralModel,
    JFactory,
    PrimaryFlux,
    profiles,
)
from gammapy.maps import MapAxis, RegionGeom, WcsGeom, WcsNDMap
from gammapy.modeling import Fit
from gammapy.modeling.models import (
    ExpCutoffPowerLawSpectralModel,TemplateSpatialModel, PowerLawNormSpectralModel,
    SkyModel,
    create_crab_spectral_model,
    Models
)
import gammapy
from pathlib import Path
def create_dm(massDM):
    

    
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
    jfactor = jfact.sum()
    j_Einasto = 7.1*10**22*u.Unit("GeV2/cm5")
    
    DarkMatterAnnihilationSpectralModel.THERMAL_RELIC_CROSS_SECTION
    modelDM = DarkMatterAnnihilationSpectralModel(mass=massDM,
    channel=channel, jfactor=j_Einasto)
    return SkyModel(
        spatial_model=spatial_model, spectral_model=modelDM,
        name="model-simu"
    )
def create_masks():
    target_position = SkyCoord(0.0,0.0, unit="deg", frame="galactic")
    exclusion_region = CircleSkyRegion(
        center=SkyCoord(0.85, 0.1, unit="deg", frame="galactic"),
        radius=0.1 * u.deg,
    )
    exclusion_region2 = CircleSkyRegion(
        center=SkyCoord(358.5, -1.10, unit="deg", frame="galactic"),
        radius=1.2* u.deg,
    )
    exclusion_region3 = RectangleSkyRegion(center= SkyCoord(0.0,0.0, unit="deg", frame="galactic"), width= 8.0*u.deg,
                                          height = 1.2*u.deg )
    skydir = target_position
    geom = WcsGeom.create(
        npix=(150, 150), binsz=0.05, skydir=skydir, proj="TAN", frame="galactic"
    )
    exclusion_mask = ~geom.region_mask([exclusion_region])
    mask2 = ~geom.region_mask([exclusion_region2])
    mask3 = ~geom.region_mask([exclusion_region3])
    
    exclusion_mask &= mask2
    exclusion_mask &= mask3
    
    
    return exclusion_mask
path =  Path(os.environ["CTADATA"])/"data/baseline/gc"

irf_path = Path(os.environ["CALDB"])/"South_z20_0.5h/irf_file.fits"

s_v = []
s_v2 = []
paths = list(path.rglob("*.fits"))

datastore = DataStore.from_events_files(paths,irf_path)


irfs = load_cta_irfs(irf_path)
#Choose observational id from table given by read_fits.py
#obs_ids = [311000,310040,310969]
obs_ids = np.load("obs_ids_dc.npy")
observations = datastore.get_observations(obs_ids)

target_position = SkyCoord(0.0,0.0, unit="deg", frame="galactic")
on_region_radius = Angle("0.2 deg")
on_region = CircleSkyRegion(center=target_position, radius=on_region_radius)

exclusion_mask = create_masks()
import time

t0 = time.time()

t1 = time.time()
errs = []
total = t1-t0
n = 500
mass = np.logspace(-1,2,n)
m_v = []
for i in range(n):
    m = mass[i]
    e_low = 0.1*m
    e_up = 1.2*m
    energy_axis = MapAxis.from_energy_bounds(
         e_low, e_up, nbin=10, per_decade=True, unit="TeV", name="energy"
     )
    energy_axis_true = MapAxis.from_energy_bounds(
         e_low/2.0, e_up*2.0, nbin=12, per_decade=True, unit="TeV", name="energy_true"
     )
     
    geom = RegionGeom.create(region=on_region, axes=[energy_axis])
    dataset_empty = SpectrumDataset.create(geom=geom, energy_axis_true=energy_axis_true)
    
    dataset_maker = SpectrumDatasetMaker(
        containment_correction=True, selection=["counts", "exposure", "edisp"]
    )
    bkg_maker = ReflectedRegionsBackgroundMaker(exclusion_mask=exclusion_mask)
    #bkg_maker = ReflectedRegionsBackgroundMaker()
     
    safe_mask_masker = SafeMaskMaker(methods=["aeff-max"], aeff_percent=10)
     
    datasets = Datasets()
     
    for obs_id, observation in zip(obs_ids, observations):
        dataset = dataset_maker.run(dataset_empty.copy(name=str(obs_id)), observation)
        dataset_on_off = bkg_maker.run(dataset, observation)
        dataset_on_off = safe_mask_masker.run(dataset_on_off, observation)
        datasets.append(dataset_on_off)
     
    dataset_stacked = Datasets(datasets).stack_reduce()
    dataset_stacked.models = create_dm(m*u.Unit("TeV"))
    fit = Fit()
    result_stacked = fit.run([dataset_stacked])


    par = dataset_stacked.models.parameters.free_parameters[0]
    b = 60
    par.scan_n_values = b
    profile = fit.stat_profile(datasets=dataset_stacked, parameter=par)
    TS = profile["stat_scan"] - result_stacked.total_stat
    # Compute the corresponding statistical significance 
    x = profile["model-simu.spectral.scale_scan"]
    errs.append(result_stacked.models[0].parameters[0].error)
    
    l = int(b/2)
    upper_TS = abs(2.71-TS[l:])
    upper_x = x[l:]
    for a in range(len(upper_TS)):
        if upper_TS[a] == min(upper_TS):
            scale_int = a
    s_v.append(upper_x[scale_int])
    s_v2.append(upper_x[-1])
    m_v.append(m)
t1 = time.time()

total = t1-t0
print(total)
np.save("m_v_50.npy",m_v)
np.save("s_v_50.npy",s_v)
np.save("errorst_50.npy", errs)
np.save("s_v2_50.npy",s_v2)
