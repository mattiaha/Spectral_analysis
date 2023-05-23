

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


path =  Path(os.environ["CTADATA"])/"data/baseline/gc"

irf_path = Path(os.environ["CALDB"])/"South_z20_0.5h/irf_file.fits"


paths = list(path.rglob("*.fits"))

datastore = DataStore.from_events_files(paths,irf_path)


irfs = load_cta_irfs(irf_path)

plt.show()

#Choose observational id from table given by read_fits.py
obs_ids = [311000,310040,310969]


#use obs_ids2 if npy-file is generated from read_fits.py 
#obs_ids2 = np.load("obs_ids_dc.npy")
observations = datastore.get_observations(obs_ids)
#define on region
target_position = SkyCoord(0.0,0.0, unit="deg", frame="galactic")
on_region_radius = Angle("0.2 deg")
on_region = CircleSkyRegion(center=target_position, radius=on_region_radius)

exclusion_region = CircleSkyRegion(
    center=SkyCoord(0.85, 0.1, unit="deg", frame="galactic"),
    radius=0.1 * u.deg,
)
exclusion_region2 = CircleSkyRegion(
    center=SkyCoord(358.5, -1.10, unit="deg", frame="galactic"),
    radius=1.2* u.deg,
)
exclusion_region3 = RectangleSkyRegion(center= SkyCoord(0.0,0.0, unit="deg", frame="galactic"), width= 8.0*u.deg,
                                       height = 0.6*u.deg )
skydir = target_position
geom = WcsGeom.create(
    npix=(150, 150), binsz=0.05, skydir=skydir, proj="TAN", frame="galactic"
)
plt.figure(5)
exclusion_mask = ~geom.region_mask([exclusion_region])
mask2 = ~geom.region_mask([exclusion_region2])
mask3 = ~geom.region_mask([exclusion_region3])



exclusion_mask &= mask2
exclusion_mask &= mask3
exclusion_mask.plot()


energy_axis = MapAxis.from_energy_bounds(
    1.0, 12, nbin=12, per_decade=True, unit="TeV", name="energy"
)
energy_axis_true = MapAxis.from_energy_bounds(
    0.5, 24, nbin=30, per_decade=True, unit="TeV", name="energy_true"
)

geom = RegionGeom.create(region=on_region, axes=[energy_axis])
dataset_empty = SpectrumDataset.create(geom=geom, energy_axis_true=energy_axis_true)

dataset_maker = SpectrumDatasetMaker(
    containment_correction=True, selection=["counts", "exposure", "edisp"]
)
bkg_maker = ReflectedRegionsBackgroundMaker(exclusion_mask=exclusion_mask)

safe_mask_masker = SafeMaskMaker(methods=["aeff-max"], aeff_percent=10)

datasets = Datasets()

for obs_id, observation in zip(obs_ids, observations):
    dataset = dataset_maker.run(dataset_empty.copy(name=str(obs_id)), observation)
    dataset_on_off = bkg_maker.run(dataset, observation)
    dataset_on_off = safe_mask_masker.run(dataset_on_off, observation)
    datasets.append(dataset_on_off)

from gammapy.visualization import plot_spectrum_datasets_off_regions
observations[0].events.peek()



profile = profiles.EinastoProfile(r_s=20 * u.kpc , alpha=0.17,
rho_s=0.081* u.Unit("GeV / cm3"))
profiles.DMProfile.DISTANCE_GC = 8.5*u.kpc
profiles.DMProfile.LOCAL_DENSITY = 0.39 * u.Unit("GeV / cm3")


profile.scale_to_local_density()
print(profile.parameters.value)

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
dataset_stacked = Datasets(datasets).stack_reduce()
plt.figure()
ax = exclusion_mask.plot()
on_region.to_pixel(ax.wcs).plot(ax=ax, edgecolor="k")
plot_spectrum_datasets_off_regions(ax=ax, datasets=datasets)

dataset_stacked.models = dm_model

print(dm_model)
stacked_fit = Fit()
result_stacked = stacked_fit.run([dataset_stacked])

datasets.models = dm_model

for i in range(len(obs_ids)):
    fig = plt.figure()
    ax_spectrum, ax_residuals = datasets[i].plot_fit()
    ax_spectrum.set_ylim(0.1, 500)
    ax_spectrum.set_xlim(0.5, 20)
    ax_residuals.set_xlim(0.5,20)
    datasets[i].plot_masks(ax=ax_spectrum)

plt.figure()
ax_spectrum, ax_residuals = dataset_stacked.plot_fit()
#ax_spectrum.set_ylim(0.1, 500)
#ax_spectrum.set_xlim(0.5, 20)
#ax_residuals.set_xlim(0.5,20)

dataset_stacked.plot_masks(ax=ax_spectrum)
plt.show()
e_min, e_max = 1.5, 10
energy_edges = np.geomspace(e_min, e_max, 8) * u.TeV
# creation of the boxes and axis
fpe = FluxPointsEstimator(
    energy_edges=energy_edges, source="model-simu", selection_optional="all"
)
flux_points = fpe.run(datasets=dataset_stacked)
display(flux_points.to_table(sed_type="dnde", formatted=True))

flux_points_dataset = FluxPointsDataset(data=flux_points, models=dm_model)
flux_points_dataset.plot_fit()
