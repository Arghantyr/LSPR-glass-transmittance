# What is lspr?

_lspr_ is a Python 3 scientific package for modelling UV-Vis spectral behaviour of metallic nanoparticle-dieletric composites.

# Dependencies

To run _lspr_ smoothly you will need:
1. Python 3.4+
2. Pandas 1.0.3+
3. NumPy 1.18.1+
4. SciPy 1.4.1+
5. Statsmodels 0.11.0+
6. Matplotlib 3.1.3+

# What does lspr enable?

_lspr_ enables one to approximate some spectroscopic features of metallic nanoparticle-dielectric composites, e.g. transmittance, absorbance, scattering cross-section. The model utilizes dipole approximation of Mie solution to Maxwell equations, described in more detail in article by Olson et al.<sup>1</sup>

To perform a simulation, one needs to construct a composite, so matrix and metallic nanoparticles, so for each a set of refractive index and extinction coefficient data for simulated spectral range are necessary. The model **does not take into account possible interactions at the nanoparticle/matrix interface**, i.e. the composite absorbance is a simple sum of matrix and nanoparticle absorbances calculated for pristine materials.



<sup>1</sup> _**Optical Characterization of Single Plasmonic Nanoparticles**_ - J. Olson, S. Dominguez-Medina, A. Hoggard, L.-Y. Wang, W.-Sh. Chang, S. Link; *Chem Soc Rev. 2015 January 7; 44(1): 40-57. doi: 10.1039/c4cs00131a.*

# How to get it running?

## PyPI

```
pip install lspr
```

## Anaconda

```
conda install -c angantyr lspr
```

# Usage

More examples on how to use the lspr and specific methods can be found in the documentation.
The transmittance of 31nm silver nanoparticles homogeneously distributed in vacuum with concentration 5.8 &middot; 10<sup>13</sup> cm<sup>-3</sup>

```
import lspr

# Get the csv file with n,k values for silver in range 200-1000nm
# spaced evenly every 1nm
Ag = lspr.OpticalMaterial(name = "Silver",source = "n-k_Silver_200-1000nm_(1nm).csv")

# Get n,k values for vacuum
Vacuum = lspr.OpticalMaterial(name = "Vacuum",
                         source = "n-k_data_for_vacuum_200-1000nm_(1nm).csv")

# Create Ag nanoparticles with mean size of 31nm, standard deviation of 1nm
# and concentration of 5.8e13 cm^-3
nAg = lspr.Nanoparticle(name = "Silver nanoparticles", material=Ag,
                   size=31, concentration=5.8e13)

# Set the transmittance measurement of a 1 micron thick sample of the composite
(silver nanoparticles dispersed in vacuum); use "T_corr" for reflectance corrected result
Transmittance = lspr.SpectralMeasurement("Measure the transmittance",
                                    matrix=Vacuum, nanoparticle=nAg,
                                    thickness=3).plot_results("T_corr")
```

![The result is a preformatted plot of the reflectance corrected transmittance](https://github.com/Arghantyr/LSPR-glass-transmittance/blob/master/example.jpg)

# License
lspr is distributed under [MIT](https://choosealicense.com/licenses/mit/) license.

# Support
If issue is spotted please open an issue on the [GitHub repo of the project](https://github.com/Arghantyr/LSPR-glass-transmittance/issues). For changes, upgrades and simillar refer to the [project's wiki site](https://github.com/Arghantyr/LSPR-glass-transmittance/wiki).
