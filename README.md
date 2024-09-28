# TESS-cont: The TESS contamination tool

*TESS-cont* is a **user-friendly** Python tool to **quantify the flux fraction** coming from **nearby stars** in the **TESS photometric aperture** of **any observed target**. The package **(1)** identifies the main contaminant sources, **(2)** quantifies their individual and total flux contributions to the selected aperture (i.e. SPOC or custom), and **(3)** determines whether any of these stars could be the origin of the observed **transit** and **variability** signals. 

![Presentación sin título(2)-cropped](https://github.com/user-attachments/assets/59ef2a7f-f7db-4c9a-aa74-ec2ff71dc1a7)

The *TESS-cont* algorithm is based on building the Pixel Response Functions (PRFs) of nearby *Gaia* sources and computing their flux distributions across the TESS Target Pixel Files (TPFs) or Full Frame Images (FFIs). A more detailed description can be found in Section 2.1.2 of [this work]().

## Installation & Requirements

**Clone or download** this folder on your computer. All the dependencies are widely used, so you probably already have them installed. However, if this is not the case, you can enter the folder via the terminal and type
```
pip install -r requirements.txt
```
If you have any problems with the installation, you can drop an issue [here](https://github.com/castro-gzlz/TESS-cont/issues).

## Usage

*TESS-cont* is **easy** and **intuitive** to use. **It only requires the TIC number** of the target star to study, although it can also have a **flexible operation** based on different **optional parameters**.

#### Python through the terminal

You just need to create a configuration file *[my_config_file.ini](https://github.com/castro-gzlz/mr-plotter/blob/main/config/my_config_file.ini)* inside the [*config*](https://github.com/castro-gzlz/mr-plotter/tree/main/config) folder and then type

```
python TESS-cont.py my_config_file.ini
```
The file *my_config_file.ini* should contain all the necessary information to make the analysis, which will be saved into the *[output](https://github.com/castro-gzlz/mr-plotter/tree/main/output)* folder. In the [Configuration file](#configuration-file) section we describe **all the parameters** that can be used in the configuration files.


#### Python through Jupyter Notebook

*TESS-cont* can also be used via Jupyter Notebook (TESS-cont.ipynb). Similarly to the terminal version, you just need to select a configuration file (config_file) and run all cells.

## Usage examples

### Example 1: A TIC number is all you need

In this example, we specify the **TIC number** of TOI-5005 (TIC 282485660) to quantify the flux contributions to its SPOC photometric aperture in **Sector 65**. 

```
python TESS-cont.py TOI-5005_S65.ini
```

![example1](https://github.com/user-attachments/assets/494d54e3-7397-4812-99bd-a72c54aef721)


The left-hand plot is a **heatmap** indicating **the flux percentage from the target star** falling within each pixel. As we can see, inside the aperture, TOI-5005 contributes between **60%** and **90%** of the total flux, depending on the considered pixel. The disk sizes are scaled to the emitted fluxes. 

The right-hand plot is a **pie chart** representing the flux from the target and nearby stars **falling inside the photometric aperture**. Overall, **TOI-5005 contributes 84.4%** to the measured photometry, while **the remaining 15.6% comes from different nearby sources**. Interestingly, the second most contaminant star (**Star#2**) is a very bright source located **outside the TESS Target Pixel File (TPF)**. This example shows the **usefulness of *TESS-cont* to identify not evident contaminant sources that could be otherwise overlooked**. 

### Example 2: Possible origins of the transit (and variability) signals

We now study the star TOI-4479. As we can see, the field is **less crowded** (there are fewer white disks) but the SPOC photometry is **more contaminated**, with a total contaminant flux fraction of **30.4%** inside the aperture. 

```
python TESS-cont.py TOI-4479_S41.ini
```

![example2](https://github.com/user-attachments/assets/0825cde1-80e0-464f-a105-ef4c4f9529d1)


This star was found to have a planet candidate, TOI-4479.01, with a transit depth of 3471 ppm (parts per million) according to [ExoFOP](https://exofop.ipac.caltech.edu/tess/target.php?id=126606859), which was later confirmed to be a planet by [Esparza-Borges et al. (2022)](https://ui.adsabs.harvard.edu/abs/2022A%26A...666A..10E/abstract) with follow-up observations. We here use the *TESS-cont* **DILUTION** feature to know whether we can discard the nearby contaminant sources as being the origin of the signal found based on the TESS photometry alone. By introducing the measured transit depth, and selecting ```dilution_corr: True``` (since the TESS transit depths are already corrected for dilution), we find that the observed transit signal could have originated in the six most contaminant stars with the following transit depths: 

<sup> **Gaia_ID,TIC_ID,transit_depth(%)** <br /> 
1841177816084707584,No TIC ID found,**1.74** <br />
1841177816084683648,No TIC ID found,**6.01** <br />
1841178129618984960,No TIC ID found,**8.26** <br />
1841177717302124288,No TIC ID found,**11.27** <br />
1841178232698201216,No TIC ID found,**24.36** <br />
1841177682942388608,No TIC ID found,**27.66**  </sup>

Since those values (**2-28%**) are lower than the unphysical threshold of 100%, we cannot discard that they could be the origin of the transit signal found based on the TESS data alone. This example reveals the necessity of higher-spatial resolution observations (e.g. ground-based photometry or spectroscopy) to confidently discard those contaminant sources. 


### More examples: Custom apertures, SAP photometry correction, FFIs, and more. 

**Examples 1 and 2 represent the basic operation of *TESS-cont***. However, the package is also able to **create custom apertures to minimize contamination** from nearby sources, **correct SAP photometry from crowding** (which is very useful for targets without PDCSAP) starting from the TPFs an/or FFIs, and more! 

**During the next few days, I will be implementing more usage examples for these scenarios. Stay tuned and save the package if you find it useful for your research!**

## Configuration file


#### [MANDATORY] | Mandatory parameter

| Parameter  | Possible values | Description |
| ------------- | ------------- | ------------- |
| target | Any TIC, TOI, or star's name resoluble in Simbad | Target star identifier |


#### [OPTIONAL] | Optional parameters

| Parameter  | Possible values | Description |
| ------------- | ------------- | ------------- |
| sector | Any number| TESS sector to be analysed. Default: first with observations |
| n_sources | Any number | Number of contaminant sources to study individually. Default: 5 |
| gaia_catalog | DR2 or DR3 | Gaia catalog. Default: DR3 |
| target_name | Any name | Target name. Default: Target |
| plot_target_name | True or False | Specify the target name in the plots. Default: False|
| loc_legend | best, upper left, center, etc | Location of the legend (heatmap plot). Default: best |
| img_fmt | png or pdf | Format of output images. Default: pdf |

#### [DILUTION] | Arguments for the **DILUTION** analysis

| Parameter  | Possible values | Description |
| ------------- | ------------- | ------------- |
| td | Any number| Measured transit depth in the target star |
| dilution_corr | True or False | The introduced depth is corrected for dilution (e.g. TESS-SPOC products) or not. Default: True  |
| td_unit | ppm, ppt, per, frac | Unit of the introduced transit depth |

## Credits

## Credits

If you use *mr-plotter*, please give credit to [this work](): 

```
TBC
```
and add the following sentence within the acknowledgements section:

> This publication made use of \texttt{TESS-cont} (\url{https://github.com/castro-gzlz/TESS-cont}), which also used \texttt{tpfplotter} \citep{2020A&A...635A.128A} and \texttt{TESS\_PRF} \citep{2022ascl.soft07008B}.

```
@ARTICLE{2020A&A...635A.128A,
       author = {{Aller}, A. and {Lillo-Box}, J. and {Jones}, D. and {Miranda}, L.~F. and {Barcel{\'o} Forteza}, S.},
        title = "{Planetary nebulae seen with TESS: Discovery of new binary central star candidates from Cycle 1}",
      journal = {\aap},
     keywords = {planetary nebulae: general, techniques: photometric, binaries: general, Astrophysics - Solar and Stellar Astrophysics},
         year = 2020,
        month = mar,
       volume = {635},
          eid = {A128},
        pages = {A128},
          doi = {10.1051/0004-6361/201937118},
archivePrefix = {arXiv},
       eprint = {1911.09991},
 primaryClass = {astro-ph.SR},
       adsurl = {https://ui.adsabs.harvard.edu/abs/2020A&A...635A.128A},
      adsnote = {Provided by the SAO/NASA Astrophysics Data System}
}
```

```
@software{2022ascl.soft07008B,
       author = {{Bell}, Keaton J. and {Higgins}, Michael E.},
        title = "{TESS\_PRF: Display the TESS pixel response function}",
 howpublished = {Astrophysics Source Code Library, record ascl:2207.008},
         year = 2022,
        month = jul,
          eid = {ascl:2207.008},
       adsurl = {https://ui.adsabs.harvard.edu/abs/2022ascl.soft07008B},
      adsnote = {Provided by the SAO/NASA Astrophysics Data System}
}
```






