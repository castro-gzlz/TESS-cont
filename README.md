# TESS-cont: The TESS contamination tool


<p>

  <img src="https://github.com/user-attachments/assets/b0409e7c-daa9-4dcf-b8bb-53a59b4878d2"
       align="left"
       width="320"
       style="margin-right: 20px; margin-bottom: 10px;">
  *TESS-cont* is a Python tool to **compute the flux contamination** from **nearby stars** in the **TESS aperture** of **any observed target**.
  
  The package **(1)** identifies the main contaminant *Gaia* sources, **(2)** computes their individual and total flux contributions to the aperture, and **(3)** determines whether any of these stars could be responsible for the observed **transit** and/or **activity** signals. The standard plotting output includes:  **(1)** heatmaps showing the flux fraction coming from the target star, and **(2)** pie charts highlighting the most contaminant stars inside (and outside!) the aperture.

</p>

*TESS-cont* is a **user-friendly** tool to **compute the flux contamination** from **nearby stars** in the **TESS aperture** of **any observed target**. 
The package **(1)** identifies the main contaminant *Gaia* sources, **(2)** computes their individual and total flux contributions to the aperture, 
and **(3)** determines whether any of these stars could be responsible for the observed **transit** and/or **activity** signals.





## Installation & Requirements

**Clone or download** this folder on your computer. All the dependencies are widely used, so you probably already have them installed. However, if this is not the case, you can enter the folder via the terminal and type
```
pip install -r requirements.txt
```
If you have any problems with the installation, you can drop an issue [here](https://github.com/castro-gzlz/TESS-cont/issues).

## Usage

*TESS-cont* is **very easy** to use. **It only requires the TIC/TOI** number (or Simbad name) of the target star to study, although it can also have a **flexible operation** based on different [optional parameters](#optional--optional-parameters).

#### Python through the terminal

You just need to have a configuration file (e.g. [*TOI-5005_S65.ini*](https://github.com/castro-gzlz/TESS-cont/blob/main/config/TOI-5005_S65.ini)) inside the [*config*](https://github.com/castro-gzlz/TESS-cont/tree/main/config) folder and type
```
python TESS-cont.py TOI-5005_S65.ini
```
In the [Configuration file](#configuration-file) section we describe **all the parameters** that can be used.

#### Python through Jupyter Notebook

*TESS-cont* can also be used via Jupyter Notebook ([TESS-cont.ipynb](https://github.com/castro-gzlz/TESS-cont/blob/main/TESS-cont.ipynb)). Similarly to the terminal version, you just need to select a configuration file and run all cells.

## Usage examples

### Example 1: A TIC number is all you need

In this example, we specify the **TIC number** of TOI-5005 (TIC 282485660) to quantify the flux contributions to its SPOC photometric aperture in **Sector 65**. 

```
python TESS-cont.py TOI-5005_S65.ini
```

![example1](https://github.com/user-attachments/assets/a7ca5fdb-b297-4f37-9a0a-a455454e2c75)


The left-hand plot is a **heatmap** indicating **the flux percentage from the target star** falling within each pixel. As we can see, inside the aperture, TOI-5005 contributes between **60%** and **90%** of the total flux, depending on the considered pixel. The disk sizes are scaled to the emitted fluxes. 

The right-hand plot is a **pie chart** representing the flux from the **target** and **nearby** stars **inside the photometric aperture**. Overall, **TOI-5005 contributes 84.4%** to the measured photometry, while **the remaining 15.6% comes from different nearby sources**. Interestingly, the second most contaminant star (**Star 2**) is a very bright source located **outside the TESS Target Pixel File (TPF)**. This example shows the **usefulness of *TESS-cont* to identify not evident contaminant sources that could be otherwise overlooked**. 

### Example 2: Possible origins of the transit (and variability) signals

We now study the star TOI-4479. As we can see, the field is **less crowded** (there are fewer white disks) but the SPOC photometry is **more contaminated**, with a total contaminant flux fraction of **30.4%** inside the aperture. 

```
python TESS-cont.py TOI-4479_S41.ini
```

![example2](https://github.com/user-attachments/assets/074bf91b-e2e3-4bdf-975c-c15568d6e165)


This star was found to have a planet candidate, TOI-4479.01, with a transit depth of 3471 ppm (parts per million) according to [ExoFOP](https://exofop.ipac.caltech.edu/tess/target.php?id=126606859), which was later confirmed to be a planet by [Esparza-Borges et al. (2022)](https://ui.adsabs.harvard.edu/abs/2022A%26A...666A..10E/abstract). We here use the *TESS-cont* **DILUTION** feature **to know whether we can discard the nearby contaminant sources as being the origin of the signal found** based on the TESS photometry alone. By introducing the measured transit depth, and selecting ```dilution_corr: True``` (since the SPOC transit depths are already corrected for dilution), we find that the observed transit signal **could have originated in the six most contaminant stars** with the following depths: 

**Star,Gaia_ID,transit_depth(%)**   
1,1841177816084707584,**1.74**  <br /> 
2,1841177816084683648,**6.10**  <br /> 
3,1841178129618984960,**8.256**   <br /> 
4,1841177717302124288,**11.30**  <br /> 
5,1841178232698201216,**24.35**  <br /> 
6,1841177682942388608,**27.69**

Since those values (**2-28%**) are lower than the unphysical threshold of 100%, we cannot discard that they could be the origin of the transit signal found based on the TESS data alone (e.g. [Castro-González et al. 2020](https://ui.adsabs.harvard.edu/abs/2020MNRAS.499.5416C/abstract); [de Leon et al. 2021](https://ui.adsabs.harvard.edu/abs/2021MNRAS.508..195D/abstract)). This example reveals the necessity of **higher-spatial resolution observations** (e.g. ground-based photometry or spectroscopy) to confidently **discard those contaminant sources**. 


## Computational cost: Approximate (fast) vs accurate (slower) method

By default, *TESS-cont* interpolates the PRFs to **all the pixels with *Gaia* sources, which might take several minutes**. To **streamline this process**, we have implemented a **fast approximate method** that can be easily activated through ```method_prf: approximate``` within the [OPTIONAL](#optional--optional-parameters) section. This method **interpolates the PRF only once** (in the middle of the TPF/FFI), before locating it in its corresponding position, assuming that its shape does not vary much across the nearby pixels. 

This approximate method **typically provides very similar results** to the default ```method_prf: accurate``` method, so it can be **useful to have a first hint** of the contamination level affecting your target (especially in highly crowded fields). However, **we encourage to use the accurate method for final analyses/publications**. We have coded *TESS-cont* so that the PRFs of stars in common pixels are only computed once, and hence **even in highly crowded fields the computational cost should not surpass ~5 minutes**. 

## Other uses, contamination metrics, and precautions

**Other uses**. *TESS-cont* can be also used to **generate custom apertures** based on the computed pixel-by-pixel contamination. We can select a certain threshold (e.g. 80%) of flux coming from the target star, and generate and save an aperture that meets such a threshold. This feature is currently not documented, but you can drop me a message and I'll be happy to help.

**Contamination metrics**. A **by-product** of the *TESS-cont* operation is the computation of the *CROWDSAP* and *FLFRCSAP* metrics. These are automatically saved in the [metrics.dat](https://github.com/castro-gzlz/TESS-cont/blob/main/metrics.dat) file. We encourage ensuring that there are no major differences with the official TESS metrics. If this were the case, it could be probably related to a discrepancy between the DR3 and DR2 *Gaia* catalogues. 

**Precautions**. By default, *TESS-cont* uses the *Gaia* DR3 catalogue. However, the TIC catalogue (and SPOC PDCSAP) is stacked to *Gaia* DR2. Therefore, to analyze the TESS PDCSAP photometry the **DR2 catalog should be selected** as an [OPTIONAL](#optional--optional-parameters) parameter: ```gaia_catalog: DR2```. **We highly encourage running *TESS-cont* based on the DR2 AND DR3 catalogues to ensure that there are no major differences**. If there were, it would be recommended to **use the DR3 contamination metrics** to correct the PDCSAP/SAP photometry from crowding as explained [here](https://heasarc.gsfc.nasa.gov/docs/tess/UnderstandingCrowding.html). 

## Configuration file


#### [MANDATORY] | Mandatory parameter

| Parameter  | Possible values | Description |
| ------------- | ------------- | ------------- |
| target | Any TIC, TOI, or star's name resoluble in Simbad | Target star identifier |


#### [OPTIONAL] | Optional parameters

| Parameter  | Possible values | Description |
| ------------- | ------------- | ------------- |
| sector | Any number| TESS sector. **Default**: first with observations |
| method_prf | accurate or approximate | Method to compute the PRFs. **Default**: accurate |
| search_radius | Any number | Search radius of *Gaia* sources (in arcsec). **Default**: 200 |
| n_sources | Any number | Contaminant sources to study individually. **Default**: 5 |
| gaia_catalog | DR2 or DR3 | Gaia catalog. **Default**: DR3 |
| target_name | Any name | Target name. **Default**: Target |
| plot_target_name | True or False | Specify the target name in the plot. **Default**: False|
| plot_all_gaia | True or False | Plot all *Gaia* nearby sources. **Default**: True |
| plot_percentages | True or False | Write the target's flux percentages. **Default**: True |
| loc_legend | best, upper left, etc | Location of the legend (heatmap plot). **Default**: best |
| scale_factor | Any number | Scale factor for the stars. **Default**: 4000 |
| scale_heatmap | natural or log | Scale of the heatmap color code. **Default**: natural |
| tpf_or_tesscut | tpf or tesscut | TPF or FFI tesscut. **Default**: tpf |
| cutout_size| Number, Number | Size of the FFI tesscut. **Default**: 11,11 |
| img_fmt | pdf, png, or pdfpng | Format of output images. **Default**: pdfpng |

#### [DILUTION] | Arguments for the **DILUTION** analysis

| Parameter  | Possible values | Description |
| ------------- | ------------- | ------------- |
| td | Any number| Measured transit depth in the target's aperture |
| td_unit | ppm, ppt, per, frac | Unit of the measured transit depth |
| dilution_corr | True or False | The measured depth is corrected for dilution (e.g. SPOC). **Default**: True  |

## Credits

If you use *TESS-cont*, please give credit to [this work](https://ui.adsabs.harvard.edu/abs/2024arXiv240918129C/abstract): 

```
@ARTICLE{2024A&A...691A.233C,
       author = {{Castro-Gonz{\'a}lez}, A. and {Lillo-Box}, J. and {Armstrong}, D.~J. and {Acu{\~n}a}, L. and {Aguichine}, A. and {Bourrier}, V. and {Gandhi}, S. and {Sousa}, S.~G. and {Delgado-Mena}, E. and {Moya}, A. and {Adibekyan}, V. and {Correia}, A.~C.~M. and {Barrado}, D. and {Damasso}, M. and {Winn}, J.~N. and {Santos}, N.~C. and {Barkaoui}, K. and {Barros}, S.~C.~C. and {Benkhaldoun}, Z. and {Bouchy}, F. and {Brice{\~n}o}, C. and {Caldwell}, D.~A. and {Collins}, K.~A. and {Essack}, Z. and {Ghachoui}, M. and {Gillon}, M. and {Hounsell}, R. and {Jehin}, E. and {Jenkins}, J.~M. and {Keniger}, M.~A.~F. and {Law}, N. and {Mann}, A.~W. and {Nielsen}, L.~D. and {Pozuelos}, F.~J. and {Schanche}, N. and {Seager}, S. and {Tan}, T. -G. and {Timmermans}, M. and {Villase{\~n}or}, J. and {Watkins}, C.~N. and {Ziegler}, C.},
        title = "{TOI-5005 b: A super-Neptune in the savanna near the ridge}",
      journal = {\aap},
     keywords = {techniques: photometric, techniques: radial velocities, planets and satellites: composition, planets and satellites: detection, planets and satellites: individual: TOI-5005 b, stars: individual: TOI 5005 (TIC 282485660), Astrophysics - Earth and Planetary Astrophysics},
         year = 2024,
        month = nov,
       volume = {691},
          eid = {A233},
        pages = {A233},
          doi = {10.1051/0004-6361/202451656},
archivePrefix = {arXiv},
       eprint = {2409.18129},
 primaryClass = {astro-ph.EP},
       adsurl = {https://ui.adsabs.harvard.edu/abs/2024A&A...691A.233C},
      adsnote = {Provided by the SAO/NASA Astrophysics Data System}
}
```
and add the following sentence within the acknowledgements section:

> This work made use of \texttt{TESS-cont} (\url{https://github.com/castro-gzlz/TESS-cont}), which also made use of \texttt{tpfplotter} \citep{2020A&A...635A.128A} and \texttt{TESS-PRF} \citep{2022ascl.soft07008B}.

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






