# TESS-cont: The TESS contamination tool

*TESS-cont* is a **user-friendly** Python tool to **quantify the flux contribution** from **nearby sources** in the **TESS photometry** of any observed star. The *TESS-cont* algorithm **(1)** identifies the main contaminant sources, **(2)** quantifies their individual and total flux contributions to the selected aperture (i.e. SPOC or custom), and **(3)** determines whether any of these stars could be the origin of the observed transit or variability signals. 

![Presentación sin título(2)-cropped](https://github.com/user-attachments/assets/59ef2a7f-f7db-4c9a-aa74-ec2ff71dc1a7)

The *TESS-cont* operation is based on building the Pixel Response Functions (PRFs) of nearby *Gaia* sources and computing their flux distributions across the Target Pixel Files (TPFs) or Full Frame Images (FFIs). A more detailed description can be found in Section 2.1.2 of [this work]().

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

### Example 1: The simplest case. A TIC number and a sector to study TOI-5005

In this example, we introduce the **TIC number** of TOI-5005 (282485660) to quantify the different flux contributions to its SPOC photometric aperture in Sector 65. By default, *TESS-cont* considers all nearby *Gaia* sources within a radius of 200 arcseconds, although this number can be customized in the [OPTIONAL] section.

```
python TESS-cont.py 
```

The left-hand plot is a **heatmap** indicating **the flux percentage from the target star** falling within each pixel. As we can see, inside the aperture, TOI-5005 contributes between 60% and 90% of the total flux, depending on the considered pixel. 

The right-hand plot is a **pie chart** representing the total flux percentage from the target and nearby stars **inside the photometric aperture**. Overall, **TOI-5005 contributes 84.4%** to the measured photometry, while **the remaining 15.6% comes from different nearby sources**. Of this 15.6% contaminant flux, 43% comes from small contributions of little-contaminant stars, and the remaining 57% comes from the five stars highlighted in both figures. These stars are numbered according to their total flux contribution. Interestingly, the second most contaminant star (Star#2) is a very bright source located outside the TESS TPF. This example illustrates the **usefulness of *TESS-cont* to identify contaminants that could be otherwise overlooked**. 


## Configuration file



