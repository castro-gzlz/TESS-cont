# TESS-cont: The TESS contamination tool

*TESS-cont* is a **user-friendly** Python tool to **quantify the flux contribution** from **nearby sources** in the **TESS photometry** of any observed star. The *TESS-cont* algorithm **(1)** **identifies** the main contaminant sources, **(2)** **quantifies** their individual and total contributions to the selected aperture (i.e. SPOC or custom), and **(3)** **determines** whether any of these sources could be the origin of the observed transit or variability signals. 

![Presentación sin título(2)-cropped](https://github.com/user-attachments/assets/59ef2a7f-f7db-4c9a-aa74-ec2ff71dc1a7)

The *TESS-cont* operation is based on building the Pixel Response Functions (PRFs) of nearby *Gaia* sources and computing their flux distributions across the Target Pixel Files (TPFs) or Full Frame Images (FFIs). A more detailed description can be found in Section 2.1.2 of [this work]().

## Installation & Requirements

Just **clone or download** this folder on your computer. All the dependencies are widely used, so you probably already have them installed. However, if this is not the case, you can just enter the folder via the terminal and type
```
pip install -r requirements.txt
```
If you have any problems with the installation, you can drop an issue [here](https://github.com/castro-gzlz/mr-plotter/issues).
