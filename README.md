# TESS-cont: The TESS contamination tool

*TESS-cont* is a **user-friendly** Python tool to **quantify the flux contribution** from **nearby sources** in the **TESS photometry** of any observed star. The *TESS-cont* algorithm **(1)** **identifies** the main contaminant sources, **(2)** **quantifies** their individual and total contributions to the selected aperture (i.e. SPOC or custom), and **(3)** **determines** whether any of these sources could be the origin of the observed transit or variability signals. 

![Presentación sin título(1)-cropped](https://github.com/user-attachments/assets/b1638eee-2012-41e1-be04-014458d30a2b)


To do so, *TESS-cont* builds the Pixel Response Functions (PRFs) of each nearby *Gaia* source and computes the flux distributions across the Target Pixel Files (TPFs) or Full Frame Images (FFIs). A more detailed explanation of its operation can be found in Section 2.1.2 of [this work]().
