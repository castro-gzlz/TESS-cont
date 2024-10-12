#!/usr/bin/env python
# coding: utf-8

# In[ ]:


#@|++++++++++++++++++++++++++++++++++++++++++++++++
#@|-----TESS-cont: The TESS contamination tool-----
#@|++++++++++++++++++++++++++++++++++++++++++++++++


# In[ ]:


import os
import PRF
import sys
import warnings
import argparse
import numpy as np
import pandas as pd
from tqdm import tqdm
from astropy import wcs
import lightkurve as lk
from astropy.io import ascii
from matplotlib import patches
import matplotlib.pyplot as plt
from astropy.table import Table
from colorsys import hsv_to_rgb
from astroquery.mast import Catalogs
from configparser import ConfigParser
import matplotlib.gridspec as gridspec
from matplotlib.colorbar import Colorbar
import astropy.visualization as stretching
from matplotlib.patches import ConnectionPatch
from astropy.coordinates import SkyCoord, Angle
from matplotlib.collections import PathCollection
from matplotlib.legend_handler import HandlerPathCollection
from astropy.visualization.mpl_normalize import ImageNormalize


# In[ ]:


#@|Uncomment this line for the tess-cont.ipynb version
#@|#####--Configuration file--#########
#config_file = 'TOI-5005_S65.ini'
#@|####################################


# In[ ]:


#@|++++Uncomment these lines for the TESS-cont.py version+++++++
parser = argparse.ArgumentParser()
parser.add_argument('config_file')
args = parser.parse_args()
config_file = args.config_file

warnings.filterwarnings("ignore")
print("")
print("@|---------------------@|")
print("        TESS-cont       ")
print("@|---------------------@|")
print("")


# In[ ]:


#@|Read the configuration file 'config.ini' into a config_object
config_path = 'config/'+config_file
config_object = ConfigParser()
config_object.read(config_path)


# In[ ]:


#@|-----------------------------------------
#@|------------Sections---------------------
#@|-----------------------------------------

MANDATORY = config_object['MANDATORY']
try:
    OPTIONAL = config_object['OPTIONAL']
except:
    pass
try:
    APERTURE = config_object['APERTURE']
except:
    pass
try:
    DILUTION = config_object['DILUTION']
except:
    pass


# In[ ]:


#@|-------------------
#@|MANDATORY argument
#@|-------------------
target = MANDATORY['target']


# In[ ]:


#@|-------------------
#@|OPTIONAL arguments
#@|-------------------

#@|TESS sector to be analysed
try:
    sector = OPTIONAL['sector']
except:
    pass

#@|target name (for the plots and output folders)
try:
    target_name = OPTIONAL['target_name']
except:
    target_name = str(target)
    
#@|search radius for Gaia DR3 targets
try:
    search_radius = float(OPTIONAL['search_radius'])
except:
    search_radius = 200 #arcsec
    
try:
    tpf_or_tesscut = OPTIONAL['tpf_or_tesscut']
except:
    tpf_or_tesscut = 'tpf'
    
try:
    cutout_size = OPTIONAL['cutout_size']
    cutout_size = (int(cutout_size.split(',')[0]),                    int(cutout_size.split(',')[0]))
except:
    cutout_size = (11,11) #@|similar to a tpf
    
try:
    method_prf = OPTIONAL['method_prf']
except:
    method_prf = 'accurate'
    
#@|legend location of the heatmap

try:
    loc_legend = OPTIONAL['loc_legend']
except:
    loc_legend = 'best'
    
#@|number of contaminant sources to consider individually
#@|(for the pie chart and heatmap)    

try:
    n_sources = int(OPTIONAL['n_sources'])
except:
    n_sources = 5
    
try:
    img_fmt = OPTIONAL['img_fmt']
except:
    img_fmt = 'pdfpng'
    
try:
    save_metrics = OPTIONAL['save_metrics'] == 'True'
except:
    save_metrics = True
    
try:
    gaia_catalog = OPTIONAL['gaia_catalog']
except:
    gaia_catalog = 'DR3'
    
try:
    plot_target_name = OPTIONAL['plot_target_name'] == 'True'
except:
    plot_target_name = False
    
    
#@|---------(HEATMAP arguments)-------------

#@|plot target star 
try:
    plot_target = OPTIONAL['plot_target'] == 'True'
except:
    plot_target = True

#@|plot main contaminant sources
try:
    plot_main_contaminants = OPTIONAL['plot_main_contaminants'] == 'True'
except:
    plot_main_contaminants = True

#@|plot all Gaia sources
try:
    plot_all_gaia = OPTIONAL['plot_all_gaia'] == 'True'
except:
    plot_all_gaia = True
    
    
#@|Overplot the flux ratio from the target star

try:
    plot_percentages =  OPTIONAL['plot_percentages'] == 'True'
except:
    plot_percentages = True
    
#@|scale factor for the stars. Disk area scales with flux emission
try:
    scale_factor = float(OPTIONAL['scale_factor'])
except:
    scale_factor = 4000
    
    
#@|-----------------------------------------------------------------------------------------------
#@|APERTURE arguments | not documented. For details, please email me at acastro@cab.inta-csic.es
#@|------------------------------------------------------------------------------------------------  

#@|aperture: pipeline, threshold_target_flux, or threshold_median_flux
try:
    aperture = APERTURE['aperture']
except:
    aperture = 'pipeline'
    
#@|minimum flux coming from the target star in a given pixel (only if aperture: threshold_target_flux )
try:
    threshold_target = float(APERTURE['threshold_target'])
except:
    threshold_target = 0.7 #@|(values between 0-1)

#@|sigmas above the median tpf flux (only if aperture: threshold_median_flux)
try:
    threshold_median = float(APERTURE['threshold_median'])
except:
    threshold_median = 3  #sigma

#@|save the aperture in a .csv file.
try:
    save_aper = APERTURE['save_aper'] == 'True'
except:
    save_aper = False
    
#@|-------------------------
#@|DILUTION arguments
#@|--------------------------

#@|observed transit depth of the planet candidate    
try:
    td = float(DILUTION['td'])
    transit_depth_analysis = True
except:
    transit_depth_analysis = False
    
#@|is the transit depth already corrected for dilution (e.g. SPOC td)?    
try:
    dilution_corr = DILUTION['dilution_corr'] == 'True'
except:
    dilution_corr = True

#@|transit depth unit    
try:
    td_unit = DILUTION['td_unit']  #ppm, ppt, per, frac
except:
    td_unit = 'frac'


# In[ ]:


#@|############################
#@|we create the output folder
#@|############################
try:
    os.mkdir('output/'+target_name)
except:
    pass


# In[ ]:


#@|------------------------------------------------------------------------------------------
#@|we download the tpf (or tesscut) of the target to build the pixel response function (PRF).
#@|------------------------------------------------------------------------------------------

#@|This will provide us with the following information:
#@|the camera (tpf.camera), and the CCD (tpf.ccd) of the tpf.
#@|the columns and rows of reference (tpf.column, tpf.row). Note that the PRF is different in each pixel.
#@|the aperture mask of the SPOC pipeline (in case it has one): tpf.pipeline_mask.
#@|the shape of the tpf in order to create the PRF with the same dimensions: tpf.shape[1:3].

#@|++++++++++++++++++++++++++++++++++++
#@|download the Target Pixel File (TPF)
#@|++++++++++++++++++++++++++++++++++++
if tpf_or_tesscut == 'tpf':
    try:
        #search_result = lk.search_targetpixelfile('TIC '+str(tic), sector = int(sector))
        search_result = lk.search_targetpixelfile(str(target), sector = int(sector))
    except NameError: search_result = lk.search_targetpixelfile('TIC '+str(tic))
    tpf = search_result.download()
    tic = tpf.targetid
    if len(search_result) == 0:
        try:
            print(f'Error: There is not a target pixel file for TIC {tic} (sector {sector}).            You can try with a tesscut on the FFIs!')
        except NameError: print(f'Error: There is not a target pixel file for TIC {tic}.            You can try with a tesscut on the FFIs!')
        sys.exit()

#@|+++++++++++++++++++++++++++++++++++++++++++++++++++++++
#@|download a tesscut on the TESS full-frame images (FFIs)
#@|+++++++++++++++++++++++++++++++++++++++++++++++++++++++
if tpf_or_tesscut == 'tesscut':
    try:
        search_result = lk.search_tesscut(str(target), sector = int(sector))
    except NameError: lk.search_tesscut(str(target))
    tpf = search_result.download(cutout_size = cutout_size)
    tic = tpf.targetid
    if len(search_result) == 0:
        try:
            print(f'Error: An error has occoured downloading the tesscut of TIC {tic} (sector {sector}).')
        except NameError: print(f'Error: An error has occoured downloading the tesscut of TIC {tic}.')
        sys.exit()


# In[ ]:


#@|modified from tpfplotter (J.Lillo-Box)
#@|https://github.com/jlillo/tpfplotter

def get_gaia_data(ra, dec, search_radius=250):
    c1 = SkyCoord(ra, dec, frame='icrs', unit='deg')
    from astroquery.vizier import Vizier
    Vizier.ROW_LIMIT = -1
    if gaia_catalog == 'DR3':
        gaia_cat, catID = "I/355/gaiadr3", gaia_catalog
    if gaia_catalog == 'DR2':
        gaia_cat, catID = "I/345/gaia2", gaia_catalog

    result = Vizier.query_region(c1, catalog=[gaia_cat],
                                 radius=Angle(search_radius, "arcsec"))
    try:
        result = result[gaia_cat]
    except:
        print('This target is not in Gaia '+catID)
        print('Exiting without finishing...')
        sys.exit()

    no_targets_found_message = ValueError('Either no sources were found in the query region '
                                          'or Vizier is unavailable')
    too_few_found_message = ValueError('No sources found closer than 1 arcsec to TPF coordinates')
    if result is None:
        raise no_targets_found_message
    elif len(result) == 0:
        raise too_few_found_message
        
    return result 


# In[ ]:


#@|---------------------------------------------------------------------
#@|we build an astroy table with all the Gaia sources within a radius R.
#@|---------------------------------------------------------------------
print(f'Searching for nearby Gaia {gaia_catalog} stars ...')
table = get_gaia_data(tpf.ra, tpf.dec, search_radius = search_radius)
table = table[~table['Gmag'].value.mask]  #@|we discard those targets with a 'nan' G magnitude
print(f'There are {len(table)} stars surrounding {target_name} within a radius of {search_radius} arcsec')


# In[ ]:





# In[ ]:


#@|from tpfplotter (J.Lillo-Box)
#@|https://github.com/jlillo/tpfplotter
def get_dr2_id_from_tic(tic):
    # Get the Gaia sources
    result = Catalogs.query_object('TIC'+tic, radius=.005, catalog="TIC")
    IDs = result['ID'].data.data
    k = np.where(IDs == tic)[0][0]
    GAIAs = result['GAIA'].data.data
    #Gaiamags = result['GAIAmag'].data.data

    GAIA_k = GAIAs[k]
    #Gaiamag_k = Gaiamags[k]

    if GAIA_k == '':
        GAIA_k = np.nan
    #return GAIA_k, Gaiamag_k
    return int(GAIA_k)


# In[ ]:


#@|---------------------------------------------------------------------
#@|we obtain which index in 'table' corresponds to our target (idx_target)
#@|---------------------------------------------------------------------
if tpf_or_tesscut == 'tesscut':
    gaia_dr3_target = get_dr2_id_from_tic(str(tic[4:]))
if tpf_or_tesscut == 'tpf':
    gaia_dr3_target = get_dr2_id_from_tic(str(tic))
idx_target = np.where(table['Source'] == gaia_dr3_target)[0][0]
#print(idx_target)


# In[ ]:





# In[ ]:


#@|we create a new column within 'table' called 'flux', in which we estimate, based on their G magnitudes, 
#@|the flux of each star with respect to the flux of our target: flux = f_star/f_target. 
#@|to do so, we used the following relation:
#@|f_star/f_target = 100**((m_target-m_star)/5)    
#@|from http://burro.case.edu/Academics/Astr221/Light/magscale.html

flux_col = np.zeros(len(table))
for i in range(len(table)):
    flux_col[i] = 100**((table['Gmag'][idx_target] - table['Gmag'][i]) / 5)
table['flux'] = flux_col


# In[ ]:


#@|we get the coordinates (RA, Dec) of the Gaia DR3 nearby targets ('coords'), and we use the wcs to convert 
#@|them into pixel coordinates inside the tpf (pixel_coords).
#@|the pixel coordinates are used to:
#@|1) get the colnum and rownum of each source (to get the proper PRF)
#@|2) overplot the sources over the PRF plot.
print(f'Extracting the Gaia {gaia_catalog} coordinates of the nearby targets and converting them into pixel coordinates ... ')
coords = []
pixel_coords = []

#@|proper motion correction 

if gaia_catalog == 'DR3':
    t_reference =  2457389            # 2016-01-01 12:00:00.000 | Lindegren et al. (2021)
if gaia_catalog == 'DR2':
    t_reference = 2457389 - 182.625   # 2015-07-02 21:00:00.000 | Lindegren et al. (2021)
    
t_inc = (tpf.time[0].jd - t_reference) / 365  #year

for i in tqdm(range(len(table))):
    
    #coords_ac = SkyCoord(table[i]['RA_ICRS'], table[i]['DE_ICRS'], unit="deg") # OLD --- no pm correction ---
    coords_ac = SkyCoord(table[i]['RA_ICRS']+np.nan_to_num(table['pmRA'].value.data[i]) / 3600000 * t_inc,                         table[i]['DE_ICRS']+np.nan_to_num(table['pmDE'].value.data[i]) / 3600000 * t_inc,                         unit="deg")  # defaults to ICRS frame | includes pm correction
    pixel_coords_ac = wcs.utils.skycoord_to_pixel(coords_ac, tpf.wcs, origin = 0, mode='all')
    coords.append(coords_ac)
    pixel_coords.append(pixel_coords_ac)


# In[ ]:


#@|-------------------------
#@|----we build the PRF-----
#@|-------------------------

resampled = np.zeros(tpf.shape[1:3]) #@|TOTAL PRF. It will contain the flux contributions of all Gaia sources.
resampled_list = []  #@|list of PRFs. Each item is the PRF of each Gaia target in 'table'.

#@|these are fixed values common to all targets
cam = tpf.camera
ccd = tpf.ccd
sector = tpf.sector

#print(f'Building the Point Response Functions (PRF) of each Gaia target ... ({method_prf} method)')

if method_prf == 'approximate':
    #@|In the approximate method, we estimate the prf in the middle of the TPF ONLY ONCE, and use
    #@|the obtained distribution for all targets (assuming that its shape won't change much)
    prf = PRF.TESS_PRF(cam,ccd,sector,tpf.column+tpf.shape[2]/2, tpf.row+tpf.shape[1]/2)
    print('PRF built in the middle of the TPF (approximate method)')
    
    #@|we locate the computed PRF in each star's location
    print('Resampling for the heatmap plot ...')
    for i in tqdm(range(len(table))): 
        
        #@|we resample the PRF into the shape of the original tpf
        try:
            resampled_ac = prf.locate(pixel_coords[i][0], pixel_coords[i][1], table['flux'][i], tpf.shape[1:3])
        except:
            #@|this exception occours when there are Gaia sources too far away 
            resampled_ac = np.zeros(tpf.shape[1:3]) 
        resampled_list.append(resampled_ac)
        resampled = resampled + resampled_ac
    
if method_prf == 'accurate':
    #@|In the accurate method, we estimate the PRF for each individual target in each pixel location
    #@|This method is generally slower than the approximate method
    colrow_array = np.zeros((len(table), 2))
    for i in range(len(table)): 
    
        colnum = int(tpf.column + pixel_coords[i][0])
        rownum = int(tpf.row + pixel_coords[i][1])

        colrow_array[i] = colnum, rownum
        
    ###############
    
    colrow_unique = np.unique(colrow_array, axis = 0)  
    prf_array = np.zeros(len(colrow_unique)).astype('object')

    print('Building the PRFs in each TESS pixel with nearby Gaia sources ... (this might take a while) ')
    for i in tqdm(range(len(prf_array))): 
        prf_array[i] = PRF.TESS_PRF(cam,ccd,sector, colrow_unique[i][0], colrow_unique[i][1])
        
    ################
    
    #@|we locate the computed PRFs in each star's location
    print('Resampling for the heatmap plot ...')
    for i in tqdm(range(len(table))): 
        colnum = int(tpf.column + pixel_coords[i][0])
        rownum = int(tpf.row + pixel_coords[i][1])

        idx = np.where((colrow_unique.T[0] == colnum) & (colrow_unique.T[1] == rownum))[0][0]

        #@|we resample the PRF into the shape of the original tpf
        try:
            resampled_ac = prf_array[idx].locate(pixel_coords[i][0], pixel_coords[i][1],                                                  table['flux'][i], tpf.shape[1:3])
        except:
            #@|this exception occours when there are Gaia sources too far away 
            resampled_ac = np.zeros(tpf.shape[1:3]) 
            
        resampled_list.append(resampled_ac)
        resampled = resampled + resampled_ac


# In[ ]:


#@|-----------------CROWDSAP pixel by pixel--------------------#@|
CROWDSAP_pixel_by_pixel = resampled_list[idx_target] / resampled
#@|------------------------------------------------------------#@|


# In[ ]:


#@|-----------------------------------------#@|
#@|---We select/create the proper apeture---#@|
#@|-----------------------------------------#@|

if aperture == 'pipeline':
    aperture_mask = tpf.pipeline_mask
    if len(np.where(aperture_mask==True)[0]) == 0:
        print(f'The target TIC {tic} does not have a pipeline aperture in sector {sector}. Please modify your config.ini file so that **aperture: threshold_target_flux** or **aperture: threshold_median_flux** in order to create your own aperture (see the documentaion for details on each aperture creation method).')
        sys.exit()
        
if aperture == 'threshold_median_flux':
    aperture_mask = tpf.create_threshold_mask(threshold=threshold_median)
    
if aperture == 'threshold_target_flux':
    aperture_mask = CROWDSAP_pixel_by_pixel > threshold_target  
    
#@|we save te aperture in a .csv file#@|
if save_aper:
    if aperture == 'threshold_target_flux':
        pd.DataFrame(aperture_mask).to_csv(f'output/{target_name}/{target_name}_S{sector}_aperture_{aperture}_{threshold_target}.csv',                                        index = False)
    if aperture == 'threshold_median_flux':
        pd.DataFrame(aperture_mask).to_csv(f'output/{target_name}/{target_name}_S{sector}_aperture_{aperture}_{threshold_median}.csv',                                        index = False)   
        
    #@|open as (e.g) aperture_mask = np.array(pd.read_csv('TIC_282485660_S12_aperture_threshold_median_flux_3.csv'))


# In[ ]:


#@|--------------
#@|GET 'FLFRCSAP'
#@|---------------

#@|'FLFRCSAP' is the flux fraction of the target star inside the photometric aperture, compared to the total flux
#@|emited by the target star.'FLFRCSAP'only depends on the target star itself.

FLFRCSAP = np.sum(resampled_list[idx_target][aperture_mask])
###print(f'The FLFRCSAP of TIC {tic} in Sector {sector} is {FLFRCSAP}.')
#@|total flux of the target star outside the apertue
#np.sum(resampled_list[idx_target][~tpf.pipeline_mask])


#@|--------------
#@|GET 'CROWDSAP'
#@|---------------

#@|'CROWDSAP' is the flux fraction of the target star inside the photometric aperture, compared to the total flux
#@|inside the aperture coming from all the sources.'CROWDSAP' depends on the target stars and all nearby sources.

CROWDSAP = np.sum(resampled_list[idx_target][aperture_mask]) / np.sum(resampled[aperture_mask])
###print(f'The CROWDSAP of TIC {tic} in Sector {sector} is {CROWDSAP}.')


if save_metrics == True:    

    f = open("metrics.dat", "a+") 
    with open("metrics.dat", "r") as f:
        f_read = f.read()
        
    with open("metrics.dat", "a+") as f:

        if 'config_file' in f_read:
            pass
        if'config_file' not in f_read:
            f.write('config_file,CROWDSAP,FLFRCSAP \n')
            
        if config_file in f_read:
            ###print(f'The CROWDSAP and FLFRCSAP metrics were already computed for {config_file}')
            pass
        if config_file not in f_read:
            f.write(f'{config_file},{CROWDSAP},{FLFRCSAP}'+"\n")

        
        


# In[ ]:


#@|Which targets have the main flux contribution to the aperture? 
#@|we compute the CROWDSAPs of all targtes, and select the highest ones
CROWDSAP_list = []
for i in range(len(table)):
    
    CROWDSAP_ac = np.sum(resampled_list[i][aperture_mask]) / np.sum(resampled[aperture_mask])
    CROWDSAP_list.append(CROWDSAP_ac)
    
CROWDSAP_arr = np.array(CROWDSAP_list)
    
    #print(i)


# In[ ]:


#@|we select the targets that contaminate the aperture ('idxs_contam')
#@|that is, 'idxs_contam' contains all the indexes except that of the target star
idxs_contam = np.where(CROWDSAP_arr!=CROWDSAP)[0] 


# In[ ]:


#@|we compute the contamination ratio of the contaminant sources, with respect to the overall contamination.
relative_contam = CROWDSAP_arr[idxs_contam] / np.sum(CROWDSAP_arr[idxs_contam])


# In[ ]:


#@|we sort the contaminant sources, from more to less contaminant
relative_contam_sorted = sorted(relative_contam, reverse = True)
#contaminant_sorted


# In[ ]:


if n_sources == 0:
    raise Exception("Sorry! We are studying flux contamination, so the number of contaminant souces n_sources cannot be 0.") 
if n_sources < 0:
    raise Exception('Sorry! The number of contaminant sources must be a positive number.')
if type(n_sources) != int:
    raise Exception('Sorry! The number of contaminant sources must be a positive integer.')
relative_contam = relative_contam_sorted[:n_sources]
relative_contam_rest = relative_contam_sorted[n_sources:]
relative_contam.extend([np.sum(relative_contam_rest)]) #@|we include the remaining contamination from 'Other' stars


# In[ ]:


#@|we get the indexes of the n_sources more contaminant sources (on 'table')
#@|Note: np.argsot does not have an argument 'reverse', so the highest PDCSAPS are at the end of the array

idxs_crowdsap_sorted = np.argsort(CROWDSAP_arr)
crowdsap_sorted = np.sort(CROWDSAP_arr)
idxs_without_CROWDSAP = np.where(crowdsap_sorted != CROWDSAP)

idxs_crowdsap_sorted = idxs_crowdsap_sorted[idxs_without_CROWDSAP]
crowdsap_sorted  = crowdsap_sorted[idxs_without_CROWDSAP]

#@|############################################################
idxs_selected_contaminant_sources = idxs_crowdsap_sorted[-n_sources:]
#@|###########################################################
#CROWDSAP_arr[idxs_selected_contaminant_sources] / np.sum(CROWDSAP_arr[idxs_contam]) #@|just a check we are doing
#@|things right.

idxs_selected_contaminant_sources = idxs_selected_contaminant_sources[::-1] #@|we reverse it so that the first 
#@|idexes correspond to the most contaminant sources


# In[ ]:


#@|we create the automatic label for the pie chart
#@|criterion. Star#1 is the most contaminant within the aperture, Star#2 the second most contaminant, etc.
bar_labels = []
gaia_names_selected = []
pixel_coords_selected = []
for i,index in enumerate(idxs_selected_contaminant_sources):
    
    bar_labels.append(f'Star {i+1}')      #@|this is the good one
    #@|gaia names
    gaia_names_selected.append(table['Source'][index])
    #@|coordinates
    pixel_coords_selected.append(pixel_coords[index])

bar_labels.append('Other')


# In[ ]:


#@|This is to generate fancy color palettes
#@|from https://stackoverflow.com/questions/876853/generating-color-ranges-in-python
#from colorsys import hsv_to_rgb  (moved up)
def get_hex_color_list(num_colors=5, saturation=0.4, value=1.0):
    hex_colors = []
    hsv_colors = [[float(x / num_colors), saturation, value] for x in range(num_colors)]
    
    for hsv in hsv_colors:
        hsv = [int(x * 255) for x in hsv_to_rgb(*hsv)]
    
        #Formatted as hexadecimal string using the ':02x' format specifier
        hex_colors.append(f"#{hsv[0]:02x}{hsv[1]:02x}{hsv[2]:02x}")
    
    return hex_colors


# In[ ]:


#@|################
#@|+++Pie chart++++
#@|################

print(f'Generating the pie chart plot of {target_name} (Sector {sector}) ...')

# make figure and assign axis objects
#fig, (ax1, ax2) = plt.subplots(1, 2,  width_ratios=[3, 1], figsize=(9, 5))
fig, (ax1, ax2) = plt.subplots(1, 2, gridspec_kw={'width_ratios': [1,6]}, figsize=(6.93, 5.5))
fig.subplots_adjust(top = 1.41, bottom = 0, right = 1, left = 0, 
        hspace = 0, wspace = 0)


#@|------PIE CHART-----#@|

# pie chart parameters
pie_ratios = [1-CROWDSAP, CROWDSAP]
#pie_labels = ['Nearby', 'Target']     
pie_labels = ['', '           ']  #this
pie_colors = ['lightgrey','#25BDB0']
explode = [0.5, 0]
# rotate so that first wedge is split by the x-axis
angle = -180 * pie_ratios[0]
#wedges, *_ = ax1.pie(pie_ratios, autopct='%1.1f%%', textprops = {'size': '21'}, startangle=angle,
                     #labels=pie_labels, explode=explode, colors = pie_colors, radius = 6)
wedges, *_ = ax1.pie(pie_ratios, startangle=angle, explode=explode,                      colors = pie_colors, radius = 6, labels = pie_labels)

#@|-----------------------Pie chart percentages inside the pie------------------------
#@|+++Nearby stars+++
if pie_ratios[0] >= 0.9995:
    ax1.text(1.60, -0.2, f'Nearby ({"{:1.2f}".format(pie_ratios[0]*100)}%)', fontsize = 17)
if 0.000995 < pie_ratios[0] < 0.0995:
    ax1.text(1.20, -0.2, f'Nearby ({"{:1.1f}".format(pie_ratios[0]*100)}%)', fontsize = 17)
if pie_ratios[0] < 0.000995:
    ax1.text(0.80, -0.2, f'Nearby ({"{:1.2f}".format(pie_ratios[0]*100)}%)', fontsize = 17)    
if 0.0995 < pie_ratios[0] < 0.9995:
    ax1.text(1.45, -0.2, f'Nearby ({"{:1.1f}".format(pie_ratios[0]*100)}%)', fontsize = 17)
#@|+++Target star+++
if pie_ratios[1] >= 0.9995:
    ax1.text(-5.2, -0.2, f'Target ({"{:1.2f}".format(pie_ratios[1]*100)}%)', fontsize = 17)
if 0.000995 < pie_ratios[1] < 0.0995:
    ax1.text(-4.9, -0.2, f'Target ({"{:1.1f}".format(pie_ratios[1]*100)}%)', fontsize = 17)
if pie_ratios[1] < 0.000995:
    ax1.text(-5.125, -0.2, f'Target ({"{:1.2f}".format(pie_ratios[1]*100)}%)', fontsize = 17)    
if 0.0995 < pie_ratios[1] < 0.9995:
    ax1.text(-5.125, -0.2, f'Target ({"{:1.1f}".format(pie_ratios[1]*100)}%)', fontsize = 17)
#@|------------------------------------------------------------------------------------


#@|------BAR CHART-----#@|

# bar chart parameters
bar_ratios = relative_contam
bar_labels = bar_labels
bar_colors = get_hex_color_list(num_colors = n_sources+1, saturation = 0.5, value = 1.0)
bottom = 1
width = 0.25

# Adding from the top matches the legend.
labels_ordered = []
for j, (height, label) in enumerate(sorted([*zip(bar_ratios, bar_labels)])):
    bottom -= height
    #bc = ax2.bar(0, height, width, bottom=bottom, color=bar_colors[j], label=label,
                 #alpha=0.1 + 0.25 * j) #@|old. to remove
    bc = ax2.bar(0, height, width, bottom=bottom, color=bar_colors[j], label=label) #@|updated!
    ax2.bar_label(bc, labels=[f"{height:0.0%}"], label_type='center', fontsize = 17)
    
    labels_ordered.append(label) #this is for the tpf plot
labels_ordered.reverse() #most contaminant sources first
idx_other = np.where(np.array(labels_ordered) == 'Other')[0][0] #idxs of 'Other' to skip the color code in the next plot

#ax1.set_title('Contaminant flux', fontsize = 14)
#ax2.legend(bbox_to_anchor=(0.462,0.865), loc="upper left",  bbox_transform=fig.transFigure, fontsize = 15)
ax2.legend(bbox_to_anchor=(0.420,1.422), loc="upper left",            bbox_transform=fig.transFigure, fontsize = 18, framealpha=0.8)
ax2.axis('off')
#ax2.set_xlim(- 2.5 * width, 2.5 * width)
ax2.set_xlim(- 1 , 0.3)
ax2.set_ylim(-0.05, 1.0)

# use ConnectionPatch to draw lines between the two plots
theta1, theta2 = wedges[0].theta1, wedges[0].theta2
center, r = wedges[0].center, wedges[0].r
bar_height = sum(bar_ratios)

# draw top connecting line
x = r * np.cos(np.pi / 180 * theta2) + center[0]
y = r * np.sin(np.pi / 180 * theta2) + center[1]
con = ConnectionPatch(xyA=(-width / 2, bar_height), coordsA=ax2.transData,
                      xyB=(x, y), coordsB=ax1.transData)
con.set_color('k')
con.set_linewidth(2)
ax2.add_artist(con)

# draw bottom connecting line
x = r * np.cos(np.pi / 180 * theta1) + center[0]
y = r * np.sin(np.pi / 180 * theta1) + center[1]
con = ConnectionPatch(xyA=(-width / 2, 0), coordsA=ax2.transData,
                      xyB=(x, y), coordsB=ax1.transData)
con.set_color('k')
ax2.add_artist(con)
con.set_linewidth(2)


#|----------------------
#@|save the pie chart#@|
#|----------------------

if img_fmt == 'pdf' or img_fmt == 'pdfpng':
    plt.savefig(f'output/{target_name}/{target_name}_S{sector}_piechart.pdf',                bbox_inches = 'tight', pad_inches = 0)
if img_fmt == 'png'or img_fmt == 'pdfpng':
    plt.savefig(f'output/{target_name}/{target_name}_S{sector}_piechart.png',                bbox_inches = 'tight', dpi = 400, pad_inches = 0)
    
#print('')
print('\033[1m' + f'Your pie chart {target_name}_S{sector}_piechart.pdf/png has been successfully generated and saved'+'\033[0m')


# In[ ]:


print(f'Generating the heatmap plot of {target_name} (Sector {sector}) ...')
fig = plt.figure(figsize=(6.93, 5.5))
_, ny, nx = np.shape(tpf)
gs = gridspec.GridSpec(1,3, height_ratios=[1], width_ratios=[1,0.05,0.01])
gs.update(left=0.05, right=0.95, bottom=0.12, top=0.95, wspace=0.01, hspace=0.03)
ax1 = plt.subplot(gs[0,0])
maskcolor = 'red'

norm = ImageNormalize(stretch=stretching.LogStretch(), vmin = 0.1, vmax = 99) 
#norm = ImageNormalize(vmin = 0.1, vmax = 99)                                  

splot = plt.imshow(CROWDSAP_pixel_by_pixel*100,                   zorder=0,alpha =1, 
          extent=[tpf.column-0.5,tpf.column+nx-0.5,tpf.row+ny-0.5,tpf.row-0.5], norm = norm, cmap = 'viridis')


for i in range(aperture_mask.shape[0]):
    for j in range(aperture_mask.shape[1]):
        if aperture_mask[i, j]:
            ax1.add_patch(patches.Rectangle((j+tpf.column-0.5, i+tpf.row-0.5),
                                            1, 1, color=maskcolor, fill=True,alpha=0.2))
            ax1.add_patch(patches.Rectangle((j+tpf.column-0.5, i+tpf.row-0.5), 
                                            1, 1, color=maskcolor, fill=False,alpha=1,lw=2))
                    
                
if plot_percentages:                
                    
    for i in range(tpf.shape[1]):
        for j in range(tpf.shape[2]):
            
            if np.round(CROWDSAP_pixel_by_pixel[i, j] * 100, 1) == 0.0:
                continue

            #@|trick to avoid 100.0 values (put instead 100)
            if np.round(CROWDSAP_pixel_by_pixel[i, j] * 100, 1) == 100.0:
                text = ax1.text(j+tpf.column, i+tpf.row, str(100),
                           ha="center", va="center", color="k", zorder = 1000, fontsize = 10.5) 

            else: 
                #text = ax1.text(j+tpf.column, i+tpf.row, np.round(CROWDSAP_pixel_by_pixel[i, j] * 100, 1),
                               #ha="center", va="center", color="k", zorder = 1000, fontsize = 10.5)
                
                text = ax1.text(j+tpf.column, i+tpf.row, np.round(CROWDSAP_pixel_by_pixel[i, j] * 100, 1),
                               ha="center", va="center", color="k", zorder = 1000, fontsize = 10.5)
            

#@|#####
plt.xlabel('Pixel Column Number', fontsize=14, zorder=200)
plt.ylabel('Pixel Row Number', fontsize=14, zorder=200)
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)


#@|##################################################@|
#@|---Include the star locations within the plot----#@|
#@|##################################################@|
#@|Note. The circle sizes are scaled to the stellar fluxes (table['flux']). In particular, the 's' parameter
#@|is proportional to table['flux'], so that the fluxes are proportinal to the total area.

        
#@|our target star
if plot_target:
    plt.scatter(pixel_coords[idx_target][0]+tpf.column, pixel_coords[idx_target][1]+tpf.row,                      s = (scale_factor*table['flux'][idx_target]), c = '#25BDB0', ec = 'k', lw = 1.5,                 alpha = 0.8, zorder = 99, label = target_name)
        

#@|the N most contaminant sources
if plot_main_contaminants:
    heat_colors =  get_hex_color_list(num_colors = n_sources+1, saturation = 0.5, value = 1.0)    
    heat_colors = heat_colors[::-1]
    #heat_colors.pop(0)
    heat_colors.pop(idx_other)
    #for i,index in enumerate(idxs_selected_contaminant_sources):
    for i,index in reversed(list(enumerate(idxs_selected_contaminant_sources))):
        plt.scatter(pixel_coords_selected[i][0]+tpf.column, pixel_coords_selected[i][1]+tpf.row,                    s = (scale_factor*table['flux'][index]), c = heat_colors[i], alpha = 0.8, ec = 'k',                    zorder = 100000, label = bar_labels[i])      

#@|all the remaining Gaia DR3 sources
if plot_all_gaia:
    for i in range(len(table)):
        plt.scatter(pixel_coords[i][0]+tpf.column, pixel_coords[i][1]+tpf.row,                    s = (scale_factor*table['flux'][i]), c = 'lightgrey', ec = 'w', alpha = 0.3, zorder = 98)
        
    
plt.ylim(tpf.row+ny-0.5, tpf.row-0.5)
#plt.ylim(tpf.row+ny-0.5-1, tpf.row-0.5)       #Arrange for the 13x11 TPF, larger than the PRF arrays
plt.xlim(tpf.column-0.5, tpf.column+nx-0.5)


#@|to fix marker sizes for the legend
#@|see https://stackoverflow.com/questions/24706125/setting-a-fixed-size-for-points-in-legend
legend_marker_size = 60
def updatescatter(handle, orig):
    handle.update_from(orig)
    handle.set_sizes([legend_marker_size])
plt.legend(loc = loc_legend, handler_map={PathCollection : HandlerPathCollection(update_func=updatescatter)},           framealpha=0.8, fontsize = 12).set_zorder(10000)

#@|--------
#@|COLORBAR
#@|--------
cbax = plt.subplot(gs[0,1]) # Place it where it should be.
pos1 = cbax.get_position() # get the original position
pos2 = [pos1.x0 - 0.075, pos1.y0, pos1.width, pos1.height]
cbax.set_position(pos2) # set a new position

#cbar_ticks = np.array([10, 20, 40, 60, 80])
#cb = Colorbar(ax = cbax, mappable = splot, orientation = 'vertical',
              #ticklocation = 'right', ticks =  cbar_ticks)

cb = Colorbar(ax = cbax, mappable = splot, orientation = 'vertical',
              ticklocation = 'right')

if plot_target_name:
    cb.set_label(f'Flux ratio from {target_name} (%)', labelpad=10, fontsize=14)
else:
    cb.set_label('Flux ratio from the target star (%)', labelpad=10, fontsize=14)

#@|---------------------------
#@|----save the heatmap-----#@|
#@|---------------------------


if img_fmt == 'pdf' or img_fmt == 'pdfpng':
    plt.savefig(f'output/{target_name}/{target_name}_S{sector}_heatmap.pdf',                 bbox_inches = 'tight', pad_inches = 0)
if img_fmt == 'png' or img_fmt == 'pdfpng':
    plt.savefig(f'output/{target_name}/{target_name}_S{sector}_heatmap.png',                 bbox_inches = 'tight', pad_inches = 0, dpi = 400)
    
#print('')
print('\033[1m' + f'Your heatmap {target_name}_S{sector}_heatmap.pdf/png has been successfully generated and saved'+'\033[0m')  


# In[ ]:





# In[ ]:


#@|++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#@|save the list of the 'n_sources' most contaminant sources
#@|++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


# In[ ]:


#@|we obtain the TIC IDs (tic_names_selected) of the n_sources from the Gaia DR2 IDS
#@|----This functionality is deprecated due to the limited number of TICs found-----

#tic_names_selected = []
#for i in range(n_sources):
    
    #try:
        #cat = Catalogs.query_object(f'Gaia_DR3_{gaia_names_selected[i]}', catalog="TIC")
        #idx_cat = np.where(cat['GAIA'].value.data == str(gaia_names_selected[i]))[0][0]
        #tic_names_selected.append(cat['ID'][idx_cat])
        
    #except:
        ##print(f'No TIC number found for Gaia_DR3_{gaia_names_selected[i]} ')
        #tic_names_selected.append('No TIC ID found')
        #continue


# In[ ]:





# In[ ]:


f = open(f'output/{target_name}/{target_name}_S{sector}_contaminant_sources.dat', 'w+')

header = ['Star,Gaia_ID,total_cont(%),rel_cont(%) \n']
f.writelines(header)

for i in range(n_sources):
    
    crowdsap_sorted_ac = crowdsap_sorted[-n_sources:][::-1][i] * 100 #(in %)
    total_cont_ac = relative_contam_sorted[:n_sources][i] * 100 #(in %)
    
    L = [str(i+1)+','+str(gaia_names_selected[i])+','+str(round(crowdsap_sorted_ac, 4))         +','+str(round(total_cont_ac, 4))+'\n']
    f.writelines(L)
    
#####################---rep---###############   

f = open(f'output/{target_name}/{target_name}_S{sector}_contaminant_sources.dat', 'w+')

header = ['Star,Gaia_ID,total_cont(%),rel_cont(%) \n']
f.writelines(header)

for i in range(n_sources):
    
    crowdsap_sorted_ac = crowdsap_sorted[-n_sources:][::-1][i] * 100 #(in %)
    total_cont_ac = relative_contam_sorted[:n_sources][i] * 100 #(in %)
    
    L = [str(i+1)+','+str(gaia_names_selected[i])+','+str(round(crowdsap_sorted_ac, 4))         +','+str(round(total_cont_ac, 4))+'\n']
    f.writelines(L)
    

print('\033[1m' + f'The list of the {n_sources} most contaminant sources has been saved in {target_name}_S{sector}_contaminant_sources.dat'+'\033[0m')


# In[ ]:


#@|----------------------
#@|TRANSIT DEPTH ANALYSIS
#@|----------------------


# In[ ]:


if transit_depth_analysis:
    
    if td_unit == 'ppm':
        td = td / 1e6 #transit depth (0-1)
    if td_unit == 'ppt':
        td = td / 1000 #transit depth (0-1)
    if td_unit == 'per':
        td = td / 100 #transit depth (0-1)
        
    #@|transit depth uncorrected (e.g. SPOC td obtained from PDCSAP photometry)
    if dilution_corr == True:
        td = td * CROWDSAP
        
    f = open(f'output/{target_name}/{target_name}_S{sector}_undituled_transit_depths.dat', 'w+')

    #if td_unit == 'ppm':
        #header = ['Gaia_ID,TIC_ID,transit_depth(ppm) \n']
    #if td_unit == 'ppt':
        #header = ['Gaia_ID,TIC_ID,transit_depth(ppt) \n']
    #if td_unit == 'per':
        #header = ['Gaia_ID,TIC_ID,transit_depth(%) \n']
    #if td_unit == 'frac':
        #header = ['Gaia_ID,TIC_ID,transit_depth \n']
        
    header = ['Star,Gaia_ID,transit_depth(%) \n']

    f.writelines(header)

    for i in range(n_sources):

        #if td_unit == 'ppm':
            #td_undiluted = td / crowdsap_sorted[-n_sources:][::-1][i] * 1e6

        #if td_unit == 'ppt':
            #td_undiluted = td / crowdsap_sorted[-n_sources:][::-1][i] * 1000

        #if td_unit == 'per':
            #td_undiluted = td / crowdsap_sorted[-n_sources:][::-1][i] * 100

        #if td_unit == 'frac':
            #td_undiluted = td / crowdsap_sorted[-n_sources:][::-1][i]
            
        td_undiluted = td / crowdsap_sorted[-n_sources:][::-1][i]

        L = [str(i+1)+','+str(gaia_names_selected[i])+','+str(td_undiluted*100)+'\n']
        f.writelines(L)
    
    
    #--------------------writting x2---------------------------
    
    f = open(f'output/{target_name}/{target_name}_S{sector}_undituled_transit_depths.dat', 'w+')

    #if td_unit == 'ppm':
        #header = ['Gaia_ID,TIC_ID,transit_depth(ppm) \n']
    #if td_unit == 'ppt':
        #header = ['Gaia_ID,TIC_ID,transit_depth(ppt) \n']
    #if td_unit == 'per':
        #header = ['Gaia_ID,TIC_ID,transit_depth(%) \n']
    #if td_unit == 'frac':
        #header = ['Gaia_ID,TIC_ID,transit_depth \n']
        
    header = ['Star,Gaia_ID,transit_depth(%) \n']

    f.writelines(header)

    for i in range(n_sources):

        #if td_unit == 'ppm':
            #td_undiluted = td / crowdsap_sorted[-n_sources:][::-1][i] * 1e6

        #if td_unit == 'ppt':
            #td_undiluted = td / crowdsap_sorted[-n_sources:][::-1][i] * 1000

        #if td_unit == 'per':
            #td_undiluted = td / crowdsap_sorted[-n_sources:][::-1][i] * 100

        #if td_unit == 'frac':
            #td_undiluted = td / crowdsap_sorted[-n_sources:][::-1][i]

        td_undiluted = td / crowdsap_sorted[-n_sources:][::-1][i]
        
        L = [str(i+1)+','+str(gaia_names_selected[i])+','+str(td_undiluted*100)+'\n']
        f.writelines(L)
        
        
    print('\033[1m' + f'The transit depth analysis has been saved in {target_name}_S{sector}_undituled_transit_depths.dat ')
    print('\033[0m')
    print('     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++')
    print('      Thank you for using TESS-cont :-D We hope to see you again soon!')
    print('     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++')
    print('')
        
        
else:
    
    print('')
    print('     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++')
    print('      Thank you for using TESS-cont :-D We hope to see you again soon!')
    print('     ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++')
    print('')


# In[ ]:





# In[ ]:




