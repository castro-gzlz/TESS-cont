#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 29 11:04:56 2021

Package to display how a point source should appear at any location on the 
TESS detector. Two options: TESS_PRF, referencing the actual Pixel Response
Function models on MAST (or downloaded locally). Gaussian_PRF will generate 
an Gaussian PRF model. If you want to know how a star will appear in a given 
Target Pixel File (TPF)), initialize the PRF object for the center of the PRF 
(in terms of row/column), then call the `resample` function to position the 
star at a specific location within the PRF (row/column relative to TPF corner).

@author: keatonb
"""
import numpy as np
from bs4 import BeautifulSoup
import requests
from astropy.io import fits
from scipy.interpolate import RectBivariateSpline
import os
from glob import glob
    
class TESS_PRF:
    """TESS Pixel Response Function object
    
    """
    def __init__(self,cam,ccd,sector,colnum,rownum, localdatadir = None):
        """Get TESS PRF for detector location, sector.
        
        Downloads relevant PRF files from the MAST archive by default.
        
        ***To use pre-downloaded local files, give directory containing
        subdirectories of format "cam#_ccd#/" as localdatadir, appropriate
        for sector of interest (separate for Sectors 1-3, 4+)
        
        inputs:
         - cam (int): TESS camera number
         - ccd (int): TESS ccd number
         - sector (int): TESS sector number
         - colnum (float): column number near target
         - rownum (float): row number near target
         
        """
        self.cam,self.ccd,self.sector,self.colnum,self.rownum = cam,ccd,sector,colnum,rownum
        self.prfnsamp = 9 #samples/pixel for TESS PRFs
        
        #Get PRF file info
        subdir = f'cam{int(cam)}_ccd{int(ccd)}/'
        filelist = None #local and online options
        
        if localdatadir is None:
            #Different PRFs for Sectors 1-3, and 4+
            #https://heasarc.gsfc.nasa.gov/docs/tess/observing-technical.html#point-spread-function
            if sector < 4:
                url = 'https://archive.stsci.edu/missions/tess/models/prf_fitsfiles/start_s0001/' + subdir
            else:
                url = 'https://archive.stsci.edu/missions/tess/models/prf_fitsfiles/start_s0004/' + subdir
            ext = 'fits'

            #Get list of available PRF files
            #https://stackoverflow.com/a/34718858
            def listFD(url, ext=''):
                page = requests.get(url).text
                soup = BeautifulSoup(page, 'html.parser')
                return [url + '/' + node.get('href') for node in soup.find_all('a') if node.get('href').endswith(ext)]

            filelist = [file for file in listFD(url, ext)]
        else:
            filelist = glob(os.path.join(localdatadir, subdir) + '*.fits')
        #One directory on MAST has some errant files with `phot` in filename
        filelist = [file for file in filelist if 'phot' not in file]
        
        cols = np.array([int(file[-9:-5]) for file in filelist])
        rows = np.array([int(file[-17:-13]) for file in filelist])

        #Bilinear interpolation between four surrounding PRFs
        LL = np.where((rows < rownum) & (cols < colnum))[0] #lower left
        LR = np.where((rows > rownum) & (cols < colnum))[0] #lower right
        UL = np.where((rows < rownum) & (cols > colnum))[0] #upper left
        UR = np.where((rows > rownum) & (cols > colnum))[0] #uppper right
        dist = np.sqrt((rows-rownum)**2. + (cols-colnum)**2.)
        surroundinginds = [subset[np.argmin(dist[subset])] for subset in [LL,LR,UL,UR]]
        #Following https://stackoverflow.com/a/8662355
        points = []
        for ind in surroundinginds:
            hdulist = fits.open(filelist[ind])
            prf = hdulist[0].data
            points.append((cols[ind],rows[ind],prf))
            hdulist.close()
        points = sorted(points)
        (x1, y1, q11), (_x1, y2, q12), (x2, _y1, q21), (_x2, _y2, q22) = points

        self.prf = (q11 * (x2 - colnum) * (y2 - rownum) +
                    q21 * (colnum - x1) * (y2 -  rownum) +
                    q12 * (x2 - colnum) * ( rownum - y1) +
                    q22 * (colnum - x1) * ( rownum - y1)
                    ) / ((x2 - x1) * (y2 - y1) + 0.0)
        
        ## Need to reshape PRF file for interpolation
        ## Add models just beyond pixel edges
        
        ##Size: 11x11x13x13 
        #indices: subrow index (from bottom), subcol index (from left),
        #TPF row index (from bottom), TPF col index (from left),
        reshaped = np.zeros((11,11,13,13)) 
        
        #Un-interleve PRF samples
        for i in range(9): #col number
            for j in range(9): #row number
                reshaped[i+1,j+1,:,:] = self.prf[8-i::9, 8-j::9] #col, row, cols, rows
        
        #Add columns just beyond pixel edges
        for j in range(1,10): #loop over rows
            reshaped[j,0,:,:] = np.append(reshaped[j,-2,:,1:],np.zeros((13,1)),axis=1)
            reshaped[j,-1,:,:] = np.append(np.zeros((13,1)),reshaped[j,1,:,:-1],axis=1)
        
        #Add rows just beyond pixel edges
        for i in range(0,11):
            reshaped[0,i,:,:] = np.append(reshaped[-2,i,1:,:],np.zeros((1,13)),axis=0)
            reshaped[-1,i,:,:] = np.append(np.zeros((1,13)),reshaped[1,i,:-1,:],axis=0)
        
        #Store for later use
        self.reshaped = reshaped
        
    def locate(self, sourcecol, sourcerow, fac, stampsize=(13,13)):
        """Interpolate TESS PRF at location within "interleaved" TPF
        
        sourcecol (float): col position of star (relative to TPF)
        sourcerow (float): row position of star (relative to TPF)
        stampsize (int,int): (height,width) of TPF
        
        Note: pixel positions follow the convention that integers refer to the 
        pixel center.
        """
        #Break into integer and fractional pixel
        #adding 0.5 to conform to convention
        colint = np.floor(sourcecol+0.5)
        colfract = (sourcecol+0.5) % 1
        rowint = np.floor(sourcerow+0.5)
        rowfract = (sourcerow+0.5) % 1
        
        #Sub-pixel sample locations (in each dirextion, w/ border added)
        pixelsamples = np.arange(-1/18,19.1/18,1/9)
        
        #Find four surrounding subpixel PRF models
        colbelow = np.max(np.where(pixelsamples < colfract)[0])
        colabove = np.min(np.where(pixelsamples >= colfract)[0])
        rowbelow = np.max(np.where(pixelsamples < rowfract)[0])
        rowabove = np.min(np.where(pixelsamples >= rowfract)[0])
        
        #interpolate
        points = []
        LL = self.reshaped[rowbelow,colbelow,:,:]
        points.append((pixelsamples[colbelow],pixelsamples[rowbelow],LL))
        UL = self.reshaped[rowabove,colbelow,:,:]
        points.append((pixelsamples[colbelow],pixelsamples[rowabove],UL))
        LR = self.reshaped[rowbelow,colabove,:,:]
        points.append((pixelsamples[colabove],pixelsamples[rowbelow],LR))
        UR = self.reshaped[rowabove,colabove,:,:]
        points.append((pixelsamples[colabove],pixelsamples[rowabove],UR))
    
        points = sorted(points)
        (x1, y1, q11), (_x1, y2, q12), (x2, _y1, q21), (_x2, _y2, q22) = points
        
        subsampled = (q11 * (x2 - colfract) * (y2 - rowfract) +
                      q21 * (colfract - x1) * (y2 -  rowfract) +
                      q12 * (x2 - colfract) * ( rowfract - y1) +
                      q22 * (colfract - x1) * ( rowfract - y1)
                      ) / ((x2 - x1) * (y2 - y1) + 0.0)
        #re-normalize to 1
        subsampled /= np.sum(subsampled)
        subsampled *= fac
        
        #Now must place at correct location in TPF
        tpfmodel = np.zeros(stampsize)
        
        #PRF models are 13x13 pixels
        #center of PRF is pixel (6,6)
        midprf = 6
        
        #That pixel should be (colint,rowint) in TPF
        tpfmodel[int(np.max([0,rowint-midprf])):int(np.min([stampsize[0],rowint+midprf+1])),
                 int(np.max([0,colint-midprf])):int(np.min([stampsize[1],colint+midprf+1])),] = subsampled[
            int(np.max([0,midprf-rowint])):int(np.min([2*midprf+1,midprf-rowint+stampsize[0]])),
            int(np.max([0,midprf-colint])):int(np.min([2*midprf+1,midprf-colint+stampsize[1]])),
            ]
        
        return tpfmodel
        
class Gaussian_PRF:
    """Gaussian Pixel Response Function object
    
    """
    def __init__(self,sigma=1.,npixels=15,nsamp=9):
        """Generate a synthetic, Gaussian Pixel Response Function
        
        inputs:
         - sigma (float, default 1.0): sigma of the Gaussian, in TESS pixels
         - npixels (int, default 15): width of sampled PRF in TESS pixels
         - nsamp (int, default 9): number of PRF samples per TESS pixel
        """
        self.sigma = sigma
        #force integers
        self.npixels = int(npixels)
        self.prfnsamp = int(nsamp)
        
        #center coords
        x0,y0 = (self.npixels/2.,self.npixels/2.) #in pixel units
        row,col = np.meshgrid(np.arange(self.npixels*self.prfnsamp)+0.5, 
                              np.arange(self.npixels*self.prfnsamp)+0.5)
        row /= self.prfnsamp
        col /= self.prfnsamp
        #2D Gaussian, normalized to integrate to 1
        self.prf = np.exp(-(((row-x0)/(self.sigma))**2+((col-y0)/(self.sigma))**2)/2)
        self.prf /= 2*np.pi*self.sigma**2 #normalize
        
    def locate(self, sourcecol, sourcerow, fac, stampsize=(13,13), 
             supersamplefactor = 9):
        """Sample Gaussian PRF at location within TPF
        
        sourcecol (float): col position of star (relative to TPF)
        sourcerow (float): row position of star (relative to TPF)
        stampsize (int,int): (height,width) of TPF
        supersamplefactor (int): interpolate supersampled before downsampling (default 9)
        """
        
        #Center PRF on origin
        prfcol = np.arange(-self.prf.shape[1]/2.+0.5,self.prf.shape[1]/2.+0.5)
        prfrow = np.arange(-self.prf.shape[0]/2.+0.5,self.prf.shape[0]/2.+0.5)

        #Convert to relative locations in TPF
        relprfcol = prfcol/self.prfnsamp + sourcecol
        relprfrow = prfrow/self.prfnsamp + sourcerow

        #Supersample the stamp, then downsample
        supercol = (np.arange(stampsize[1]*supersamplefactor) + 0.5) / supersamplefactor
        superrow = (np.arange(stampsize[0]*supersamplefactor) + 0.5) / supersamplefactor
        
        #Interpolate PRF values onto supersampled stamp pixels
        interppix = RectBivariateSpline(relprfrow,relprfcol,self.prf)
        interped = interppix(superrow,supercol) #Interpolate
        
        #re-normalize to 1 (acg)
        interped  /= np.sum(interped)
        interped  *= fac
        

        #Replace values outside PRF with zeros
        #@|interped[:,np.where(supercol < np.min(relprfcol))] = 0
        #@|interped[:,np.where(supercol > np.max(relprfcol))] = 0
        #@|interped[np.where(superrow < np.min(relprfrow)),:] = 0
        #@|interped[np.where(superrow > np.max(relprfrow)),:] = 0

        #Downsample to TESS TPF resolution
        output = np.zeros([int(s/supersamplefactor) for s in interped.shape])
        for i in range(supersamplefactor):
            for j in range(supersamplefactor):
                output += interped[i::supersamplefactor,j::supersamplefactor]
        output /= supersamplefactor**2.
        return output
