# SRBDRP_GrainSizeDistribution
Codes to measure the grain size distribution from thin sections or µXCT images. 

<div align="center">
    <img width=500 src="https://github.com/nattsr/SRBDRP_GrainSizeDistribution/blob/master/ReadmeFiles/Figure_MeasurementExample.jpg" 
    alt="Process" title="Grain Size Distribution Measurement Process"</img>
</div>

## Overview
Grain Size Distribution is one of the basic measurements for sediment classification. The conventional methods for grain size distribution include the sieve method, the laser diffraction method, and the point-count method. We aimed to develop a robust computer code that simulates these conventional methods. The code can measure grain size distribution on 2-D and 3-D binary images using a watershed algorithm to extract out individual grains, and using principal component algorithms to find the principal axes. The outputs include grain radius for different principal axes, grain volume, grain surface area, principal axes inclinations and azimuths, and the number of contacts for each grain. The calculated distribution can be volume-based, frequency-based, or grid-based. 

<div align="center">
    <img width=750 src="https://github.com/nattsr/SRBDRP_GrainSizeDistribution/blob/master/ReadmeFiles/Figure_3DDiffMethod.jpg" 
    alt="Process" title="Different approaches to calculate the grain size distribution"</img>
</div>

## Requirements
MATLAB is required to run the program. 

## Getting Started
The main function computeGSD.m takes an (nx,ny) or (nx,ny,nz) uint8 matrix, 2-D or 3-D binary 
image of porespace (1 = grain, 0 = pore). The code outputs are grainCentroid, grainRadius, grainAzimuth, grainInclination, 
grainVolume, nContact, grainSurfaceArea.

After computing grain properties. The distribution can be computed using computeHistFB.m, computeHistPC.m, or computeHistVB.m. The inputs should be matrix of diameter which is the result from computeGrainDiameter.m.

## Published Research Studies using computeGSD.m
The publication is in submission.
