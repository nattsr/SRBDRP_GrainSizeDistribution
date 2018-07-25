function [histVB, binCenter, statGSD] ...
    = computeHistVB( grainProp, grainVolume, binEdge )
%computeHistVB computes volume-based histogram of the different properties of GSD 
%   Input Arguments
%   - grainProp     : a (nGrain*1) double vector, grain properties such as
%                     maximum diameter, nContact, and etc. 
%   - grainVolume   : a (nGrain*1) double vector, volume of each grain or
%                     other weight
%                     *** Beware that the output from computeGSD is the
%                     voxel, please multiply grainVolume*resolution.^3
%   - binEdge       : (optional) a vector, histogram bin edge
%
%   Output Arguments
%   - histVB        : a (nBin*1) vector, probability distribution function (PDF) 
%   - binCenter     : a (nBin*1) vector, bin center in mm as specified in the code
%   - statGSD       : a struct, containing grain size distribution
%                     statistics
%
%   Notes
%   - Beware that the output from computeGSD is the voxel, 
%     please multiply by the resolution in micron for maximum grain radius and
%     multiply by resolution.^3 for grain volume
%     grainProp       = grainProp.*dx./1000;

%   Revision 1: May 2018 Nattavadee Srisutthiyakorn

%% Program
if (~exist('binEdge', 'var'))
    % Prefixed Histogram bin in mm due to the plot nature in log scale 
    % for grain diameter
    binEdge = fliplr([4.0000 3.3600 2.8300 2.3800 2.0000 1.6800 1.4100 1.1900 1.0000 ...
           0.8500 0.7100 0.6000 0.5000 0.4200 0.3500 0.2970 0.2500 0.2100 ...
           0.1770 0.1490 0.1250 0.1050 0.0880 0.0740 0.0620 0.0530 0.0440 ...
           0.0370 0.0310 0.0260 0.0220 0.0190 0.0160 0.0130 0.0110 0.0093 ...
           0.0078 0.0065 0.0055 0.0046 0.0039 0.0033 0.0028 0.0023 0.0019 ...
           0.0016 0.0014 0.0012 0.0010]);
end

% Sort data into histogram bin
nEdge           = length(binEdge);
nBin            = nEdge - 1;
histVB          = zeros(nEdge - 1, 1);
binCenter(1)    = binEdge(1);

for iEdge = 1:nEdge - 1
    [index{iEdge}] = find(and(grainProp <= binEdge(iEdge + 1),...
                             grainProp  >  binEdge(iEdge)));
    binCenter(iEdge)   = (binEdge(iEdge + 1) + binEdge(iEdge))./2; 
end

% Find the weight/volume/area of each bin in order to plot the percentage
for iBin = 1:nBin
    histVB(iBin)  = sum(grainVolume(index{iBin}))./sum(grainVolume);
end

% Normalized to get the PDF
nGrainCheck = sum(histVB);
histVB = histVB./nGrainCheck;

% Export the statistics
statGSD.mean = sum(grainProp.*grainVolume)./sum(grainVolume);   
statGSD.std  = std(grainProp, grainVolume);

% Find the cumulative density function 
% histCDF = cumsum((histVF))./sum(histVF);



