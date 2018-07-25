function [histPC, binCenter, statGSD, idxSelectedGrain] ...
    = computeHistPC( grainProp, grainCentroid, gridPC, binEdge )
%computeHistPC compute frequency-based histogram of the different properties of GSD 
%   Input Arguments
%   - grainProp     : a (nGrain*1) double vector, grain properties such as
%                     maximum diameter, nContact, and etc. 
%   - grainCentroid : a (nGrain*2) or (nGrain*3) integer matrix, xy or xyz 
%                     location of the grain in voxel
%   - gridPC        : (optional) a vector, for creating a grid for point
%                     count
%   - binEdge       : (optional) a vector, histogram bin edge
%
%   Output Arguments
%   - histFB        : a (nBin*1) vector, probability distribution function (PDF) 
%   - binCenter     : a (nBin*1) vector, bin center in mm as specified in the code
%   - statGSD       : a struct, containing grain size distribution
%                     statistics
%   - idxSelectedGrain : a (nSelectedGrains*1), the index of the grain that
%                        is closest to the grid.
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
if (~exist('gridPC', 'var'))
    maxPixel    = max(max(grainCentroid));
    maxNumGrain    = length(grainProp);
    nSampling   = round((maxNumGrain./10).^(1/3)); % The total of nSampling*nSampling samples for 2-D and
                      % and nSampling*nSampling*nSampling samples for 3-D will be taken from point count
    gridPC      = linspace(1, maxPixel, nSampling);
    gridPC      = round(gridPC);
end 
   
% Create the grid
[x, y, z]       = meshgrid(gridPC);
gridCentroid    = [x(:) y(:) z(:)];
idxSelectedGrain = dsearchn(grainCentroid, gridCentroid);

% Select only the index
grainProp       = grainProp(idxSelectedGrain);

% Sort data into histogram bin
nEdge           = length(binEdge);
nBin            = nEdge - 1;
histPC          = zeros(nEdge - 1, 1);
binCenter(1)    = binEdge(1);

for iEdge = 1:nEdge - 1
    [index{iEdge}] = find(and(grainProp <= binEdge(iEdge + 1),...
                             grainProp  >  binEdge(iEdge)));
    binCenter(iEdge)   = (binEdge(iEdge + 1) + binEdge(iEdge))./2; 
end

for iBin = 1:nBin
    histPC(iBin)  = length(index{iBin});
end

% Normalized to get the PDF
nGrainCheck = sum(histPC);
histPC = histPC./nGrainCheck;

% Export the statistics
statGSD.mean = mean(grainProp);
statGSD.std  = std(grainProp);

% Find the cumulative density function 
% histCDF = cumsum((histFB))./sum(histFB);



