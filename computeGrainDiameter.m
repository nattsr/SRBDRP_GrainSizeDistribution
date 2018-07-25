function [grainDiameter] = computeGrainDiameter(grainRadius, measureOption)
%computeGrainDiameter compute grain diameter
%   Input Arguments
%   - grainRadius   : a (nGrain*6) or (nGrain*4) double matrix,
%                     radius of each grain [r1 r2 r3 r4 r5 r6] or
%                     [r1 r2 r3 r4] in Micron
%   - measureOption : a string, compute the distribution using
%                     "max" grain size
%                     "mean" grain size
%
%   Output Arguments
%   - grainDiameter : a (nGrain*1) double matrix,
%                     a max diameter of each grain

%   Revision 1: May 2018 Nattavadee Srisutthiyakorn



%%
if (~exist('measureOption', 'var'))
    measureOption = "max";
end



[~, nRadius]= size(grainRadius);
% Add to get the grain diameter instead of the radius.
if nRadius == 6
    grainDiameter(:,1) = grainRadius(:,1) + grainRadius(:,2);
    grainDiameter(:,2) = grainRadius(:,3) + grainRadius(:,4);
    grainDiameter(:,3) = grainRadius(:,5) + grainRadius(:,6);
elseif nRadius == 4
    grainDiameter(:,1) = grainRadius(:,1) + grainRadius(:,2);
    grainDiameter(:,2) = grainRadius(:,3) + grainRadius(:,4);
end

% Use the find the average grain size or max grain size for sieving
if measureOption == "mean"
    grainDiameter  = mean(grainDiameter,2);
elseif measureOption == "max"
    grainDiameter  = max(grainDiameter,[],2);
end


end

