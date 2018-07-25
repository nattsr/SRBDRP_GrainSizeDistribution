function [ grainCentroid, grainRadius, grainAzimuth, grainInclination, ...
           grainVolume, nContact, grainSurfaceArea ] ...
    = computeGSD( image, minThres, bc, qcPlot )
%computeGSD compute grain size distribution 
%   Input Arguments
%   - image      : an (nx*ny) or (nx*ny*nz) uint8 matrix, 2-D or 3-D binary 
%                  image of porespace (1 = grain, 0 = pore)
%   - minThres   : an integer, a threshold to suppress all minima in the 
%                  intensity image whose depth is less than this number
%   - bc         : an integer, boundary condition
%                  0: impose no boundary condition 
%                  1: remove grains that are close to the boundary 
%   - qcPlot     : an integer, plot the QC image
%                  1: show plot
%
%   Output Arguments
%   - grainCentroid    : a (nGrain*2) or (nGrain*3) integer matrix, xy or xyz 
%                        location of the grain in voxel
%   - grainRadius      : a (nGrain*4) or (nGrain*6) double matrix, radius
%                        of each grain in voxel. The vector is based on the 
%                        principal component analysis 
%   - grainAzimuth     : a (nGrain*4) or (nGrain*6) double matrix,
%                        azimuth on each axis of radius based on 
%                        spherical cooridnate (ISO physics).
%   - grainInclination : (For 3D), a (nGrain*6) double matrix,
%                        inclination on each axis of radius based on 
%                        spherical cooridnate (ISO physics).
%   - grainVolume      : a (nGrain*1) vector, the volume of grain in voxel
%   - nContact         : a (nGrain*1) vector, the number of contact of each
%                        grain
%   - grainSurfaceArea : a (nGrain*1) vector, the surface area of grain in voxel

%   Revision 3: Nov 2017 Natt Srisutthiyakorn - Add spherical
%                                               coordinate/qc3D
%   Revision 2: Oct 2017 Natt Srisutthiyakorn - Add grain contact




%% Program
% Default parameters
if nargin < 2
    minThres    = 1;
    bc          = 1;
    qcPlot      = 0;
end

% Determine whether it is 2-D or 3-D image
[imSize(1), imSize(2), imSize(3)]    = size(image);



%% Watershed algorithm
% Create the distance map image (100+ s for 1024^3 voxels)
disp('Step 1: Create distance map - find the distance of each point to the nearest solid')
imageDistGrain          = bwdist(~image);
imageDistGrain          = -imageDistGrain;
imageDistGrain(~image)  = -Inf;

% Impose height threshold (500+ s for 1024^3 voxels) if applicable
if minThres > 1
    imageDistGrain = imhmin(imageDistGrain,minThres);
end

% Watershed algorithm (3500+ s for 1024^3 voxels)
disp('Step 2: Apply watershed algorithm')
imageGrainIdx = watershed(imageDistGrain);

% QC pore (component 0 is boundary and now pore)
imageGrainIdx(~image) = 0;

% Find region properties (100+ s for 1024^3 voxels) 
disp('Step 3: Find center of mass and volume')
stats       = regionprops('table', imageGrainIdx, 'Centroid', 'Area');
componentNo = unique(imageGrainIdx);



% Get grain properties (exclude 0)
grainNo             = componentNo(componentNo ~= 0);
grainNo             = grainNo(grainNo ~= 1);
grainNo             = grainNo(grainNo ~= 2); % Noise in 1-2?
nGrain              = length(grainNo);
allGrainCentroid    = round(stats.Centroid(grainNo,:));
allGrainVolume      = round(stats.Area(grainNo,:));
nTotalSize          = length(imageGrainIdx(:));
if imSize(3) == 1     % 2D-------------------------------------------------
    grainRadius      = zeros(nGrain,4);
    grainAzimuth     = zeros(nGrain,4);
    grainInclination = zeros(nGrain,4);
elseif imSize(3) > 1  % 3D-------------------------------------------------
    grainRadius      = zeros(nGrain,6);
    grainAzimuth     = zeros(nGrain,6);
    grainInclination = zeros(nGrain,6);
end

% QC
%imageGrainIdx(:,:,1)
%[vol, idx] = max(allGrainVolume)

% QC Plot
if qcPlot
    if imSize(3) == 1 % 2D-------------------------------------------------
        figure
        subplot(1,3,1)
        imagesc(image);
        axis equal; xlim([0 imSize(2)]); ylim([0 imSize(1)]);
        colormap(flipud(gray))
        title('Original')

        subplot(1,3,2)
        imagesc(imageGrainIdx)
        title('Watershed')
        axis equal; xlim([0 imSize(2)]); ylim([0 imSize(1)]);

        subplot(1,3,3)
        imagesc(imageGrainIdx)
        title('Measurement')
        axis equal; xlim([0 imSize(2)]); ylim([0 imSize(1)]);  
    elseif imSize(3) > 1  % 3D---------------------------------------------
        figure
        colors = jet(nGrain);
        colors = colors(randperm(length(colors)),:); % Randomize the color 
        
        for iGrain = 1:nGrain
            idxSingleGrain  = (imageGrainIdx(:) == grainNo(iGrain));
            imageGrain      = zeros(nTotalSize,1);
            imageGrain(idxSingleGrain) = 1;
            imageGrain      = reshape(imageGrain,[imSize(1), imSize(2), imSize(3)]);
            
            fv = isosurface(imageGrain,0);
            patch(fv,'EdgeColor','none','facecolor',colors(iGrain,:));
            box on;
            view(45,45);axis equal
            alpha(0.2)
            hold on
        end
    end
end    
    
 

%% Find the surface area of each grain using the isosurface
grainSurfaceArea = zeros(nGrain,1);
tic
for iGrain = 1:nGrain
    if iGrain/100 == round(iGrain/100)
    disp(['Step 4: Measure grain surface area (', num2str(iGrain),...
        '/',num2str(nGrain), ')'])
    end
    idxSingleGrain      = (imageGrainIdx(:) == grainNo(iGrain));
    imageGrain    = zeros(nTotalSize,1);
    imageGrain(idxSingleGrain) = 1;  
    
    if imSize(3) > 1
        imageGrain      = reshape(imageGrain,[imSize(1), imSize(2), imSize(3)]);
        % Extract the surface surface
        fv = isosurface(imageGrain,0);
        % p =  patch(fv,'facecolor','cyan','EdgeColor','none'); - QC plots
        % verts = get(p, 'Vertices');
        % faces = get(p, 'Faces');
        % close all;
        
        % Find the surface area
        vertices = fv.vertices;
        faces = fv.faces;
        a = vertices(faces(:, 2), :) - vertices(faces(:, 1), :);
        b = vertices(faces(:, 3), :) - vertices(faces(:, 1), :);
        c = cross(a, b, 2);
        grainSurfaceArea(iGrain) = 1/2 * sum(sqrt(sum(c.^2, 2)));
    end
    
end
toc

    

%% Find the number of contact
nContact = zeros(nGrain,1);

for iGrain = 1:nGrain
    if iGrain/100 == round(iGrain/100)
    disp(['Step 5: Measure number of contact(', num2str(iGrain),...
        '/',num2str(nGrain), ')'])
    end
    idxSingleGrain      = (imageGrainIdx(:) == grainNo(iGrain));
    imageGrain    = zeros(nTotalSize,1);
    imageGrain(idxSingleGrain) = 1;  
    
    if imSize(3) == 1
        imageGrain    = reshape(imageGrain,[imSize(1), imSize(2)]);
        % Dilate the grain first to get the outer boundary.
        imageGrain = bwmorph(imageGrain,'dilate');
        imageGrain = bwmorph(imageGrain,'dilate');
        % Find boundary of a grain
        grainBound  = bwboundaries(imageGrain,'noholes');
        linearInd   = sub2ind([imSize(1),imSize(2)], grainBound{1}(:,1), grainBound{1}(:,2));
        contactIdx  = unique(imageGrainIdx(linearInd));
    elseif imSize(3) > 1
        contactIdx      = [];
        imageGrain      = reshape(imageGrain,[imSize(1), imSize(2), imSize(3)]);
        for iz = 1:imSize(3)
            % Dilate the grain first to get the outer boundary.
            imageGrain(:,:,iz) = bwmorph(imageGrain(:,:,iz),'dilate');
            imageGrain(:,:,iz) = bwmorph(imageGrain(:,:,iz),'dilate');
            % Find boundary of a grain
            grainBound  = bwboundaries(imageGrain(:,:,iz),'noholes');
            if length(grainBound) == 1
                linearInd   = sub2ind([imSize(1),imSize(2)], grainBound{1}(:,1), grainBound{1}(:,2));
                contactIdx  = [contactIdx; imageGrainIdx(linearInd)]; 
            end
        end     
    end
    
    % Find the unique index that is not 0 (the pore) and the itself.
    contactIdx  = unique(contactIdx);
    contactIdx  = contactIdx(contactIdx ~= 0);
    contactIdx  = contactIdx(contactIdx ~= grainNo(iGrain));
    nContact(iGrain) = length(contactIdx);
end



%% Measure the grain size    
for iGrain = 1:nGrain
    if iGrain/100 == round(iGrain/100)
    disp(['Step 6: Measure grain size (', num2str(iGrain),...
        '/',num2str(nGrain), ')'])
    end
    idxSingleGrain      = (imageGrainIdx(:) == grainNo(iGrain));
    imageGrain    = zeros(nTotalSize,1);
    imageGrain(idxSingleGrain) = 1;
    
    
    if imSize(3) == 1 % 2D-------------------------------------------------
        imageGrain    = reshape(imageGrain,[imSize(2), imSize(1)]);
        [idxXX, idxYY] = find(imageGrainIdx == grainNo(iGrain));
        
        % Obtain principal direction
        [coeff] = pca([idxXX idxYY]);
        
        nVec = size(coeff,2);
        X0 = allGrainCentroid(iGrain,1);
        Y0 = allGrainCentroid(iGrain,2);

        for iVec = 1:nVec
            clear temp*

            % Find index of a straight line in both direction on principal axis
            a = coeff(1,iVec);
            b = coeff(2,iVec);

            if a >= b
                tempX1 = [X0:imSize(2)]';
                tempX2 = [X0:-1:1]';
                nL1    = ones(length(tempX1),1);
                nL2    = ones(length(tempX2),1);
                tempY1 = nL1.*round(b./a.*(tempX1-X0) + Y0);
                tempY2 = nL2.*round(b./a.*(tempX2-X0) + Y0);    
            elseif b > a
                tempY1 = [Y0:imSize(1)]';
                tempY2 = [Y0:-1:1]';
                nL1    = ones(length(tempY1),1);
                nL2    = ones(length(tempY2),1);
                tempX1 = nL1.*round(a./b.*(tempY1-Y0) + X0);
                tempX2 = nL2.*round(a./b.*(tempY2-Y0) + X0);        
            end

            % Initialization of lines for measurement
            nLine1 = length(tempX1);
            nLine2 = length(tempX2);
            tempLine1 = zeros(nLine1,1);
            tempLine2 = zeros(nLine2,1);

            for iLine1 = 1:nLine1
                try % instead of checking boundary
                tempLine1(iLine1,1) = imageGrain(tempY1(iLine1),...
                                                       tempX1(iLine1));
                end
            end

            for iLine2 = 1:nLine2
                try
                tempLine2(iLine2,1) = imageGrain(tempY2(iLine2),...
                                                       tempX2(iLine2));
                end
            end

            % Find the boundary wihtin the lines
            bound1 = find(tempLine1 == 0, 1, 'first');
            if isempty(bound1)
                bound1 = nLine1;
            end

            bound2 = find(tempLine2 == 0, 1, 'first');
            if isempty(bound2)
                bound2 = nLine2;
            end

            % Calculate grain size to the boundary
            dX1 = tempX1(bound1) - tempX1(1);
            dY1 = tempY1(bound1) - tempY1(1);
            dX2 = tempX2(bound2) - tempX2(1);
            dY2 = tempY2(bound2) - tempY2(1);
            
            % Compute 2 radius from PCA at a time
            grainRadius(iGrain,2.*(iVec-1) + 1) ...
                = sqrt(dX1^2 + dY1^2);
            grainRadius(iGrain,2.*(iVec-1) + 2) ...
                = sqrt(dX2^2 + dY2^2);

            grainAzimuth(iGrain,2.*(iVec-1) + 1) ...
                = atan(dX1/dY1)*180/pi;
            grainAzimuth(iGrain,2.*(iVec-1) + 2) ...
                = atan(dX2/dY2)*180/pi;
            
            if qcPlot
                if imSize(3) == 1 
                hold on
                plot(tempX1(1:bound1),tempY1(1:bound1),'r','LineWidth',2)
                plot(tempX2(1:bound2),tempY2(1:bound2),'r','LineWidth',2)
                end
            end
            
        end
    
        
    elseif imSize(3) > 1 % 3D----------------------------------------------
        imageGrain    ...
            = reshape(imageGrain,[imSize(1), imSize(2), imSize(3)]);

        [idxXX, idxYY, idxZZ] = find(imageGrainIdx == grainNo(iGrain));

        % principal component analysis to get principal direction (each column
        % = one principal component
        [coeff] = pca([idxXX idxYY idxZZ]);
        nVec = size(coeff,2);

        % Centroid
        X0 = allGrainCentroid(iGrain,1);
        Y0 = allGrainCentroid(iGrain,2);
        Z0 = allGrainCentroid(iGrain,3);

        for iVec = 1:nVec
            clear temp*

            % Obtain the value of secondary and tertiary direction
            a = coeff(1,iVec);
            b = coeff(2,iVec);
            c = coeff(3,iVec);

            [~,idxMax] = max([a,b,c]);

            if idxMax == 1;
                tempX1 = [X0:imSize(1)]';
                tempX2 = [X0:-1:1]';
                nL1    = ones(length(tempX1),1);
                nL2    = ones(length(tempX2),1);
                tempY1 = nL1.*round(b./a.*(tempX1-X0) + Y0);
                tempY2 = nL2.*round(b./a.*(tempX2-X0) + Y0);    
                tempZ1 = nL1.*round(c./a.*(tempX1-X0) + Z0);
                tempZ2 = nL2.*round(c./a.*(tempX2-X0) + Z0);   
            elseif idxMax == 2;
                tempY1 = [Y0:imSize(2)]';
                tempY2 = [Y0:-1:1]';
                nL1    = ones(length(tempY1),1);
                nL2    = ones(length(tempY2),1);
                tempX1 = nL1.*round(a./b.*(tempY1-Y0) + X0);
                tempX2 = nL2.*round(a./b.*(tempY2-Y0) + X0); 
                tempZ1 = nL1.*round(c./b.*(tempY1-Y0) + Z0);
                tempZ2 = nL2.*round(c./b.*(tempY2-Y0) + Z0);             
            elseif idxMax == 3;
                tempZ1 = [Z0:imSize(3)]';
                tempZ2 = [Z0:-1:1]';
                nL1    = ones(length(tempZ1),1);
                nL2    = ones(length(tempZ2),1);
                tempX1 = nL1.*round(a./c.*(tempZ1-Z0) + X0);
                tempX2 = nL2.*round(a./c.*(tempZ2-Z0) + X0); 
                tempY1 = nL1.*round(c./c.*(tempZ1-Z0) + Y0);
                tempY2 = nL2.*round(c./c.*(tempZ2-Z0) + Y0);               
            end

            % Initialization of lines for measurement
            nLine1 = length(tempX1);
            nLine2 = length(tempX2);
            tempLine1 = zeros(nLine1,1);
            tempLine2 = zeros(nLine2,1);

            for iLine1 = 1:nLine1
                try % instead of checking boundary
                    tempLine1(iLine1,1) ...
                        = imageGrain(tempY1(iLine1),...
                                           tempX1(iLine1),...
                                           tempZ1(iLine1));
                end
            end

            for iLine2 = 1:nLine2
                try
                    tempLine2(iLine2,1) ...
                        = imageGrain(tempY2(iLine2),...
                                           tempX2(iLine2),...
                                           tempZ2(iLine2));
                end
            end

            % Find the boundary wihtin the lines
            bound1 = find(tempLine1 == 0, 1, 'first');
            if isempty(bound1)
                bound1 = nLine1;
            end

            bound2 = find(tempLine2 == 0, 1, 'first');
            if isempty(bound2)
                bound2 = nLine2;
            end

            % Calculate grain size to the boundary
            dX1 = tempX1(bound1) - tempX1(1);
            dY1 = tempY1(bound1) - tempY1(1);
            dZ1 = tempZ1(bound1) - tempZ1(1);
            dX2 = tempX2(bound2) - tempX2(1);
            dY2 = tempY2(bound2) - tempY2(1);
            dZ2 = tempZ2(bound2) - tempZ2(1);
            
            grainRadius(iGrain,2.*(iVec-1) + 1) ...
                = sqrt((dX1)^2 + (dY1)^2 + (dZ1)^2);
            grainRadius(iGrain,2.*(iVec-1) + 2) ...
                = sqrt((dX2)^2 + (dY2)^2 + (dZ2)^2);
                 
            grainAzimuth(iGrain,2.*(iVec-1) + 1) ...
                = acos(dZ1/grainRadius(iGrain,2.*(iVec-1) + 1));
            grainAzimuth(iGrain,2.*(iVec-1) + 2) ...
                = acos(dZ2/grainRadius(iGrain,2.*(iVec-1) + 2));
            
            grainInclination(iGrain,2.*(iVec-1) + 1) ...
                = atan(dY1/dX1)*180/pi;
            grainInclination(iGrain,2.*(iVec-1) + 2) ...
                = atan(dY2/dX2)*180/pi;

        end
    end
end
   


%% Output
% Clean 0 data
idxNonZero          = all(grainRadius,2);
grainRadius         = grainRadius(idxNonZero,:);
grainCentroid       = allGrainCentroid(idxNonZero,:);
grainVolume         = allGrainVolume(idxNonZero,:);
grainAzimuth        = grainAzimuth(idxNonZero,:);
grainInclination    = grainInclination(idxNonZero,:);
nContact            = nContact(idxNonZero,:);
grainSurfaceArea    = grainSurfaceArea(idxNonZero,:);

%% Boundary condition Exclude grains at the boundary 
if bc
    radius      = max(grainRadius,[],2);
    nGrain      = length(radius);
    idxAccept   = [];
    idxReject   = [];

    if imSize(3) == 1 %2D------------------------------------------------------
        % Define all the corner points
        c1 = [1, 1];
        c2 = [imSize(1), 1];
        c3 = [1, imSize(1)];
        c4 = [imSize(1), imSize(1)];

        % Calculate the distance from a point to a line for each centroid
        for iGrain = 1:nGrain
            p = grainCentroid(iGrain,:);
            distanceLine1 = abs(det([c2 - c1; p - c2])/sqrt(sum((c2 - c1).^2)));
            distanceLine2 = abs(det([c3 - c1; p - c3])/sqrt(sum((c3 - c1).^2)));
            distanceLine3 = abs(det([c4 - c3; p - c3])/sqrt(sum((c4 - c3).^2)));
            distanceLine4 = abs(det([c4 - c2; p - c2])/sqrt(sum((c4 - c2).^2)));
            if (distanceLine1 > radius(iGrain) && distanceLine2 > radius(iGrain) &&...
                    distanceLine3 > radius(iGrain) && distanceLine4 > radius(iGrain))
                idxAccept = [idxAccept, iGrain];
            else
                idxReject = [idxReject, iGrain];
            end

        end    

    elseif imSize(3) > 1 %3D---------------------------------------------------
        % Define the corner points
        corner = [1 1 1; 
                  1 imSize(1) 1; 
                  1 1 imSize(1); 
                  1 imSize(1) imSize(1);
                  imSize(1) 1 1; 
                  imSize(1) imSize(1) 1; 
                  imSize(1) 1 imSize(1); 
                  imSize(1) imSize(1) imSize(1)]; 

        % Define 6 different planes from corners point
        plane = [1 2 3;
                 1 3 5;
                 3 4 7;
                 2 4 6;
                 1 2 5;
                 5 6 7];
        % Calculate the distance from a point to a line for each centroid      
        for iGrain = 1:nGrain
            point = grainCentroid(iGrain,:);

            for iPlane = 1:6
                c1 = corner(plane(iPlane,1),:);
                c2 = corner(plane(iPlane,2),:);
                c3 = corner(plane(iPlane,3),:);

                normal = cross(c1 - c2, c1 - c3);

                d = dot(normal, c1);

                % Find the closest distance from a point to plane
                distance(iGrain,iPlane) = (dot(normal, point) - d)./sqrt(dot(normal,normal));
            end

            if min(abs(distance(iGrain,:))) > radius(iGrain)
                idxAccept = [idxAccept, iGrain];
            else
                idxReject = [idxReject, iGrain];
            end

        end    

    end
    
    % Screen the data to the number 
    grainRadius         = grainRadius(idxAccept,:);
    grainCentroid       = grainCentroid(idxAccept,:);
    grainVolume         = grainVolume(idxAccept,:);
    grainAzimuth        = grainAzimuth(idxAccept,:);
    grainInclination    = grainInclination(idxAccept,:);
    nContact            = nContact(idxAccept,:);
    grainSurfaceArea    = grainSurfaceArea(idxAccept,:);
    
    
    if qcPlot
        if imSize(3) == 1 
        scatter(grainCentroid(:,1), grainCentroid(:,2),'go')
        end
    end
end

    


