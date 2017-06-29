function  [MMAX3]=detectObject(imageLoaded, clusters, c, iImg, numPart, numFilter, sx, sy, nOrient, part_sx, part_sy,...  % general parameters
    relativePartLocationRange, argmaxMethod, allTemplateAffinityMatrix, resolutionShiftLimit, rotationRange, partRotationRange, PartLocX, PartLocY,...   % parameters for local max
    detectedCroppedSavingFolder, morphedCroppedSavingFolder, boundingBoxSavingFolder, RatioDisplacementSUM3)

numResolution = length(imageLoaded.ImageMultiResolution);
numRotate = length(rotationRange);
numPartRotate = length(partRotationRange);
% compute SUM2 map
tic 
SUM2map=cell(numPartRotate, numPart,numResolution);
for iRes = 1:numResolution
    tmpSUM2=sparseFRAME_SUM2_part(single(imageLoaded.allSizex(iRes)), single(imageLoaded.allSizey(iRes)), single(numFilter), clusters(c).S2T(:), imageLoaded.MAX1mapFind(iRes,:));
    SUM2map(:,:,iRes)=reshape(tmpSUM2, [numPartRotate, numPart]);
end
disp(['finished SUM2 map of image ' num2str(iImg) ' by moving parts in cluster ' num2str(c) ': '  num2str(toc) ' seconds']);

% compute MAX2 map
tic
tmpMAX2map=cell(size(SUM2map));
MAX2LocTrace=cell(size(SUM2map));
MAX2TransformTrace=cell(size(SUM2map));

for iRes = 1:numResolution
    [tmpMAX2, tmpMAX2LocTrace, tmpMAX2TransformTrace,M2RowColShift] = ...
        mexc_sparseFRAME_MAX2( allTemplateAffinityMatrix(:), SUM2map(:,:,iRes), relativePartLocationRange, argmaxMethod, 1);
    
    tmpMAX2map(:,:,iRes) = reshape(tmpMAX2,[numPartRotate numPart]);
    MAX2LocTrace(:,:,iRes) = reshape(tmpMAX2LocTrace,[numPartRotate numPart]);   % start  from 0
    MAX2TransformTrace(:,:,iRes) = reshape(tmpMAX2TransformTrace,[numPartRotate numPart]);      % start  from 0
end

% continue to compute MAX2 map (max over resolution)
MAX2ResolutionTrace = cell(size(tmpMAX2map)); % to initialize
MAX2map = tmpMAX2map;
for iRes = 1:numResolution
    current_size = size( MAX2map{1,1,iRes} );
    for j = 1:size(MAX2map,1)
        for k = 1:size(MAX2map,2)
            map = -1e10 * ones( current_size, 'single' ); % an auxiliary variable, to find the ARGMAX resolution
            MAX2ResolutionTrace{j,k,iRes} = int32( -1 * ones( current_size ) );
            for jRes = 1:numResolution
                if abs(jRes-iRes) <= resolutionShiftLimit
                    ref = tmpMAX2map{j,k,jRes};
                    ref = imresize(ref,current_size,'nearest');
                    ind = ref > map;
                    map(ind) = ref(ind);
                    tocopy = tmpMAX2map{j,k,jRes};
                    tocopy = imresize(tocopy,current_size,'nearest');
                    MAX2map{j,k,iRes}(ind) = tocopy(ind);
                    MAX2ResolutionTrace{j,k,iRes}(ind) = jRes - 1; % start from 0
                end
            end
        end
    end
end
disp(['finished MAX2 map of image ' num2str(iImg) ' by template of cluster ' num2str(c) ': '  num2str(toc) ' seconds']);


%%%% align the object by detection
%% compute SUM3map
tic
SUM3map = cell(numRotate, numResolution);
for r = 1:length(rotationRange) % this is the rotation of the S3 template
    for iRes = 1:numResolution
%         tmpM2 = MAX2map(:,:,iRes);
%         tmpS3 = mexc_ComputeSUM3( tmpM2(:), clusters(c).S3T(r), 1, numPartRotate);
%         SUM3map(r, iRes) = {tmpS3{1} - clusters(c).logZ};  % for learning mixture model, we need log of Z.
      
        tmpM2 = MAX2map(:,:,iRes);      
       % SUM3map(r, iRes) = mexc_ComputeSUM3_logZ( tmpM2(:), clusters(c).S3T(r), 1, numPartRotate, single(clusters(c).logZ));
        SUM3map(r, iRes) = mexc_ComputeSUM3_logZ_partial( tmpM2(:), clusters(c).S3T(r), 1, numPartRotate, single(clusters(c).logZ), single(RatioDisplacementSUM3));
        
    end
end

MAX3score = single(zeros(1, numResolution));  % maximum log-likelihood score at each resolution
allFx = zeros(1, numResolution); % detected location at each resolution
allFy = zeros(1, numResolution); % detected location at each resolution
MMAX3 = -1e10;
for r = 1:length(rotationRange)
    for iResolution=1:numResolution
        [MAX3score(iResolution), ind]=max(SUM3map{r,iResolution}(:));
        [allFx(iResolution), allFy(iResolution)] = ind2sub([imageLoaded.allSizex(iResolution), imageLoaded.allSizey(iResolution)],ind);
    end
    
    [maxOverResolution, ind] = max(MAX3score);   % most likely resolution
    if (MMAX3 < maxOverResolution)
        MMAX3 = maxOverResolution;
        Mrot = rotationRange(r);
        Mind = ind;  % best resoltion
        MFx = allFx(ind);
        MFy = allFy(ind);
    end
end

bestRotInd = find(Mrot==rotationRange); % best object rotation index, starting from 1
disp(['finished SUM3 map of image ' num2str(iImg) ' by template of cluster ' num2str(c) ' and detect the object: '  num2str(toc) ' seconds']);

%%% from here, we have located the object.

%%
drawBoundingBox;   % a section of codes to draw bounding box

%% save the croped detected object
cropedImage=single(zeros(sx, sy));
Ccopy(cropedImage, single(imageLoaded.ImageMultiResolution{Mind}), MFx, MFy, floor(sx/2), floor(sy/2), sx, sy, imageLoaded.allSizex(Mind), imageLoaded.allSizey(Mind), -Mrot*pi/nOrient);

gLow = min(cropedImage(:));
gHigh = max(cropedImage(:));
img_tem = (cropedImage-gLow)/(gHigh-gLow);
imwrite(img_tem,fullfile(detectedCroppedSavingFolder,['detected-cluster-' num2str(c) '-img-' num2str(iImg,'%04d') '.png']));

%%%% save the cropped morphed image
% set default values of some output variables
bestPartRes = Mind;   % it may be changed for each part

cropedMorphedImage=single(zeros(sx, sy));
for iPart = 1:numPart
    
    r = clusters(c).S3T{bestRotInd}.selectedTransform(iPart)+1; % the transform index now starts from 1 not 0  (rotation)
    Fx = MFx + round(clusters(c).S3T{bestRotInd}.selectedRow(iPart));
    Fy = MFy + round(clusters(c).S3T{bestRotInd}.selectedCol(iPart));
    imagesize = size(MAX2map{r,iPart,Mind});
    
    if Fx >= 1 && Fx <= imagesize(1) && Fy >= 1 && Fy <= imagesize(2)
        
        tmp = MAX2ResolutionTrace{r,iPart,Mind};
        bestPartRes = tmp(Fx,Fy) + 1; % best part resolution, and the index starts from 1
        current_size = size(tmp);
        tmp = MAX2LocTrace{r,iPart,bestPartRes};
        new_size = size(tmp);
        % compute new coordinate in image with best resolution
        Fx = round(Fx*double(new_size)/current_size);
        Fy = round(Fy*double(new_size)/current_size);
        
        if Fx >= 1 && Fx <= size(tmp,1) && Fy >= 1 && Fy <= size(tmp,2);
            translationInd = tmp(Fx,Fy) + 1;
        else
            translationInd = floor(size(M2RowColShift,1)/2);
        end
        
        tmp = MAX2TransformTrace{r,iPart,bestPartRes};
        if Fx >= 1 && Fx <= size(tmp,1) && Fy >= 1 && Fy <= size(tmp,2);
            transformInd = tmp(Fx,Fy) + 1;
        else
            transformInd = floor(numPartRotate/2) + 1;
        end
        
        actualPartRotationInd = transformInd - numPartRotate*(ceil(double(transformInd)/numPartRotate)-1);
        
        Fx = floor( Fx + M2RowColShift(translationInd,1));
        Fy = floor( Fy + M2RowColShift(translationInd,2));
        
    else
        actualPartRotationInd = r;
    end
    
    % find the part location at the higher resolution
    % %                 Fx = (Fx-1 + .5);
    % %                 Fy = (Fy-1 + .5);
    % crop the feature patch that is registered to the part template
    % %                 tmpSUM1mapLearn = mexc_CropInstance(SUM1mapFind(bestPartRes,:),Fx,Fy,...
    % %                           partRotationRange(actualPartRotationInd),tScale,1,...
    % %                           partOutRow{actualPartRotationInd},partOutCol{actualPartRotationInd},...
    % %                           nOrient,1,part_sx,part_sy);
    % %
    % %                 for o = 1:numOrient
    % %                     tmpSUM1mapFind{1,o}(PartLocX(iPart)-1+(1:partSizeX),PartLocY(iPart)-1+(1:partSizeY)) = tmpSUM1mapLearn{o};
    % %                 end
    
    cropedImage=single(zeros(part_sx, part_sy));
    Ccopy(cropedImage, single(imageLoaded.ImageMultiResolution{bestPartRes}), Fx, Fy, floor(part_sx/2), floor(part_sy/2), part_sx, part_sy, ...
        imageLoaded.allSizex(bestPartRes), imageLoaded.allSizey(bestPartRes), -partRotationRange(actualPartRotationInd)*pi/nOrient);
    
    cropedMorphedImage(PartLocX(iPart)-1+(1:part_sx),PartLocY(iPart)-1+(1:part_sy))=cropedImage;
    
    
    
    % % %                 %%% crop the cropped feature patch that is registered to the part template
    % % %                 croppedPartSUM1map = cell(numFilter, 1);
    % % %                 for iFilter = 1:numFilter
    % % %                     croppedPartSUM1map{iFilter, 1} = single(zeros(part_sx, part_sy));
    % % %                 end
    % % %                 for iScale=1:nScaleGabor
    % % %                   for iF = 1: 2 % sine or cosine part
    % % %                     for orient = 1 : nOrient
    % % %                         orient1 = orient - partRotationRange(actualPartRotationInd);
    % % %                         if (orient1 > nOrient)
    % % %                             orient1 = orient1 - nOrient;
    % % %                         end
    % % %                         if (orient1 <= 0)
    % % %                             orient1 = orient1 + nOrient;
    % % %                         end
    % % %                         Ccopy(croppedPartSUM1map{orient+((iScale-1)*2+(iF-1))*nOrient, 1}, imageLoaded.SUM1mapFind{bestPartRes, orient1+((iScale-1)*2+(iF-1))*nOrient}, Fx, Fy, floor(part_sx/2), floor(part_sy/2), ...
    % % %                             part_sx, part_sy, imageLoaded.allSizex(bestPartRes), imageLoaded.allSizey(bestPartRes), -partRotationRange(actualPartRotationInd)*pi/nOrient);
    % % %                     end
    % % %                   end
    % % %                 end
    % % %                 for iScale=1:nScaleDoG
    % % %                       Ccopy(croppedPartSUM1map{nScaleGabor*nOrient*2+iScale, 1},imageLoaded.SUM1mapFind{bestPartRes,nScaleGabor*nOrient*2+iScale}, Fx, Fy, floor(part_sx/2), floor(part_sy/2),...
    % % %                           part_sx, part_sy, imageLoaded.allSizex(bestPartRes), imageLoaded.allSizey(bestPartRes), -partRotationRange(actualPartRotationInd)*pi/nOrient);
    % % %                 end
    % % %
    % % %                 for iFilter = 1:numFilter
    % % %                       croppedMorphedSUM1{iFilter,1}(PartLocX(iPart)-1+(1:part_sx),PartLocY(iPart)-1+(1:part_sy)) = croppedPartSUM1map{iFilter};
    % % %                 end
    %%%%%%%%%
    
end

% save morphed image
gLow = min(cropedMorphedImage(:));
gHigh = max(cropedMorphedImage(:));
img_tem = (cropedMorphedImage-gLow)/(gHigh-gLow);
imwrite(img_tem,fullfile(morphedCroppedSavingFolder,['morphed-cluster-' num2str(c) '-img-' num2str(iImg,'%04d') '.png']));

% save variable
% % %             save(fullfile(morphedCroppedSavingFolder,['morphed-cluster-' num2str(c) '-img-' num2str(iImg,'%04d') '.mat']), 'cropedMorphedImage', 'croppedMorphedSUM1');
save(fullfile(morphedCroppedSavingFolder,['morphed-cluster-' num2str(c) '-img-' num2str(iImg,'%04d') '.mat']), 'cropedMorphedImage');