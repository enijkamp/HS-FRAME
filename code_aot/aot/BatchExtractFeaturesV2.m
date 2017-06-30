function SPM_SUM3 = BatchExtractFeaturesV2(codebook, imagesBatch)



load 'storedExponentialModel'; % load in exponential model
%resizeFactor = 1; % resize the input images
numResolution = 3; % number of resolutions to search for in detection stage
sizeTemplatex = 144; sizeTemplatey = 144; % Object template size. size...x denote height, and size...y denotes width.
templateSize = [sizeTemplatex sizeTemplatey];
partSizeX = floor(sizeTemplatex/3); % part template size
partSizeY = floor(sizeTemplatex/3);

% manually label the image patch to initialize the template:
starting = 1; % initilize by single image learning from I{starting}
originalResolution = 3; % original resolution is the one at which the imresize factor = 1, see 11th line beneath this line
startx = 1; endx = startx + sizeTemplatex - 1; % % bounding box of the first object
starty = 1; endy = starty + sizeTemplatey - 1;

%% These are parameters associated with Morphable active basis model

% to be frequently adjusted:
numIteration = 2;  % number of iterations
partRotationRange = 2*(-2:2); % absolute part rotation (rotation of partial templates)
numPartRotate = length(partRotationRange);
maxPartRelativeRotation = 2;
resolutionShiftLimit = 1;
minRotationDif = (sin(maxPartRelativeRotation*pi/numOrient)-sin(0))^2 + (cos(maxPartRelativeRotation*pi/numOrient)-cos(0))^2 + 1e-10;
rotationRange = 2*(-1:1); % whole object rotation
numRotate = length(rotationRange);

% to be occationally adjusted
numElement = 20; % number of Gabors in active basis
locationPerturbFraction = .2; % part perturbation
locationShiftLimit = 2; % shift in normal direction = locationShiftLimit*subsample pixels
orientShiftLimit = 1; % shift in orientation
subsampleS2 = 3; subsampleM2 = 1;
partMarginX = floor(partSizeX*.75);
partMarginY = floor(partSizeY*.75);

% to stay fixed
epsilon = .1; % allowed correlation between selected Gabors
subsample = 1; % subsample in computing MAX1 and SUM2 maps
Correlation = CorrFilter(allFilter, epsilon); % correlation between filters

PartLocX0 = 1:partSizeX:templateSize(1)-partSizeX+1;
PartLocY0 = 1:partSizeY:templateSize(2)-partSizeY+1;
partLocRange = floor(sqrt(partSizeX*partSizeY)*locationPerturbFraction);
numCandPart = length(PartLocX0) * length(PartLocY0);
PartLocX = zeros(numCandPart,1);
PartLocY = zeros(numCandPart,1);
iPart = 1;
for x = PartLocX0
    for y = PartLocY0
        PartLocX(iPart) = x;
        PartLocY(iPart) = y;
        iPart = iPart + 1;
    end
end









SPM_SUM3 = [];

for img= 1:length(imagesBatch)
    
    disp([' start processing image ' num2str(img)]); tic    
    
%    tmpIm = imread(fullfile(imgPath, categoryNames{iClass},imgList(img).name));
%     
%     if size(tmpIm,3) == 3
%         tmpIm = rgb2gray(tmpIm);
%     end

    tmpIm = imagesBatch{img};
    I = imresize(single(tmpIm), [sizeTemplatex,sizeTemplatey], 'nearest');
    ImageMultiResolution = cell(1,numResolution);
    for j=1:numResolution
        resolution = .6+(j-1)*.2; % so that .8+(originalResolution-1)*.2 = 1
        ImageMultiResolution{j} = imresize(I, resolution, 'nearest');  % images at multiple resolutions
    end
    
    SUM1mapFind = ApplyFilterfft(ImageMultiResolution, allFilter,...
        localHalfx, localHalfy, thresholdFactor); % filtering images at multiple resolutions
    mexc_Sigmoid(saturation, SUM1mapFind);
    
    MAX1map = cell(size(SUM1mapFind));
    M1Trace = cell(size(SUM1mapFind));
    for iRes = 1:numResolution
        [MAX1map(iRes,:) M1Trace(iRes,:) M1RowShift M1ColShift M1OriShifted] = ...
            mexc_ComputeMAX1( numOrient, SUM1mapFind(iRes,:), locationShiftLimit,...
            orientShiftLimit, 1);
    end
    
    SUM3map = cell(length(rotationRange), numResolution, length(codebook));
    
    for c = 1:length(codebook)
        
        selectedPart = find(codebook{c}.PartOnOff);
        
        %% compute SUM2 maps for non-overlapping parts
        SUM2map = cell(numPartRotate,numCandPart,numResolution);
        S2T = cell( numPartRotate, numCandPart );
        for iPart = 1:numCandPart
            for r = 1:length(partRotationRange)
                S2T{r,iPart} = struct( 'selectedRow',single(codebook{c}.allSelectedx{iPart,r}(:)),...
                    'selectedCol', single(codebook{c}.allSelectedy{iPart,r}(:)),...
                    'selectedOri', single(codebook{c}.allSelectedOrient{iPart,r}(:)),...
                    'selectedScale', zeros(length(codebook{c}.allSelectedx{iPart,r}),1,'single'),...
                    'selectedLambda', single(codebook{c}.selectedlambda{iPart}(:)),...
                    'selectedLogZ', single(codebook{c}.selectedLogZ{iPart}(:)) );
            end
        end
        
        for iRes = 1:numResolution
            tmpS2 = mexc_ComputeSUM2( numOrient,...
                MAX1map(iRes,:), S2T(:), subsampleS2 );
            SUM2map(:,:,iRes) = reshape(tmpS2,[numPartRotate numCandPart]);
        end
        
        %% compute SUM2 maps for overlapping parts
        largerSUM2map = cell(numPartRotate,numCandPart,numResolution);
        largerS2T = cell( numPartRotate, numCandPart );
        for iPart = 1:numCandPart
            for r = 1:length(partRotationRange)
                largerS2T{r,iPart} = struct( 'selectedRow',single(codebook{c}.largerAllSelectedx{iPart,r}(:)),...
                    'selectedCol', single(codebook{c}.largerAllSelectedy{iPart,r}(:)),...
                    'selectedOri', codebook{c}.largerAllSelectedOrient{iPart,r}(:),...
                    'selectedScale', zeros(length(codebook{c}.largerAllSelectedx{iPart,r}),1,'single'),...
                    'selectedLambda', single(codebook{c}.largerSelectedlambda{iPart}(:)),...
                    'selectedLogZ', single(codebook{c}.largerSelectedLogZ{iPart}(:)));
            end
        end
        
        % compute SUM2
        for iRes = 1:numResolution
            tmpS2 = mexc_ComputeSUM2( numOrient,...
                MAX1map(iRes,:), largerS2T(:), subsampleS2 );
            largerSUM2map(:,:,iRes) = reshape(tmpS2,[numPartRotate numCandPart]);
        end
        
        % compute MAX2 maps for overlapping parts (local maximization w.r.t. translation and rotation)
        templateAffinityMatrix = cell(numPartRotate,numCandPart);
        for iPart = 1:numCandPart
            for r1 = 1:length(partRotationRange)
                angle1 = pi/numOrient * partRotationRange(r1);
                templateAffinityMatrix{r1,iPart} = [];
                jPart = iPart;
                for r2 = 1:length(partRotationRange)
                    angle2 = pi/numOrient*partRotationRange(r2);
                    if (sin(angle1) - sin(angle2))^2 + (cos(angle1)-cos(angle2))^2 <= minRotationDif
                        templateAffinityMatrix{r1,iPart} = int32( [templateAffinityMatrix{r1,iPart} r2+(jPart-1)*numPartRotate-1] );
                    end
                end
            end
        end
        
        largerMAX2map = cell(size(SUM2map));
        largerMAX2LocTrace = cell(size(SUM2map));
        largerMAX2TransformTrace = cell(size(SUM2map));
        
        for iRes = 1:numResolution
            [tmpMAX2 tmpMAX2LocTrace tmpMAX2TransformTrace...
                M2RowColShift] = mexc_ComputeMAX2( templateAffinityMatrix(:), ...
                largerSUM2map(:,:,iRes), ...
                locationPerturbFraction, ...
                int32(sqrt(partSizeX*partSizeY)*ones(numPartRotate*numCandPart,1)/subsampleS2), subsampleM2 );
            largerMAX2map(:,:,iRes) = reshape(tmpMAX2,[numPartRotate numCandPart]);
            largerMAX2LocTrace(:,:,iRes) = reshape(tmpMAX2LocTrace,[numPartRotate numCandPart]);
            largerMAX2TransformTrace(:,:,iRes) = reshape(tmpMAX2TransformTrace,[numPartRotate numCandPart]);
        end
        
        % fill in the MAX2 map for non-overlapping parts
        
        tmpMAX2map = mexc_FakeMAX2( SUM2map, largerMAX2LocTrace, largerMAX2TransformTrace,...
            templateAffinityMatrix, int32(ones(numel(SUM2map),1) * sqrt(partSizeX*partSizeY)/subsampleS2), M2RowColShift );
        
        % max over resolution
        MAX2map = tmpMAX2map;
        tmpLargerMAX2map = largerMAX2map;
        MAX2ResolutionTrace = cell(size(largerMAX2map)); % to initialize
        for iRes = 1:numResolution
            current_size = size( MAX2map{1,1,iRes} );
            for j = 1:size(MAX2map,1)
                for k = 1:size(MAX2map,2)
                    map = -1e10 * ones( current_size, 'single' ); % an auxiliary variable, to find the ARGMAX resolution
                    MAX2ResolutionTrace{j,k,iRes} = int32( -1 * ones( current_size ) );
                    for jRes = 1:numResolution
                        if abs(jRes-iRes) <= resolutionShiftLimit
                            ref = largerMAX2map{j,k,jRes};
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
        
        % MAX2map = largerMAX2map;
        
        bestS3Loc = -1; bestRes = 0; bestRot = 0;
        
        S3T = cell(length(rotationRange),1);
        for r = 1:length(rotationRange) % this is the rotation of the S3 template
            rot = rotationRange(r);
            
            % MAX2score, Fx, Fy are only for temperary storage
            MAX2score = single(zeros(1, numResolution));
            Fx = zeros(1, numResolution);
            Fy = zeros(1, numResolution);
            
            % compute SUM3 maps
            selectedTransform = zeros(length(selectedPart),1,'single');
            for j = 1:length(selectedPart)
                selectedTransform(j) = find( codebook{c}.allS3SelectedOri(r,j) == partRotationRange );
            end
            S3T{r} = struct( 'selectedRow',single(floor(.5+codebook{c}.allS3SelectedRow(r,selectedPart)/subsampleM2/subsampleS2)),...
                'selectedCol', single(floor(.5+codebook{c}.allS3SelectedCol(r,selectedPart)/subsampleM2/subsampleS2)),...
                'selectedInd', single(selectedPart) - 1,...
                'selectedTransform', selectedTransform - 1,...
                'selectedLambda', ones(length(selectedPart),1,'single'),...
                'selectedLogZ', single( 0*ones(length(selectedPart),1) ) );
            
            for iRes = 1:numResolution
                tmpM2 = MAX2map(:,:,iRes);
                SUM3map(r, iRes, c) = mexc_ComputeSUM3( tmpM2(:), S3T(r), 1, numPartRotate );
            end
        end
    end
    
    SUM3map = reshape(SUM3map, [numResolution*numRotate*length(codebook), 1]);
    SPM_SUM3 = [SPM_SUM3 ; featureSPM(SUM3map,spm_numLayers,spm_threshold)];
end


