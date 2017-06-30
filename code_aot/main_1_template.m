% clear
clear;
close all;

% code
addpath('aot');

% config
para = config();

% mex
prev = cd('aot');
compileMex();
cd(prev);

% run
tic_toc = zeros(length(para.task_ids),1);

%% Init

%% These are the parameters associated with data preprocessing, to adjust from category to category

NumCluster = para.numCluster;

if (exist('output','dir'))
    delete('output/*.*'); 
else
    mkdir('output');
end

if ~exist('working','dir'), mkdir('working');end
for c = 1:NumCluster, mkdir(['output/' num2str(c)]);end
if ~exist('template','dir'), mkdir('template');end

load 'storedExponentialModel'; % load in exponential model
%resizeFactor = 1; % resize the input images
numResolution = para.numResolution; % number of resolutions to search for in detection stage
sizeTemplatex = 144; sizeTemplatey = 144; % Object template size. size...x denote height, and size...y denotes width.
templateSize = [sizeTemplatex sizeTemplatey];
partSizeX = floor(sizeTemplatex/2); % part template size
partSizeY = floor(sizeTemplatex/2);

% manually label the image patch to initialize the template:
starting = 1; % initilize by single image learning from I{starting}
originalResolution = 3; % original resolution is the one at which the imresize factor = 1, see 11th line beneath this line 
startx = 1; endx = startx + sizeTemplatex - 1; % % bounding box of the first object
starty = 1; endy = starty + sizeTemplatey - 1;

%% These are parameters associated with Morphable active basis model

% to be frequently adjusted:
numIteration = para.numIteration;
partRotationRange = para.partRotationRange;
maxPartRelativeRotation = para.maxPartRelativeRotation;
resolutionShiftLimit = para.resolutionShiftLimit;
rotationRange = para.rotationRange;

numPartRotate = length(para.partRotationRange);
numRotate = length(para.rotationRange);
minRotationDif = (sin(para.maxPartRelativeRotation*pi/numOrient)-sin(0))^2 + (cos(para.maxPartRelativeRotation*pi/numOrient)-cos(0))^2 + 1e-10;

% to be occationally adjusted
numElement = para.numElement;
locationPerturbFraction = para.locationPerturbFraction;
locationShiftLimit = para.locationShiftLimit;
orientShiftLimit = para.orientShiftLimit;
subsampleS2 = para.subsampleS2;
subsampleM2 = para.subsampleM2;

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


for task_id = para.task_ids
    categoryName = para.categoryNames{task_id};
    big_ticID = tic;

    imageFolder = [para.dataPath categoryName]; % folder of training images  
    imageName = dir([imageFolder '/*.jpg']);
    numImage = size(imageName, 1); % number of training images 

    save('working/Config');

    %% Load in training images and initialize SUM1 maps for learning
    I = cell(1, numImage);
    for img = 1 : numImage
        tmpIm = imread([imageFolder '/' imageName(img).name]); 
        if size(tmpIm,3) == 3
            tmpIm = rgb2gray(tmpIm);
        end
        I{img} = imresize(single(tmpIm), [sizeTemplatex,sizeTemplatey], 'nearest');
    end

    %% Generate multi-resolution images and compute SUM1, MAX1 maps
    disp('start filtering training images at all resolutions'); tic
    for img = 1:numImage
        ImageMultiResolution = cell(1,numResolution);
        for j=1:numResolution
            resolution = .6+(j-1)*.2; % so that .8+(originalResolution-1)*.2 = 1
            ImageMultiResolution{j} = imresize(I{img}, resolution, 'nearest');  % images at multiple resolutions
        end
        save(sprintf('working/ImageAndFeature_%d.mat',img),'ImageMultiResolution');

        SUM1mapFind = ApplyFilterfft(ImageMultiResolution, allFilter,...
            localHalfx, localHalfy, thresholdFactor); % filtering images at multiple resolutions
        mexc_Sigmoid(saturation, SUM1mapFind);
        save(sprintf('working/ImageAndFeature_%d.mat',img),'SUM1mapFind','-append');

        MAX1map = cell(size(SUM1mapFind));
        M1Trace = cell(size(SUM1mapFind));
        for iRes = 1:numResolution
            [MAX1map(iRes,:) M1Trace(iRes,:) M1RowShift M1ColShift M1OriShifted] = ...
                mexc_ComputeMAX1( numOrient, SUM1mapFind(iRes,:), locationShiftLimit,...
                orientShiftLimit, 1 );
        end
        save(sprintf('working/ImageAndFeature_%d.mat',img),'MAX1map','M1Trace','M1RowShift','M1ColShift','M1OriShifted','-append');
    end

    disp(['filtering time: ' num2str(toc) ' seconds']);


    clusters = cell(NumCluster, 1);
    MAX3scoreAll = rand(numImage, NumCluster);   % randomly assign members to different cluster
    [~, ImageCluster] = min(MAX3scoreAll, [], 2);

    purity = zeros(numIteration, 1);
    entropy = zeros(numIteration, 1);

    for c = 1:NumCluster

        clusters{c}.imageIdx = find(ImageCluster==c);
        clusters{c}.numImage = size(clusters{c}.imageIdx, 1);

        %% Prepare output variables for learning
        clusters{c}.selectedOrient = cell(numCandPart, 1);  % orientation and location of selected Gabors
        clusters{c}.selectedx = cell(numCandPart, 1);
        clusters{c}.selectedy = cell(numCandPart, 1);
        clusters{c}.selectedlambda = cell(numCandPart, 1); % weighting parameter for scoring template matching
        clusters{c}.selectedLogZ = cell(numCandPart, 1); % normalizing constant

        clusters{c}.largerSelectedOrient = cell(numCandPart, 1);  % orientation and location of selected Gabors
        clusters{c}.largerSelectedx = cell(numCandPart, 1);
        clusters{c}.largerSelectedy = cell(numCandPart, 1);
        clusters{c}.largerSelectedlambda = cell(numCandPart, 1); % weighting parameter for scoring template matching
        clusters{c}.largerSelectedLogZ = cell(numCandPart, 1); % normalizing constant

        clusters{c}.commonTemplate = cell( numCandPart, 1 );
        for i = 1:length(clusters{c}.commonTemplate)
            clusters{c}.commonTemplate{i} = single(zeros(partSizeX, partSizeY)); % template of active basis 
        end

        deformedTemplate = cell(numImage, 1); % templates for training images 
        for img = 1 : numImage
            deformedTemplate{img} = zeros( sizeTemplatex, sizeTemplatey, 'single' );  
        end

        clusters{c}.allSelectedx = cell(numCandPart, numRotate); 
        clusters{c}.allSelectedy = cell(numCandPart, numRotate);
        clusters{c}.allSelectedOrient = cell(numCandPart, numRotate);
        clusters{c}.largerAllSelectedx = cell(numCandPart, numRotate);
        clusters{c}.largerAllSelectedy = cell(numCandPart, numRotate);
        clusters{c}.largerAllSelectedOrient = cell(numCandPart, numRotate);

        %% Initialize learning from the starting image
        locationShiftLimit0 = 0; orientShiftLimit0 = 0; 
        clusters{c}.SUM1mapLearn0 = cell(clusters{c}.numImage, numOrient);
        for imgs = 1:clusters{c}.numImage;
            load( sprintf('working/ImageAndFeature_%d.mat',clusters{c}.imageIdx(imgs)) );
            for orient = 1 : numOrient
                clusters{c}.SUM1mapLearn0{imgs, orient} = single(zeros(sizeTemplatex, sizeTemplatey));
                sizex = size(SUM1mapFind{originalResolution,orient},1);
                sizey = size(SUM1mapFind{originalResolution,orient},2);
                Ccopy(clusters{c}.SUM1mapLearn0{imgs, orient}, SUM1mapFind{originalResolution,orient}, startx-1, starty-1, 0, 0, sizeTemplatex, sizeTemplatey, sizex, sizey, 0); 
            end
        end
        disp('start from single image learning'); 
        tic

        tmpSelectedx = zeros(numElement,1);
        tmpSelectedy = zeros(numElement,1);
        tmpSelectedo = zeros(numElement,1);
        tmpSelectedlambda = zeros(numElement,1);
        tmpSelectedlogz = zeros(numElement,1);
        tmpCommonTemplate = zeros(sizeTemplatex,sizeTemplatey,'single');
        CsharedSketch(numOrient, locationShiftLimit0, orientShiftLimit0, subsample, ... % about active basis  
            numElement, clusters{c}.numImage, sizeTemplatex, sizeTemplatey, clusters{c}.SUM1mapLearn0, ... % about training images 
            halfFilterSize, Correlation, allSymbol(1, :), ... % about filters
            numStoredPoint, storedlambda, storedExpectation, storedLogZ, ... % about exponential model 
            tmpSelectedo, tmpSelectedx, tmpSelectedy, tmpSelectedlambda, tmpSelectedlogz, ... % learned parameters
            tmpCommonTemplate,deformedTemplate); % learned templates


        % split the object template into non-overlapping partial templates
        for iPart = 1:numCandPart
            ind = find( tmpSelectedx >= PartLocX(iPart) & tmpSelectedx < PartLocX(iPart) + partSizeX & ...
                    tmpSelectedy >= PartLocY(iPart) & tmpSelectedy < PartLocY(iPart) + partSizeY );
            clusters{c}.selectedOrient{iPart} = single( tmpSelectedo(ind) );
            clusters{c}.selectedx{iPart} = single( floor( tmpSelectedx(ind) - PartLocX(iPart) ) );
            clusters{c}.selectedy{iPart} = single( floor( tmpSelectedy(ind) - PartLocY(iPart) ) );
            clusters{c}.selectedlambda{iPart} = single( tmpSelectedlambda(ind) );
            clusters{c}.selectedLogZ{iPart} = single( tmpSelectedlogz(ind) );
            filename = sprintf('output/template_iter%d_part%d_%d_Cluster%d.png',0,PartLocX(iPart),PartLocY(iPart), c);
            % add a small margin to make sure all Gabor elements are displayed fully
            im = displayMatchedTemplate([partSizeX+2*halfFilterSize partSizeY+2*halfFilterSize],clusters{c}.selectedx{iPart}+halfFilterSize,...
                    clusters{c}.selectedy{iPart}+halfFilterSize,clusters{c}.selectedOrient{iPart},zeros(length(ind),1,'single'),clusters{c}.selectedlambda{iPart},allSymbol,numOrient);
            towrite = -double(im);
            towrite = uint8(255 * (towrite-min(towrite(:)))/(1e-10+max(towrite(:))-min(towrite(:))));
            if max(towrite(:)) == 0
                towrite(:) = 255;
            end
            %imwrite(towrite,filename);
        end

        % split the object template into larger overlapping partial templates (context sensitive)
        for iPart = 1:numCandPart
            ind = find( tmpSelectedx >= PartLocX(iPart)-partMarginX & tmpSelectedx < PartLocX(iPart) + partSizeX + partMarginX & ...
                    tmpSelectedy >= PartLocY(iPart) - partMarginX & tmpSelectedy < PartLocY(iPart) + partSizeY + partMarginY );
            clusters{c}.largerSelectedOrient{iPart} = single( tmpSelectedo(ind) );
            clusters{c}.largerSelectedx{iPart} = single( tmpSelectedx(ind) - PartLocX(iPart) + partMarginX );
            clusters{c}.largerSelectedy{iPart} = single( tmpSelectedy(ind) - PartLocY(iPart) + partMarginY );
            clusters{c}.largerSelectedlambda{iPart} = single( tmpSelectedlambda(ind) );
            clusters{c}.largerSelectedLogZ{iPart} = single( tmpSelectedlogz(ind) );
            filename = sprintf('output/largertemplate_iter%d_part%d_%d_Cluster%d.png',0,PartLocX(iPart),PartLocY(iPart), c);
            im = displayMatchedTemplate([partSizeX+2*partMarginX+2*halfFilterSize partSizeY+2*partMarginY+halfFilterSize*2],clusters{c}.largerSelectedx{iPart}+halfFilterSize,...
                    clusters{c}.largerSelectedy{iPart}+halfFilterSize,clusters{c}.largerSelectedOrient{iPart},zeros(length(ind),1,'single'),clusters{c}.largerSelectedlambda{iPart},allSymbol,numOrient);
            towrite = -double(im);
            towrite = uint8(255 * (towrite-min(towrite(:)))/(1e-10+max(towrite(:))-min(towrite(:))));
            if max(towrite(:)) == 0
                towrite(:) = 255;
            end
            %imwrite(towrite,filename);
        end

        RotateTemplate;

        disp(['mex-C learning time: ' num2str(toc) ' seconds']);


        clusters{c}.PartOnOff = ones(numCandPart,1); % all parts are selected initially
        clusters{c}.S3SelectedRow = zeros(1,numCandPart,'single');
        clusters{c}.S3SelectedCol = zeros(1,numCandPart,'single');
        clusters{c}.S3SelectedOri = zeros(1,numCandPart,'single');
        for iPart = 1:numCandPart
            clusters{c}.S3SelectedRow(iPart) = PartLocX(iPart) - 1 + floor(partSizeX/2);
            clusters{c}.S3SelectedCol(iPart) = PartLocY(iPart) - 1 + floor(partSizeY/2);
        end
        clusters{c}.allS3SelectedRow = zeros(numRotate,numCandPart,'single');
        clusters{c}.allS3SelectedCol = zeros(numRotate,numCandPart,'single');
        clusters{c}.allS3SelectedOri = zeros(numRotate,numCandPart,'single');
        RotateS3Template;

    end

    %% EM Learning 

    patches = cell(numImage,1);

    SUM1mapLearn = cell(numImage,numOrient);
    for i = 1:numel(SUM1mapLearn)
        SUM1mapLearn{i} = zeros(sizeTemplatex,sizeTemplatey,'single');
    end

    %{
    SUM1mapLearnPerPart = cell(numCandPart,numImage,numOrient);
    for i = 1:numel(SUM1mapLearnPerPart)
        SUM1mapLearnPerPart{i} = zeros(partSizeX,partSizeY,'single');
    end
    %}

    for it = 1:numIteration

        %% Iteration part 1 (E step): detect the object in each image using SUM3 map

        disp(['detection for step ' num2str(it)]);
        aveMAX2 = zeros(numCandPart,1); % records the average MAX2 scores for each non-overlapping partial templates

        MAX3scoreAll(:,:) = -1e10;
        
        for img=1:numImage

            disp(['    start detecting in image ' num2str(img)]); tic
            load( sprintf('working/ImageAndFeature_%d.mat',img) );
            for c = 1:NumCluster 
            
                selectedPart = find(clusters{c}.PartOnOff);
                
                %% compute SUM2 maps for non-overlapping parts
                SUM2map = cell(numPartRotate,numCandPart,numResolution);
                S2T = cell( numPartRotate, numCandPart );
                for iPart = 1:numCandPart
                    for r = 1:length(partRotationRange)
                        S2T{r,iPart} = struct( 'selectedRow',single(clusters{c}.allSelectedx{iPart,r}(:)),...
                            'selectedCol', single(clusters{c}.allSelectedy{iPart,r}(:)),...
                            'selectedOri', single(clusters{c}.allSelectedOrient{iPart,r}(:)),...
                            'selectedScale', zeros(length(clusters{c}.allSelectedx{iPart,r}),1,'single'),...
                            'selectedLambda', single(clusters{c}.selectedlambda{iPart}(:)),...
                            'selectedLogZ', single(clusters{c}.selectedLogZ{iPart}(:)) );
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
                        largerS2T{r,iPart} = struct( 'selectedRow',single(clusters{c}.largerAllSelectedx{iPart,r}(:)),...
                            'selectedCol', single(clusters{c}.largerAllSelectedy{iPart,r}(:)),...
                            'selectedOri', clusters{c}.largerAllSelectedOrient{iPart,r}(:),...
                            'selectedScale', zeros(length(clusters{c}.largerAllSelectedx{iPart,r}),1,'single'),...
                            'selectedLambda', single(clusters{c}.largerSelectedlambda{iPart}(:)),...
                            'selectedLogZ', single(clusters{c}.largerSelectedLogZ{iPart}(:)));
                    end
                end
                
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
                        selectedTransform(j) = find( clusters{c}.allS3SelectedOri(r,j) == partRotationRange );
                    end
                    S3T{r} = struct( 'selectedRow',single(floor(.5+clusters{c}.allS3SelectedRow(r,selectedPart)/subsampleM2/subsampleS2)),...
                            'selectedCol', single(floor(.5+clusters{c}.allS3SelectedCol(r,selectedPart)/subsampleM2/subsampleS2)),...
                            'selectedInd', single(selectedPart) - 1,...
                            'selectedTransform', selectedTransform - 1,...
                            'selectedLambda', ones(length(selectedPart),1,'single'),...
                            'selectedLogZ', single( 0*ones(length(selectedPart),1) ) );
                    SUM3map = cell(numResolution,1);
                    for iRes = 1:numResolution
                        tmpM2 = MAX2map(:,:,iRes);
                        SUM3map(iRes) = mexc_ComputeSUM3( tmpM2(:), S3T(r), 1, numPartRotate );
                        [MAX3 loc] = max(SUM3map{iRes}(:));
                        if MAX3 > MAX3scoreAll(img, c)
                            MAX3scoreAll(img, c) = MAX3;
                            bestS3Loc = loc;
                            bestRes = iRes; % best object resolution
                            bestRot = rot; % starting from a negative number
                        end
                    end

                end

                % argmax location of SUM3map
                if img == starting
                    bestRot = 0;
                    bestRotInd = find(bestRot==rotationRange);
                    therex = floor( (startx + endx) / 2 / subsampleS2 );
                    therey = floor( (starty + endy) / 2 / subsampleS2 );
                    bestRes = originalResolution;
                else
                    bestRotInd = find(bestRot==rotationRange);
                    therey = ceil(bestS3Loc/size(SUM3map{bestRes},1));
                    therex = bestS3Loc - (therey-1) * size(SUM3map{bestRes},1);
                end
                
                % For detection evaluation:
                
                % if positive image
                %  posScores(img) = MMAX3;
                
                % if large natural image
                % use SUM3map
                %for j = 1:length(SUM3map)
                %	tmp = SUM3map{j}(1:4:end,1:4:end);
                %	negScores = [negScore;tmp(:)];
                %end

                if MAX3scoreAll(img,c)==max(MAX3scoreAll(img,:)) 
                    % copy the detected patch
                    for iPart = 1:numCandPart
                        r = find( clusters{c}.allS3SelectedOri(bestRotInd,iPart) == partRotationRange ); % the index of part rotation
                        Fx = therex + floor(.5+clusters{c}.allS3SelectedRow(bestRotInd,iPart)/subsampleM2/subsampleS2);
                        Fy = therey + floor(.5+clusters{c}.allS3SelectedCol(bestRotInd,iPart)/subsampleM2/subsampleS2); % sub-sampled position
                        imagesize = size(MAX2map{r,iPart,bestRes}); % subsampled image size
                        if Fx >= 1 && Fx <= imagesize(1) && Fy >= 1 && Fy <= imagesize(2)
                            tmp = MAX2map{r,iPart,bestRes};
                            aveMAX2(iPart) = aveMAX2(iPart) + tmp(Fx,Fy) / numImage;
                            
                            tmp = MAX2ResolutionTrace{r,iPart,bestRes};
                            bestPartRes = tmp(Fx,Fy) + 1; % best part resolution
                            current_size = size(tmp);
                            
                            tmp = largerMAX2LocTrace{r,iPart,bestPartRes};
                            new_size = size(tmp);
                            Fx = floor(.5+Fx*double(new_size)/current_size);
                            Fy = floor(.5+Fy*double(new_size)/current_size);
                            
                            if Fx >= 1 && Fx <= size(tmp,1) && Fy >= 1 && Fy <= size(tmp,2)
                                translationInd = tmp(Fx,Fy) + 1;
                            else
                                translationInd = floor(size(M2RowColShift,1)/2);
                            end
                            
                            tmp = largerMAX2TransformTrace{r,iPart,bestPartRes};
                            if Fx >= 1 && Fx <= size(tmp,1) && Fy >= 1 && Fy <= size(tmp,2)
                                transformInd = tmp(Fx,Fy) + 1;
                            else
                                transformInd = floor(numPartRotate/2) + 1;
                            end
                            
                            actualPartRotationInd = transformInd - numPartRotate*(ceil(double(transformInd)/numPartRotate)-1);
                            Fx = floor( Fx + M2RowColShift(translationInd,1) * sqrt(partSizeX*partSizeY)/subsampleS2 );
                            Fy = floor( Fy + M2RowColShift(translationInd,2) * sqrt(partSizeX*partSizeY)/subsampleS2 );
                            actualPartRotation = partRotationRange(actualPartRotationInd);
                        else
                            actualPartRotationInd = r;
                            actualPartRotationInd = find(r==partRotationRange);
                        end


                        denseX = -floor(partSizeX/2) + (1:partSizeX);
                        denseY = -floor(partSizeY/2) + (1:partSizeY);
                        count = 0;
                        inRow = zeros(length(denseX)*length(denseY),1,'single');
                        inCol = zeros(length(denseX)*length(denseY),1,'single');
                        for y = denseY
                            for x = denseX
                                count = count+1;
                                inRow(count) = x;
                                inCol(count) = y;
                            end
                        end
                        tScale = 0; rScale = 1; cScale = 1; inO = zeros(numel(inRow),1,'single'); inS = zeros(numel(inRow),1,'single');
                        [outRow, outCol] = ...
                            mexc_TemplateAffineTransform(tScale,rScale,cScale,...
                                actualPartRotation,inRow,inCol,inO,inS,numOrient);

                        % find the part location at the higher resolution
                        Fx = (Fx-1 + .5) * subsampleS2 * subsampleM2;
                        Fy = (Fy-1 + .5) * subsampleS2 * subsampleM2;
                        
                        % crop the feature patch that is registered to the part template
                        tmpSUM1mapLearn = mexc_CropInstance(SUM1mapFind(bestPartRes,:),Fx,Fy,...
                            actualPartRotation,tScale,1,...
                            outRow,outCol,...
                            numOrient,1,partSizeX,partSizeY);

                        for o = 1:numOrient
                            SUM1mapLearn{img,o}(PartLocX(iPart)-1+(1:partSizeX),PartLocY(iPart)-1+(1:partSizeY)) = tmpSUM1mapLearn{o};
                        end
                    
                    end
                    ImageCluster(img) = c;
                end
            end

            % overlay the template
            
            imageSizeAtBestObjectResolution = size( ImageMultiResolution{bestRes} );
            matchedSym = zeros( imageSizeAtBestObjectResolution );
            for iPart = 1:numCandPart
                
                % render the template for each part separately, then overlay the rendered images
                gaborXX = [];
                gaborYY = [];
                gaborOO = [];
                gaborMM = [];
            
                r = find( clusters{ImageCluster(img)}.allS3SelectedOri(bestRotInd,iPart) == partRotationRange ); % the index of part rotation
                % part location
                Fx = therex + floor(.5+clusters{ImageCluster(img)}.allS3SelectedRow(bestRotInd,iPart)/subsampleS2);
                Fy = therey + floor(.5+clusters{ImageCluster(img)}.allS3SelectedCol(bestRotInd,iPart)/subsampleS2); % sub-sampled position
                imagesize = size(largerMAX2LocTrace{r,iPart,bestRes});
                if Fx > 0 && Fx <= imagesize(1) && Fy > 0 && Fy <= imagesize(2)
                    tmp = MAX2ResolutionTrace{r,iPart,bestRes};
                    bestPartRes = tmp(Fx,Fy) + 1; % best part resolution
                    current_size = size(tmp);
                    
                    tmp = largerMAX2LocTrace{r,iPart,bestPartRes};
                    new_size = size(tmp);
                    Fx = floor(.5+Fx*double(new_size)/current_size);
                    Fy = floor(.5+Fy*double(new_size)/current_size);
                    
                    if Fx >= 1 && Fx <= size(tmp,1) && Fy >= 1 && Fy <= size(tmp,2);
                        translationInd = tmp(Fx,Fy) + 1;
                    else
                        translationInd = floor(size(M2RowColShift,1)/2);
                    end
                    
                    tmp = largerMAX2TransformTrace{r,iPart,bestPartRes};
                    if Fx >= 1 && Fx <= size(tmp,1) && Fy >= 1 && Fy <= size(tmp,2);
                        transformInd = tmp(Fx,Fy) + 1;
                    else
                        transformInd = floor(numPartRotate/2) + 1;
                    end
                    
                    actualPartRotationInd = transformInd - numPartRotate*(ceil(double(transformInd)/numPartRotate)-1);
                    Fx = floor( Fx + M2RowColShift(translationInd,1) * sqrt(partSizeX*partSizeY)/subsampleS2 );
                    Fy = floor( Fy + M2RowColShift(translationInd,2) * sqrt(partSizeX*partSizeY)/subsampleS2 );
                end
                
                % Gabor basis elements locations
                for j = 1:length( clusters{ImageCluster(img)}.selectedx{iPart} )
                    gaborX = floor(.5 + (Fx+.5)*subsampleS2 + clusters{ImageCluster(img)}.allSelectedx{iPart,actualPartRotationInd}(j));
                    gaborY = floor(.5 + (Fy+.5)*subsampleS2 + clusters{ImageCluster(img)}.allSelectedy{iPart,actualPartRotationInd}(j));
                    gaborO = clusters{ImageCluster(img)}.allSelectedOrient{iPart,actualPartRotationInd}(j);
                    if gaborX > 0 && gaborX <= size(M1Trace{bestPartRes,1},1) && gaborY > 0 && gaborY <= size(M1Trace{bestPartRes,1},2)
                        trace = M1Trace{bestPartRes,gaborO+1}(gaborX,gaborY) + 1;
                        dx = M1RowShift{gaborO+1}(trace);
                        dy = M1ColShift{gaborO+1}(trace);
                        shiftedo = M1OriShifted{gaborO+1}(trace);
                        gaborX = floor(.5 + gaborX + single(dx));
                        gaborY = floor(.5 + gaborY + single(dy));
                        gaborO = single(shiftedo);
                    end
                    gaborXX = [gaborXX;gaborX];
                    gaborYY = [gaborYY;gaborY];
                    gaborOO = [gaborOO;gaborO];
                    if gaborX > 0 && gaborX <= size(M1Trace{bestPartRes,1},1) && gaborY > 0 && gaborY <= size(M1Trace{bestPartRes,1},2)
                        val = sqrt( SUM1mapFind{bestPartRes,gaborO+1}(gaborX,gaborY) );
                        val = max(0, val-0.3);
                    else
                        val = 0;
                    end
                    gaborMM = [gaborMM; val];
                end
                tmpMatchedSym = displayMatchedTemplate(size(ImageMultiResolution{bestPartRes}),gaborXX,...
                    gaborYY,gaborOO,zeros(length(gaborXX),1,'single'),gaborMM,allSymbol,numOrient);
                tmpMatchedSym = double( imresize(tmpMatchedSym,imageSizeAtBestObjectResolution,'bilinear') );
                matchedSym = max(matchedSym,tmpMatchedSym);
            end
            
            srcImg = ImageMultiResolution{bestRes};
            

            filename = sprintf('output/matchedS3T%d_iter%d.png',img,it);
            towrite = -double(matchedSym);
            towrite = uint8(255 * (towrite-min(towrite(:)))/(1e-10+max(towrite(:))-min(towrite(:))));
            imwrite(towrite,filename);
            filename = sprintf('output/matchedS3image%d.png',img);
            towrite = double(srcImg);
            towrite = uint8(255 * (towrite-min(towrite(:)))/(1e-10+max(towrite(:))-min(towrite(:))));
            imwrite(towrite,filename);


            denseX = -floor(sizeTemplatex/2) + (1:sizeTemplatex);
            denseY = -floor(sizeTemplatey/2) + (1:sizeTemplatey);
            count = 0;
            inRow = zeros(length(denseX)*length(denseY),1,'single');
            inCol = zeros(length(denseX)*length(denseY),1,'single');
            for y = denseY
                for x = denseX
                    count = count+1;
                    inRow(count) = x;
                    inCol(count) = y;
                end
            end
            tScale = 0; rScale = 1; cScale = 1; inO = zeros(numel(inRow),1,'single'); inS = zeros(numel(inRow),1,'single');
            [outRow, outCol] = ...
                mexc_TemplateAffineTransform(tScale,rScale,cScale,...
                bestRot,inRow,inCol,inO,inS,numOrient);
            
            patches(img) = mexc_CropInstance(ImageMultiResolution(bestRes),floor((therex+.5)*subsampleS2),floor((therey+.5)*subsampleS2),...
                bestRot,0,1,...
                outRow,outCol,...
                1,1,sizeTemplatex,sizeTemplatey);
        end
        towrite = displayImages(patches,5,60,60);

        for c = 1:NumCluster
            [~, ImageCluster] = max(MAX3scoreAll, [], 2);
            clusters{c}.imageIdx = find(ImageCluster==c);
            clusters{c}.numImage = size(clusters{c}.imageIdx, 1);
        end
        
        %% Iteration part 2 (M step): re-learn the template
        
        for c = 1:NumCluster
            
            disp(['multiple image learning for step ' num2str(it) ' on Cluster' num2str(c)]);
            tic
            
            tmpSelectedx = zeros(numElement,1);
            tmpSelectedy = zeros(numElement,1);
            tmpSelectedo = zeros(numElement,1);
            tmpSelectedlambda = zeros(numElement,1);
            tmpSelectedlogz = zeros(numElement,1);
            deformedTemplate = cell(numImage,1);
            tmpCommonTemplate = zeros( sizeTemplatex, sizeTemplatey, 'single' );
            for i = 1:numImage
                deformedTemplate{i} = zeros( sizeTemplatex, sizeTemplatey, 'single' );
            end

            CsharedSketch(numOrient, locationShiftLimit, orientShiftLimit, subsample, ... % about active basis  
            numElement, clusters{c}.numImage, sizeTemplatex, sizeTemplatey, SUM1mapLearn(clusters{c}.imageIdx, :), ... % about training images 
            halfFilterSize, Correlation, allSymbol(1, :), ... % about filters
            numStoredPoint, storedlambda, storedExpectation, storedLogZ, ... % about exponential model 
            tmpSelectedo, tmpSelectedx, tmpSelectedy, tmpSelectedlambda, tmpSelectedlogz, ... % learned parameters
            tmpCommonTemplate, deformedTemplate); % learned templates
            
            filename = sprintf('output/%d/template_iter%d.png',c,it);
            towrite = -double(tmpCommonTemplate);
            towrite = uint8(255 * (towrite-min(towrite(:)))/(1e-10+max(towrite(:))-min(towrite(:))));
            if max(towrite(:)) == 0
                towrite(:) = 255;
            end
            imwrite(towrite,filename);
            
            % split the object template into non-overlapping partial templates
            for iPart = 1:numCandPart
                ind = find( tmpSelectedx >= PartLocX(iPart) & tmpSelectedx < PartLocX(iPart) + partSizeX & ...
                        tmpSelectedy >= PartLocY(iPart) & tmpSelectedy < PartLocY(iPart) + partSizeY );
                clusters{c}.selectedOrient{iPart} = single( tmpSelectedo(ind) );
                clusters{c}.selectedx{iPart} = single( floor( tmpSelectedx(ind) - PartLocX(iPart) ) );
                clusters{c}.selectedy{iPart} = single( floor( tmpSelectedy(ind) - PartLocY(iPart) ) );
                clusters{c}.selectedlambda{iPart} = single( tmpSelectedlambda(ind) );
                clusters{c}.selectedLogZ{iPart} = single( tmpSelectedlogz(ind) );
                filename = sprintf('output/%d/template_iter%d_part%d_%d.png',c,it,PartLocX(iPart),PartLocY(iPart));
                % add a small margin to make sure all Gabor elements are displayed fully
                im = displayMatchedTemplate([partSizeX+2*halfFilterSize partSizeY+2*halfFilterSize],clusters{c}.selectedx{iPart}+halfFilterSize,...
                    clusters{c}.selectedy{iPart}+halfFilterSize,clusters{c}.selectedOrient{iPart},zeros(length(ind),1,'single'),clusters{c}.selectedlambda{iPart},allSymbol,numOrient);
                towrite = -double(im);
                towrite = uint8(255 * (towrite-min(towrite(:)))/(1e-10+max(towrite(:))-min(towrite(:))));
                if max(towrite(:)) == 0
                    towrite(:) = 255;
                end
                %imwrite(towrite,filename);
            end

            % split the object template into larger overlapping partial templates (context sensitive)
            for iPart = 1:numCandPart
                ind = find( tmpSelectedx >= PartLocX(iPart)-partMarginX & tmpSelectedx < PartLocX(iPart) + partSizeX + partMarginX & ...
                        tmpSelectedy >= PartLocY(iPart) - partMarginY & tmpSelectedy < PartLocY(iPart) + partSizeY + partMarginY );
                clusters{c}.largerSelectedOrient{iPart} = single( tmpSelectedo(ind) );
                clusters{c}.largerSelectedx{iPart} = single( tmpSelectedx(ind) - PartLocX(iPart) + partMarginX );
                clusters{c}.largerSelectedy{iPart} = single( tmpSelectedy(ind) - PartLocY(iPart) + partMarginY );
                clusters{c}.largerSelectedlambda{iPart} = single( tmpSelectedlambda(ind) );
                clusters{c}.largerSelectedLogZ{iPart} = single( tmpSelectedlogz(ind) );
                filename = sprintf('output/%d/largertemplate_iter%d_part%d_%d.png',c,it,PartLocX(iPart),PartLocY(iPart));
                im = displayMatchedTemplate([partSizeX+2*partMarginX+2*halfFilterSize partSizeY+2*partMarginY+halfFilterSize*2],clusters{c}.largerSelectedx{iPart}+halfFilterSize,...
                    clusters{c}.largerSelectedy{iPart}+halfFilterSize,clusters{c}.largerSelectedOrient{iPart},zeros(length(ind),1,'single'),clusters{c}.largerSelectedlambda{iPart},allSymbol,numOrient);
                towrite = -double(im);
                towrite = uint8(255 * (towrite-min(towrite(:)))/(1e-10+max(towrite(:))-min(towrite(:))));
                if max(towrite(:)) == 0
                    towrite(:) = 255;
                end
                %imwrite(towrite,filename);
            end
            
            RotateTemplate;
        
        end
        save(sprintf('template/hab_task%d_iter%d.mat', task_id, it), 'clusters');
    end
    tic_toc(task_id) = toc(big_ticID);
end

disp('done.');

