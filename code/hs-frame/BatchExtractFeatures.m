function [SPM_SUM2, SPM_SUM3] = BatchExtractFeatures(imgCell,codebook, para, idx)

    %%%%% parameters for multiple selection (cell size) or deformable parts 
    
    spm_numLayers = 3;
    spm_threshold = 0;
    read_from_data = true;
    
    %%%%%%%%%%%%%%%%%%%
    nPartCol=para.nPartCol;
    nPartRow=para.nPartRow;
    part_sx=para.part_sx;
    part_sy=para.part_sy;
    gradient_threshold_scale=0.8;     %  we use addaptive threshold in the multiple selection.  The threshold is equal to (gradient_threshold_scale * maximum gradient).
    numPart=nPartCol*nPartRow; % number of deformable parts
    argmaxMethod=1;  % 1: local max in squared region, 0: cross-shaped rgion for parts
    relativePartRotationRange = para.relativePartRotationRange;  % relative rotation range for parts
    relativePartLocationRange= para.relativePartLocationRange; % default 1
    resolutionShiftLimit = para.resolutionShiftLimit;
    sx=nPartRow*part_sx;   % template size x
    sy=nPartCol*part_sy;   % template size y
    
    
    rotateShiftLimit = 5;%para.rotateShiftLimit; 4 % template rotation  from -rotateShiftLimit to rotateShiftLimit, eg. (1)-2:2 if rotateShiftLimit=2 (2)0 is without rotation
    rotationRange = -rotateShiftLimit:rotateShiftLimit;
    numRotate = length(rotationRange);
    RatioDisplacementSUM3 = para.ratioDisplacementSUM3;   % default=0. Compute all values in SUM3 map
    
    nOrient = 16;
    locationShiftLimit=para.locationShiftLimit;  % arg-max
    orientShiftLimit=para.orientShiftLimit;    % arg-max
    %%%%%%
    
    %%%%%%% parameters about resolution
    numResolution = para.numResolution; % number of resolutions to search for in detection stage, the middle level is the original size
   
    
    
    scaleStepSize = 0.1; % you can tune this either 0.1 or 0.2
    originalResolution = round(numResolution/2);    % original resolution is the one at which the imresize factor = 1.  the middle one is by default
    displace = 0;
    originalResolution = originalResolution - displace;    % shift the index of the original resolution to the left.  We don't want too many shrinkage resolution.

    %%% local normalization parameters
    isLocalNormalize= para.isLocalNormalize; %false; %
    isSeparate=false;
    localNormScaleFactor=2; %0.5, 1, 2, 3, or 4, we suggest 2;
    thresholdFactor=0.01;

    GaborScaleList = para.GaborScaleList;
    DoGScaleList = para.DoGScaleList;

    if isLocalNormalize
        DoGScaleList=[];
    end

    nScaleDoG=length(DoGScaleList);
    nScaleGabor=length(GaborScaleList);

    %%%%%%%%%%%%%%% auxiliary variables for hierachical template
    partRotations=cell(numRotate, 1);
    partRotationRange = [];   % save the possible range of part rotation
    for ii = 1:numRotate
        tmpRotation = rotationRange(ii)+relativePartRotationRange;
        partRotations{ii} = tmpRotation;
        partRotationRange = union( partRotationRange, tmpRotation );
    end
    numPartRotate = length(partRotationRange);

    % record the top-left location of each part in a large template
    PartLocX0 = 1:part_sx:sx-part_sx+1;
    PartLocY0 = 1:part_sy:sy-part_sy+1;

    % another representation of the top-left location
    PartLocX = zeros(numPart,1);
    PartLocY = zeros(numPart,1);
    iPart = 1;
    for x = PartLocX0
        for y = PartLocY0
            PartLocX(iPart) = x;
            PartLocY(iPart) = y;
            iPart = iPart + 1;
        end
    end
    
    % compute affinity matrix
    minRotationDif = (sin(1*pi/nOrient)-sin(0))^2 + (cos(1*pi/nOrient)-cos(0))^2 + 1e-10;  % 1e-10

    allTemplateAffinityMatrix = cell(numPartRotate, numPart);
    templateAffinityMatrix = cell(numPartRotate,numPart);
    for iPart = 1:numPart
        for r1 = 1:length(partRotationRange)
            angle1 = pi/nOrient * partRotationRange(r1);
            templateAffinityMatrix{r1,iPart} = [];
            jPart = iPart;
            for r2 = 1:length(partRotationRange)
                angle2 = pi/nOrient*partRotationRange(r2);
                if (sin(angle1) - sin(angle2))^2 + (cos(angle1)-cos(angle2))^2 <= minRotationDif
                    templateAffinityMatrix{r1,iPart} = int32( [templateAffinityMatrix{r1,iPart} r2-1] );
                end
            end
        end
    end
    for iPart = 1:numPart
        startInd = (iPart-1)*numPartRotate;
        for r2 = 1:length(partRotationRange)
            allTemplateAffinityMatrix{r2, iPart} = templateAffinityMatrix{r2, iPart}(:)+startInd;
        end
    end
    %%% create filter bank
    filters=[];
    for iScale=1:nScaleGabor
        f = MakeFilter(GaborScaleList(iScale),nOrient);
        for i=1:nOrient
            f_r{i} =single(real(f{i}));
            f_i{i} =single(imag(f{i}));
        end

        filters = [filters f_r f_i];
    end

    for iScale=1:nScaleDoG

        f0 = single(dog(DoGScaleList(iScale),0));
        filters = [filters f0];

    end

    numFilter = length(filters);
    halfFilterSizes = zeros(size(filters));
    for iF = 1:numFilter
        halfFilterSizes(iF)=(size(filters{iF},1)-1)/2;
    end
    overAllArea = sum((sx-2*halfFilterSizes).*(sy-2*halfFilterSizes));
    Corr = CorrFilterFrame(filters); %calculate correlation among filters

    
    %%
    i=0;
    for iClass = 1:para.numCategory
        for c=1:para.numCluster

            i=i+1;
            cluster_single=codebook(i,:);  

            [cluster_single.S2T, cluster_single.S3T] = hierachicalTemplate(numPart, part_sx, part_sy, sx, sy, rotateShiftLimit, nOrient, numRotate, cluster_single.template, nScaleGabor, partRotationRange,PartLocX, PartLocY);
            codebook(i,:) = cluster_single;
        end
    end
    %%
    
    SPM_SUM2 = [];
    SPM_SUM3 = [];

    for iImg = 1:length(imgCell)
        
        if read_from_data && exist(['./output/map' num2str(idx(iImg)) '.mat'],'file')
            load(['./output/map' num2str(idx(iImg)) '.mat']);
            disp(['====> Found cache data for image ' num2str(iImg)]);
        else

            disp(['======> start filtering and maxing image ' num2str(iImg)]);

            img = imgCell{iImg};
            if ndims(img)==3
                img=rgb2gray(img);
            end
            img = im2single(img);
            allSizex = zeros(1, numResolution);
            allSizey = zeros(1, numResolution);
            ImageMultiResolution = cell(1, numResolution);
            for resolution=1:numResolution
                resizeFactor = 1.0 + (resolution - originalResolution)*scaleStepSize;
                img2 = imresize(img, resizeFactor, 'nearest');  % images at multiple resolutions
                img2 = img2-mean(img2(:));
                img2 = img2/std(img2(:))*sqrt(10);
                ImageMultiResolution{resolution} = img2;
                [sizex, sizey] = size(ImageMultiResolution{resolution});
                allSizex(resolution) = sizex;
                allSizey(resolution) = sizey;
            end

            % compute MAX1 map for images in different resolutions
            tic
            [SUM1mapFind, MAX1mapFind] = applyfilterBank_MultiResolution_sparseV4(ImageMultiResolution, filters, halfFilterSizes, nOrient,...
                locationShiftLimit,orientShiftLimit,isLocalNormalize,isSeparate,localNormScaleFactor,thresholdFactor,nScaleGabor,nScaleDoG, sqrt(10));  % if using local normalization, the last parameter is important.
            numPartRotate = length(partRotationRange);

            SUM3map = cell(numRotate, numResolution, length(codebook));
            SUM2map = cell(numPartRotate, numPart, numResolution, length(codebook));
            SUM2map_temp = cell(numPartRotate, numPart, numResolution);

            for c = 1:length(codebook)
                % compute SUM2 map

                current_size = 0;
                for iRes = 1:numResolution
                    tmpSUM2=sparseFRAME_SUM2_part(single(allSizex(iRes)), single(allSizey(iRes)), single(numFilter), codebook(c).S2T(:), MAX1mapFind(iRes,:));
                    if current_size==0, current_size = size(tmpSUM2{1});end
                    for npr = 1:length(tmpSUM2)
                            tmpSUM2{npr} = imresize(tmpSUM2{npr},current_size,'nearest');
                    end
                    SUM2map_temp(:,:,iRes)=reshape(tmpSUM2, [numPartRotate, numPart]);
                    SUM2map(:,:,iRes,c)=reshape(tmpSUM2, [numPartRotate, numPart]);
                end
                % compute MAX2 map
                tmpMAX2map=cell(size(SUM2map_temp));
                MAX2LocTrace=cell(size(SUM2map_temp));
                MAX2TransformTrace=cell(size(SUM2map_temp));

                for iRes = 1:numResolution
                    [tmpMAX2, tmpMAX2LocTrace, tmpMAX2TransformTrace,M2RowColShift] = ...
                        mexc_sparseFRAME_MAX2( allTemplateAffinityMatrix(:), SUM2map_temp(:,:,iRes), relativePartLocationRange, argmaxMethod, 1);

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
                                    ind = ref > map;
                                    map(ind) = ref(ind);
                                    MAX2map{j,k,iRes}(ind) = ref(ind);
                                    MAX2ResolutionTrace{j,k,iRes}(ind) = jRes - 1; % start from 0
                                end
                            end
                        end
                    end
                end
                %%%% align the object by detection
                %% compute SUM3map
                for r = 1:numRotate % this is the rotation of the S3 template
                    for iRes = 1:numResolution
                %         tmpM2 = MAX2map(:,:,iRes);
                %         tmpS3 = mexc_ComputeSUM3( tmpM2(:), codebook(c).S3T(r), 1, numPartRotate);
                %         SUM3map(r, iRes) = {tmpS3{1} - codebook(c).logZ};  % for learning mixture model, we need log of Z.
                        tmpM2 = MAX2map(:,:,iRes);      
                    % SUM3map(r, iRes) = mexc_ComputeSUM3_logZ( tmpM2(:), codebook(c).S3T(r), 1, numPartRotate, single(codebook(c).logZ));
                        SUM3map(r, iRes, c) = mexc_ComputeSUM3_logZ_partial( tmpM2(:), codebook(c).S3T(r), 1, numPartRotate, single(codebook(c).logZ), single(RatioDisplacementSUM3));
                    end
                end
            end % codebook

            SUM2map = reshape(SUM2map, [numResolution*numPartRotate*numPart*length(codebook), 1]);
            SUM3map = reshape(SUM3map, [numResolution*numRotate*length(codebook), 1]);
        end
        SPM_SUM2 = [SPM_SUM2 ; featureSPM(SUM2map,spm_numLayers,spm_threshold)];
        SPM_SUM3 = [SPM_SUM3 ; featureSPM(SUM3map,spm_numLayers,spm_threshold)];
        
        %save(['./output/map' num2str(idx(iImg)) '.mat'], 'SUM2map', 'SUM3map');
        disp(['finished MAX3 map of image ' num2str(iImg) ' : ' num2str(toc) ' seconds']);
    end % pic
        
end % end of function


