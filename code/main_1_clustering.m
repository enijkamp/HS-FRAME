% clear
clear
close all;

% code
addpath('hs-frame');

% config
para = config();

% mex
compileMex();

% run
task_ids = para.task_ids;
noWorkers = para.noWorkers;
tic_toc = zeros(length(task_ids), 5);

for task_id = task_ids

    categoryName = para.categoryNames{task_id};
    
    for seed = 1:1

        big_ticID = tic;
        
        disp(['we will start the clustering task ' num2str(task_id) ' : ' categoryName]);
        
        %%%% CPU paralleled for computing SUM2, MAX2, SUM3, and MAX3
        
        runParallel=1;  % Set variable to know when to run parallel
        
        % set seed for random numbers generation
        rng(seed);
        
        % parameters about cropping for relearning
        %toUseCroppedImg=false;
        
        %%%% EM clustering parameters
        numEMIteration = para.numEMIteration;
        numCluster=para.numCluster;
        %%%%%
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% paraemters for two-stage learning algorithm
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        sigsq = 10;
        numSketch = para.numWavelet; %400 sparsification
        numTopFeatureToShow = 70; % not used so far, please ignore
        lambdaLearningRate_MP = 0.1/sqrt(sigsq);  % 0.1
        epsilon = 0.03;
        L = 10;
        nIter = 80;
        numSample = 3; % how many HMC calls for each learning iteration
        isSaved = 1;
        
        %%%%% parameters for multiple selection (cell size) or deformable parts %%%%%%%%%%%%%%%%%%%
        nPartCol = para.nPartCol;
        nPartRow = para.nPartRow;
        part_sx = para.part_sx;
        part_sy = para.part_sy;
        
        gradient_threshold_scale = 0.8; %  we use adaptive threshold in the multiple selection. The threshold is equal to (gradient_threshold_scale * maximum gradient).
        numPart = nPartCol*nPartRow; % number of deformable parts
        
        argmaxMethod = 1;  % 1: local max in squared region, 0: cross-shaped rgion for parts
        relativePartRotationRange = para.relativePartRotationRange;  % relative rotation range for parts
        relativePartLocationRange= para.relativePartLocationRange; % default 1
        resolutionShiftLimit = para.resolutionShiftLimit;
        %%%%%
        
        %%%%% parameters for large template
        sx = nPartRow*part_sx;   % template size x
        sy = nPartCol*part_sy;   % template size y
        rotateShiftLimit = para.rotateShiftLimit;   % template rotation from -rotateShiftLimit to rotateShiftLimit, eg. (1)-2:2 if rotateShiftLimit=2 (2)0 is without rotation
        rotationRange = -rotateShiftLimit:rotateShiftLimit;
        numRotate = length(rotationRange);
        RatioDisplacementSUM3 = para.ratioDisplacementSUM3;   % default=0. Compute all values in SUM3 map
        
        nOrient = 16;
        locationShiftLimit = para.locationShiftLimit;  % arg-max
        orientShiftLimit = para.orientShiftLimit;    % arg-max
        %%%%%%
        
        %%%%%%% parameters about resolution
        numResolution = para.numResolution;  % number of resolutions to search for in detection stage, the middle level is the original size
        scaleStepSize = 0.1; % you can tune this either 0.1 or 0.2
        originalResolution = round(numResolution/2);    % original resolution is the one at which the imresize factor = 1.  the middle one is by default
        displace = 0;
        originalResolution = originalResolution - displace;    % shift the index of the original resolution to the left.  We don't want too many shrinkage resolution.
        %%%%%%%
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        interval = 1; % 5 control how many iterations we need to wait to select next wavlet
        numWavelet = para.numWavelet;   % default 300, 370
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        
        %%%%% parameters for Gibbs sampler
        threshold_corrBB = 0;    % threshold for correlation matrix  0.01
        lower_bound_rand = 0.001;   % lower bound of random number
        upper_bound_rand = 0.999;   % upper bound of random number
        c_val_list=-25:3:25;  % range of c
        lambdaLearningRate_boosting = 0.1/sqrt(sigsq);  % 0.1
        %%%%
        
        %%%%%%%%%%%%%% multiple chains
        useMultiChain = true;
        nTileRow = 10; %nTileRow \times nTileCol defines the number of paralle chains
        nTileCol = 10;
        if useMultiChain == false
            nTileRow = 1;
            nTileCol = 1;
        end
        %%%%%%%%%%%%%
        
        %%% local normalization parameters
        isLocalNormalize = para.isLocalNormalize; %false; %
        isSeparate = false;
        localNormScaleFactor = 2; %0.5, 1, 2, 3, or 4, we suggest 2;
        thresholdFactor = 0.01;
        %%%%%%%%%%%%%%%%%%%%%%%%%%
        
        inPath = [para.dataPath categoryName];
        cachePath = ['./output/' para.name '/feature'];
        templatePath = ['./output/' para.name '/template'];
        resultPath = ['./output/' para.name '/result_task_' num2str(task_id) '_seed_' num2str(seed)];
        
        if ~exist(cachePath,'dir')
            mkdir(cachePath)
        else
            rmdir(cachePath,'s')
            mkdir(cachePath)
        end
        if ~exist(templatePath,'dir'),mkdir(templatePath),end
        if ~exist(resultPath,'dir')
            mkdir(resultPath)
            mkdir(fullfile(resultPath,'img'));
        else
            
            %     rmdir(resultPath,'s')
            mkdir(resultPath)
            mkdir(fullfile(resultPath,'img'));
        end
        
        
        %% Step 0: prepare filter, training images,and filter response on images
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
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% end of auxiliary variables
        
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
        
        %%%%%%%% optional: for drawing sketch in the learning process
        [allFilterR, ~, filterSymbol] = MakeFilterBank(GaborScaleList, DoGScaleList, nOrient); % filter bank for function "drawSketch"
        half = zeros(1, nOrient*nScaleGabor+nScaleDoG);
        for i=1:(nOrient*nScaleGabor+nScaleDoG)
            half(i) = (size(allFilterR{1, i}, 1)-1)/2;
        end
        %%%%%%%%%%
        
        files = dir(fullfile(inPath,'*.jpg'));
        numImage=length(files);
        disp(['start filtering']);
        
        for iImg = 1:numImage
            copyfile(fullfile(inPath,files(iImg).name),fullfile(resultPath,'img',files(iImg).name));
            imageOriginal = imread(fullfile(inPath,files(iImg).name));
            img = imresize(imageOriginal,[sx, sy]);
            if ndims(img) == 3
                img = rgb2gray(img);
            end
            img = im2single(img);
            
            allSizex = zeros(1, numResolution);
            allSizey = zeros(1, numResolution);
            ImageMultiResolution = cell(1, numResolution);
            
            for resolution=1:numResolution
                resizeFactor = 1.0 + (resolution - originalResolution)*scaleStepSize;
                img2 = imresize(img, resizeFactor, 'nearest');  % images at multiple resolutions
                img2 = img2-mean(img2(:));
                img2 = img2/std(img2(:))*sqrt(sigsq);
                
                ImageMultiResolution{resolution} = img2;
                
                [sizex, sizey] = size(ImageMultiResolution{resolution});
                allSizex(resolution) = sizex;
                allSizey(resolution) = sizey;
            end
            
            % compute MAX1 map for images in different resolutions
            disp(['======> start filtering and maxing image ' num2str(iImg)]);
            tic
            [SUM1mapFind, MAX1mapFind] = applyfilterBank_MultiResolution_sparseV4(ImageMultiResolution, filters, halfFilterSizes, nOrient,...
                locationShiftLimit,orientShiftLimit,isLocalNormalize,isSeparate,localNormScaleFactor,thresholdFactor,nScaleGabor,nScaleDoG, sqrt(sigsq));  % if using local normalization, the last parameter is important.
            
            mapName = fullfile(cachePath,['SUMMAXmap-image' num2str(iImg) '.mat']);
            save(mapName, 'imageOriginal', 'ImageMultiResolution','SUM1mapFind', 'MAX1mapFind','allSizex', 'allSizey');
            
            disp(['filtering time: ' num2str(toc) ' seconds']);
            
            
            % %     mapName = fullfile(cachePath,['SUMMAXmap-image' num2str(iImg) '.mat']);
            % %     current_file_name=files(iImg).name;
            % %     save(mapName, 'M1','current_file_name');  % only save the MAX1
        end
        
        
        %% Prepare variables for EM
        
        %%%%%%% initialization of alignment
        clusters=struct('imageIndex',cell(numCluster,1),'cropImage', cell(numCluster,1),'rHat',[],'template',[],'logZ',[], 'S2T', [], 'S3T', []);  % structure to store information of cluster

        MAX3scoreAll = rand(numImage, numCluster);   % randomly assign members to different cluster
        
        for c = 1:numCluster
            
            clusters(c).imageIndex=[];
            
            clusters(c).cropImage={};
            
            %%% initialize the observed statistics by setting zeros
            for iFilter = 1:numFilter
                clusters(c).rHat{iFilter}=zeros(sx, sy,'single');
            end
            
            t = 0; % index of image in the cluster, as well as the number of images in cluster
            for iImg = 1:numImage
                tic
                [~, ind] = max(MAX3scoreAll(iImg, :));
                if ind ~= c
                    continue;  % skip image that does not belong to cluster c
                end
                %clusteredImageIdx{c}=[clusteredImageIdx{c},iImg]; % collect the id for each cluster
                clusters(c).imageIndex=[clusters(c).imageIndex, iImg];
                
                t = t + 1;  % number of training images
                
                imageLoaded = load(fullfile(cachePath,['SUMMAXmap-image' num2str(iImg)]));
                
                % we initialize the alignment by cropping patch form original resolution image, with 0 level rotation and center (img_x/2. img_y/2)
                rot_init = 0;    % 0 level
                ind_init = originalResolution;  % original resolution
                Fx_init = floor(allSizex(originalResolution)/2);  % center x
                Fy_init = floor(allSizey(originalResolution)/2);   % center y
                
                cropedImage = single(zeros(sx, sy));
                Ccopy(cropedImage, single(imageLoaded.ImageMultiResolution{ind_init}), Fx_init, Fy_init, floor(sx/2), floor(sy/2), sx, sy, allSizex(ind_init), allSizey(ind_init), -rot_init*pi/nOrient);
                
                % optinal: output cropped images for iteration 0
                savingFolder0=fullfile(resultPath,['iteration0/morphedCropped/']);
                if ~exist(savingFolder0)
                    mkdir(savingFolder0);
                end
                
                gLow = min(cropedImage(:));
                gHigh = max(cropedImage(:));
                img_tem = (cropedImage-gLow)/(gHigh-gLow);
                imwrite(img_tem,fullfile(savingFolder0,['morphed-cluster-' num2str(c) '-img-' num2str(iImg,'%04d') '.png']));
                
                cropedImage = cropedImage - mean(cropedImage(:));
                cropedImage = cropedImage/std(cropedImage(:))*sqrt(sigsq);
                
                clusters(c).cropImage=[clusters(c).cropImage, double(cropedImage)];
                
                % compute feature map to learn
                [~, MAX1] = applyfilterBank_MultiResolution_sparseV4({cropedImage}, filters, halfFilterSizes, nOrient, locationShiftLimit,orientShiftLimit,...
                    isLocalNormalize,isSeparate,localNormScaleFactor,thresholdFactor,nScaleGabor,nScaleDoG, sqrt(sigsq));  % if using local normalization, the last parameter is important.
                
                % sum over the observed statistics (within cluster)
                for iFilter = 1:numFilter
                    clusters(c).rHat{iFilter}=clusters(c).rHat{iFilter}+MAX1{iFilter};
                end
                
                disp(['cropping time for image ' num2str(t) ' in cluster ' num2str(c) ': ' num2str(toc) ' seconds']);
            end% iImage
            
            disp(['Cluster ' num2str(c) ' has ' num2str(t) ' members.']);
            
            % average the observed statistics
            for iFilter = 1:numFilter
                clusters(c).rHat{iFilter}=clusters(c).rHat{iFilter}/t;
            end
            
        end
        it = 0;
        ShowClusteringAssignment;
        
        %% EM iteration
        for it = 1 : numEMIteration
            
            %%%% M-step
            disp(['M-step of iteration ' num2str(it)]);
            
            %% learn model for each cluster
            disp('Learning: learning sparse FRAME model for each cluster:');
            
            savingFolder=fullfile(resultPath,['iteration' num2str(it) '/']);
            if ~exist(savingFolder)
                mkdir(savingFolder);
            end
            
            % cpu paralleled learning
            if ~runParallel
                noWorkers = 1;
            end
            
           isOpen = isempty(gcp('nocreate'));
           if runParallel
              delete(gcp('nocreate'));
              parpool(noWorkers);
           end
           if ~runParallel && isOpen
              delete(gcp('nocreate'));
           end
            
            % assign tasks to the paralleled workers
            numAssignmentList1=ones(1,noWorkers)*floor(numCluster/noWorkers);
            Assignments_left=mod(numCluster,noWorkers);
            numAssignmentList2=[ones(1,Assignments_left), zeros(1,noWorkers-Assignments_left)];
            numAssignmentList=numAssignmentList1 + numAssignmentList2;
            endIndex=cumsum(numAssignmentList);
            startIndex=[1,endIndex(1,1:end-1)+1];
            
            ticID=tic;
            
            clear template currSample logZ deformedTemplate;
            
            temp_result=cell(numCluster,3);  % colum 1 is template, colum 2 is currSample, colum 3 is logZ
            spmd
                temp_result=codistributed(temp_result, codistributor1d(1,numAssignmentList));
                
                for c = startIndex(labindex):endIndex(labindex)
                    if ~isempty(clusters(c).imageIndex)
                        clusterSavingFolder=fullfile(savingFolder,['cluster' num2str(c) '/']);
                        if ~exist(clusterSavingFolder, 'dir')
                            mkdir(clusterSavingFolder);
                        end
                        
                        switch para.method
                            case 'one_stage'
                                
                                [template, currSample, logZ]=sparseFRAMElearnGibbs_multipleSelection(filters, nScaleGabor, nScaleDoG, nOrient, filterSymbol, half, clusters(c).rHat, ...
                                    sx, sy, halfFilterSizes, locationShiftLimit, nTileRow, nTileCol, lambdaLearningRate_boosting, numFilter, numWavelet, interval, Corr, threshold_corrBB, c_val_list, ...
                                    lower_bound_rand, upper_bound_rand, nPartCol, nPartRow, part_sx, part_sy, gradient_threshold_scale, clusterSavingFolder);
                                
                            case 'two_stage'
                                
                                disp('start filter selection');
                                [template, deformedTemplate]= filtersSelection_em(clusters(c).cropImage, GaborScaleList, DoGScaleList, nOrient, numSketch, numTopFeatureToShow,...
                                    locationShiftLimit, orientShiftLimit, sx, sy, clusterSavingFolder);
                                template.selectedLambdas=single(zeros(1,template.numSelected));
                                
                                disp('start learning Frame model');
                                [template,currSample, logZ]=sparseFRAMElearn(template, nIter, filters, clusters(c).rHat, sx, sy, halfFilterSizes, ...
                                    locationShiftLimit, nTileRow,nTileCol,epsilon,L,lambdaLearningRate_MP, numSample, nPartCol, nPartRow, part_sx, part_sy, isSaved, clusterSavingFolder);
                                
                            otherwise
                                error('Unkown method.');
                        end
                        
                        %%%%%%%%%%%%%%%%%
                        temp_result(c,1)={template};
                        temp_result(c,2)={currSample};
                        temp_result(c,3)={single(logZ)};
                    end
                end
                
            end
            
            temp_result=gather(temp_result);
            
            if isempty(gcp('nocreate')),delete(gcp('nocreate'));end
            
            % collect results and split template into several moving parts.
            for c = 1:numCluster
                if ~isempty(clusters(c).imageIndex)
                    clusters(c).template=temp_result{c,1};
                    clusters(c).currSample=temp_result{c,2};
                    clusters(c).logZ=temp_result{c,3};
                    
                    tic
                    [clusters(c).S2T, clusters(c).S3T] = hierachicalTemplate(numPart, part_sx, part_sy, sx, sy, rotateShiftLimit, nOrient, numRotate, clusters(c).template, nScaleGabor, partRotationRange,PartLocX, PartLocY);
                    disp(['spliting the learned template for cluster ' num2str(c) ': ' num2str(toc) ' seconds']);
                end
            end
            
            disp(['Total learning time of sparse FRAME for ' num2str(numCluster) ' clusters takes ' num2str(toc(ticID)) ' seconds']);
            
            
            %%%% E-step
            disp(['E-step of iteration ' num2str(it)]);
            
            detectedCroppedSavingFolder=fullfile(savingFolder, 'detectedCropped/');
            if ~exist(detectedCroppedSavingFolder)
                mkdir(detectedCroppedSavingFolder);
            end
            
            morphedCroppedSavingFolder=fullfile(savingFolder, 'morphedCropped/');
            if ~exist(morphedCroppedSavingFolder)
                mkdir(morphedCroppedSavingFolder);
            end
            
            boundingBoxSavingFolder=fullfile(savingFolder, 'boundingBox/');
            if ~exist(boundingBoxSavingFolder)
                mkdir(boundingBoxSavingFolder);
            end
            
            %%
            %%%%%%%%%%%%%%%% cpu paralleled learning
            noWorkers = para.noWorkers;
            if ~runParallel
                noWorkers = 1;
            end
            % Set matlabpool as needed
            
            isOpen = isempty(gcp('nocreate'));
%        isOpenCorr = matlabpool('size') == noWorkers;
            if runParallel
              delete(gcp('nocreate'));
              parpool(noWorkers);
            end
            if ~runParallel && isOpen
              delete(gcp('nocreate'));
            end
            
            % assign tasks to the paralleled workers
            numAssignmentList1=ones(1,noWorkers)*floor(numImage/noWorkers);
            Assignments_left=mod(numImage,noWorkers);
            numAssignmentList2=[ones(1,Assignments_left), zeros(1,noWorkers-Assignments_left)];
            numAssignmentList=numAssignmentList1 + numAssignmentList2;
            endIndex=cumsum(numAssignmentList);
            startIndex=[1,endIndex(1,1:end-1)+1];
            
            MAX3scoreAll = zeros(numImage, numCluster);
            spmd
                
                MAX3scoreAll = codistributed(MAX3scoreAll, codistributor1d(1, numAssignmentList));
                for iImg= startIndex(labindex):endIndex(labindex)
                    %for iImg = 1 : numImage
                    mapName = fullfile(cachePath,['SUMMAXmap-image' num2str(iImg)]);
                    imageLoaded = load(mapName);
                    
                    for c= 1:numCluster
                        
                        MAX3scoreAll(iImg, c) = detectObject(imageLoaded, clusters, c, iImg, numPart, numFilter, sx, sy, nOrient, part_sx, part_sy,...  % general parameters
                            relativePartLocationRange, argmaxMethod, allTemplateAffinityMatrix, resolutionShiftLimit, ...   % parameters for local max
                            rotationRange, partRotationRange, PartLocX, PartLocY, detectedCroppedSavingFolder, morphedCroppedSavingFolder, boundingBoxSavingFolder, RatioDisplacementSUM3);
                        
                    end % c
                end
                
            end
            MAX3scoreAll = gather(MAX3scoreAll);
            %%%%%%%%%%%%%%%%%
            
            %%
            disp('Cropping images and preparing feature maps for next iteration of learning:');
            
            copy_MAX3scoreAll=MAX3scoreAll;
            %%% preparation: collect training images for each cluster
            for c = 1:numCluster
                
                %clusteredImageIdx{c}=[];
                clusters(c).imageIndex=[];
                clusters(c).cropImage={};
                
                %%% initialize the observed statistics by setting zeros
                for iFilter = 1:numFilter
                    clusters(c).rHat{iFilter}=zeros(sx, sy,'single');
                end
                
                
                t = 0; % index of image in the cluster, as well as the number of images in cluster
                for iImg = 1:numImage
                    tic
                    [~, ind]=max(copy_MAX3scoreAll(iImg, :));
                    if ind~=c
                        continue;  % skip image that does not belong to cluster c
                    end
                    %clusteredImageIdx{c}=[clusteredImageIdx{c},iImg]; % collect the id for each cluster
                    clusters(c).imageIndex=[clusters(c).imageIndex, iImg];
                    
                    t = t + 1;
                    
                    %             if toUseCroppedImg
                    % load morphed cropped images for training
                    imageLoaded = load(fullfile(morphedCroppedSavingFolder,['morphed-cluster-' num2str(c) '-img-' num2str(iImg,'%04d')]), 'cropedMorphedImage');
                    cropedImage = imageLoaded.cropedMorphedImage;
                    
                    cropedImage = cropedImage - mean(cropedImage(:));
                    cropedImage = cropedImage/std(cropedImage(:))*sqrt(sigsq);
                    
                    % for re-learnig in matching pursuit
                    clusters(c).cropImage=[clusters(c).cropImage, double(cropedImage)];
                    
                    % compute feature map to learn
                    [~,MAX1] = applyfilterBank_MultiResolution_sparseV4({cropedImage}, filters, halfFilterSizes, nOrient, locationShiftLimit,orientShiftLimit,...
                        isLocalNormalize,isSeparate,localNormScaleFactor,thresholdFactor,nScaleGabor,nScaleDoG, sqrt(sigsq));  % if using local normalization, the last parameter is important.
                    %             else
                    %
                    %               % load morphed cropped images for training
                    %               imageLoaded = load(fullfile(morphedCroppedSavingFolder,['morphed-cluster-' num2str(c) '-img-' num2str(iImg,'%04d')]), 'croppedMorphedSUM1');
                    %               MAX1 = localmax(imageLoaded.croppedMorphedSUM1, nOrient, nScaleGabor, nScaleDoG, locationShiftLimit, orientShiftLimit);
                    %
                    %             end
                    
                    % sum over the observed statistics (within cluster)
                    for iFilter = 1:numFilter
                        clusters(c).rHat{iFilter}=clusters(c).rHat{iFilter}+MAX1{iFilter};
                    end
                    
                    disp(['cropping time for image ' num2str(t) ' in cluster ' num2str(c) ': ' num2str(toc) ' seconds']);
                end% iImage
                
                disp(['Cluster ' num2str(c) ' has ' num2str(t) ' members.']);
                
                % average the observed statistics
                for iFilter = 1:numFilter
                    clusters(c).rHat{iFilter}=clusters(c).rHat{iFilter}/t;
                end
                
            end
            
            ShowClusteringAssignment;  % show the cluster member before learning

            template_name = sprintf([templatePath '/template_task%d_seed%d_iter%d.mat'], task_id, seed, it);
            save(template_name, 'clusters', 'MAX3scoreAll');
            
        end %it

        tic_toc(task_id, seed) = toc(big_ticID);
    end
end

disp('done.');
tic_toc

