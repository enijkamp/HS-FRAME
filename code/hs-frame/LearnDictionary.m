function codeBook = LearnDictionary(imgCell, rootSavingFolder, para)


% activation: [imageID x y resolution templateInd score]


%% preparation

%ParameterCodeImage;

disp(['===============> Getting parameters, mex-c code, and training images']);

% if (exist('./output','dir'))%experimental results
%     delete('./output/*.*'); 
% else
%     mkdir('output');
% end


working_add=para.working_add;
if exist(working_add)
   rmdir(working_add,'s');
   mkdir(working_add); % working directory to store .mat files. 
else
   mkdir(working_add); % working directory to store .mat files. 
end

output_add=para.output_add;
output_add=[rootSavingFolder output_add];
if exist(output_add)
   rmdir(output_add,'s');
   mkdir(output_add); %  directory to store code book results. 
else
   mkdir(output_add); %  directory to store code book results. 
end

output_sparseFrame=para.output_sparseFrame;
output_sparseFrame=[rootSavingFolder output_sparseFrame];
if exist(output_sparseFrame)
   rmdir(output_sparseFrame,'s');
   mkdir(output_sparseFrame); % working directory to store synthesized images files. 
else
   mkdir(output_sparseFrame); % working directory to store synthesized images files. 
end


%% parameters for codebook
templateSize = para.templateSize;
% resizeFactor = [160 160]; 
minHeightOrWidth = para.minHeightOrWidth;
numRandomStart = para.numRandomStart;
numIter = para.numIter; % number of EM iterations
numCluster = para.numCluster; % number of data clusters
numSketch=para.numSketch; % number of Gabors in active basis at the first scale
numElementPointer = zeros(1, 2); % actual number selected, we use a pointer just to pass the number 
DELTA = 3; 
S2Thres = numSketch*DELTA; % cut-off value for detected instances of the template
% S2Thres =  -3.960e5;
% S2Thres2 =  -3.950e5;
rotateShiftLimit = para.rotateShiftLimit;
rotationRange = -rotateShiftLimit:rotateShiftLimit; % allowed global rotation of the template
locationPerturbationFraction = para.locationPerturbationFraction; % the size of neighborhood for MAX2 pooling, as well as surround supression
locationPerturbationFraction_final = para.locationPerturbationFraction_final; % (used in later EM iterations) the size of neighborhood for MAX2 pooling, as well as surround supression
SUM2mapBoundaryFraction = para.SUM2mapBoundaryFraction;

%% parameters for the outer square bounding box
templateSize = single(templateSize);
partSize = floor(sqrt(templateSize(1)*templateSize(2))); % alias for the template size (radius)
resizeTrainingImages = para.resizeTrainingImages; 
% constantImageArea = 200^2; % for resizing images (keep a constant area)
%% parameters for active basis
maxNumClusterMember = 50; % (no need to change) maximum number of training examples in each cluster used in re-learning the template
subsampleS2 = 1;  % subsampling step size for computing SUM2 maps
S1softthres = .0; % soft thresholding cutting-off for S1 maps
epsilon = .1; % allowed correlation between selected Gabors 
subsample = 1; subsampleM1 = 1; % subsample in computing MAX1 maps


numTopFeatureToShow=50;  % not effective
locationShiftLimit = para.locationShiftLimit; % shift in location 
orientShiftLimit = para.orientShiftLimit; % shift in orientation 
%% parameters for detection (controls scaling of templates)
resolutionGap = para.resolutionGap; % gap between consecutive resolutions in detection
numExtend = para.numExtend ; % number of gaps extended both below and above zero
numResolution = para.numResolution;
originalResolution = para.originalResolution; % original resolution is the one at which the imresize factor = 1
allResolution = para.allResolution;
%% parameters for filters (Gabor and DoG) 

DoGScaleList=para.DoGScaleList;
GaborScaleList=para.GaborScaleList;   % add scales of Gabor filters here


nScaleGabor=para.nScaleGabor;
nScaleDoG=para.nScaleDoG;
numOrient = para.numOrient;  % number of orientations
numFilter = para.numFilter;  % sin and cosin separately


%% parameters for HMC
sigsq=para.sigsq;
lambdaLearningRate = para.lambdaLearningRate;
nIter = para.nIter; % the number of iterations for learning lambda
epsilon = para.epsilon; % step size for the leapfrog
L = para.L; % leaps for the leapfrog
numSample = para.numSample; % how many HMC calls for each learning iteration
nTileRow = para.nTileRow;%12; %nTileRow \times nTileCol defines the number of paralle chains
nTileCol = para.nTileCol;%12;

%% parameters for normalization
isGlobalNormalization=para.isGlobalNormalization; % for the input training image. the images substracts the mean and is divided by standard deviation

isLocalNormalize=para.isLocalNormalize;
isSeparate=para.isSeparate;
localHalfx=para.localHalfx;
localHalfy=para.localHalfy;
localNormScaleFactor=para.localNormScaleFactor;
thresholdFactor=para.thresholdFactor;  


% % localOrNot = 1; % if we use local normalization or not. If not, set it to -1 
% % localHalfx = 20; localHalfy = 20; % the half range for local normalization, has to be quite large
% % windowNormalizeOrNot = -1; % whether normalize within the scanning window in detection 
% % if (localOrNot>0)
% %     windowNormalizeOrNot = -1;
% % end % if we use local normalization, we should not use window normalization in detection
% % thresholdFactor = .01;  % divide the response by max(average, maxAverage*thresholdFactor)
%% read in positive images
sizeTemplatex = templateSize(1);
sizeTemplatey = templateSize(2);
halfTemplatex = floor(sizeTemplatex/2);
halfTemplatey = floor(sizeTemplatey/2);
% imageFolder = 'positiveImage'; % folder of training images  
% imageName = dir([imageFolder '/*.jpg']);
numImage = numel(imgCell); % number of training images 
numImageTrain = numImage;
Ioriginal = cell(1, numImage);
numCI=zeros(numImage,numCluster);       % number of templates per image
numPOS=0;                               % number of positive images
totalACT=0;                             % total number of activated templates
totalL=0;                               % total log-likelihood
totalR=zeros(numImage,1);               % total area covered by the selected templates
totalSumR=zeros(numImage,1); 


for img = 1 : numImage
    tmpIm = imgCell{img};
%      tmpIm = imread([imageFolder '/' imageName(img).name]); 
    if size(tmpIm,3) == 3
        tmpIm = im2single(rgb2gray(tmpIm));
    else
        tmpIm = im2single(tmpIm);
    end

    sx = size(tmpIm,1); 
    sy = size(tmpIm,2);
    if resizeTrainingImages
       resizeFactor = minHeightOrWidth/min(sx,sy);
       tmpIm = imresize( tmpIm, resizeFactor, 'bilinear' );
    end
    
%     if isGlobalNormalization
%        tmpIm=tmpIm-mean(tmpIm(:));
%        tmpIm = tmpIm/std(tmpIm(:))*sqrt(sigsq);
%     end
    
    Ioriginal{img} = tmpIm;
    J0 = Ioriginal{img};
    J = cell(1, numResolution);
    
    J_raw = cell(1, numResolution); %% new !!
    
    allSizex = zeros(1, numResolution); 
    allSizey = zeros(1, numResolution);
    for r=1:numResolution
       img_scaled = imresize(J0, allResolution(r), 'nearest');  % scaled images
       
       J_raw{r} = img_scaled; %% new !!
       
       if isGlobalNormalization
           img_scaled = img_scaled-mean(img_scaled(:));
           img_scaled = img_scaled/std(img_scaled(:))*sqrt(sigsq);
       end
       
       J{r} = img_scaled;
       [sizex, sizey] = size(J{r}); 
       allSizex(r) = sizex; 
       allSizey(r) = sizey; 
    end
    multipleResolutionImageName = [working_add '/multipleResolutionImage' num2str(img)]; 
    save(multipleResolutionImageName, 'J','J_raw','allSizex','allSizey');
end

count = 0;
% generate the set of geometric transformations for each template
for templateScaleInd = 0:0 % large scale change
    for rotation = rotationRange 
        for rowScale = 2.^[0] % small scale change and aspect ratio change
            for colScale = 2.^[0]
                count = count + 1;
                templateTransform{count} = [templateScaleInd rowScale colScale rotation];
            end
        end
    end
end
nTransform = count;

% prepare morph-back mappings for geometric transforms
largestPartSizeX = templateSize(1); largestPartSizeY = templateSize(2);
denseX = -floor(largestPartSizeX/2) + (1:largestPartSizeX);
denseY = -floor(largestPartSizeY/2) + (1:largestPartSizeY);
count = 0;
inRow = zeros(length(denseX)*length(denseY),1,'single');
inCol = zeros(length(denseX)*length(denseY),1,'single');
inO = zeros(numel(inRow),1,'single'); 
inS = zeros(numel(inRow),1,'single');
for y = denseY
    for x = denseX
        count = count+1;
        inRow(count) = x;
        inCol(count) = y;
    end
end
outRow = cell(nTransform,1);
outCol = cell(nTransform,1);
for iT = 1:nTransform
	tScale = templateTransform{iT}(1); 
	rScale = templateTransform{iT}(2);
	cScale = templateTransform{iT}(3);
	[outRow{iT}, outCol{iT}] = ...
		mexc_TemplateAffineTransform(tScale,rScale,cScale,...
		    templateTransform{iT}(4),inRow,inCol,inO,inS,numOrient);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%MakeFiltersBank;

% create filters bank

filters=para.filters;
halfFilterSizes = para.halfFilterSizes;


storeFiltersBankName = [working_add '/storeFiltersBank.mat'];
save(storeFiltersBankName, 'filters', 'halfFilterSizes', 'numFilter', ...
                    'nScaleGabor', 'nScaleDoG');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

CreatingFeaturesMaps;

%% begin clustering with window scanning in E step

% initialize: generate random SUM2 maps and thus random cluster members
mixing = zeros(numCluster,1); % number of examples in each cluster
aveLogL = zeros(numCluster,1); % average log likelihood in each cluster
bestOverallScore = -inf;
bestModel=cell(numCluster,1);
%buf_length = 0;
for iRS = 1:numRandomStart
    
    %% initial step
    activations = []; % 3 by N matrix, where N is an unknown large number
    for i = 1:numImage
        % compute SUM2
        SUM1MAX1mapName = [working_add '/SUM1MAX1map' 'image' num2str(i)];
        load(SUM1MAX1mapName, 'SUM1map');
        SUM2map = cell(nTransform*numCluster,numResolution);
        for iRes = 1:numResolution
		    width = floor(size(SUM1map{iRes,1},2)/subsampleS2);
		    height = floor(size(SUM1map{iRes,1},1)/subsampleS2);
		    for j = 1:nTransform*numCluster
		        SUM2map{j,iRes} = single( rand( height, width ) );
		    end
        end
        
        %save('SUM2map_initialized.mat','SUM2map');
        %load('SUM2map_initialized.mat','SUM2map');
        
        % exclude the near-boundary region
        for ii = 1:numel(SUM2map)
            SUM2map{ii}(:,[1:floor(templateSize(2)/SUM2mapBoundaryFraction) end-floor(templateSize(2)/SUM2mapBoundaryFraction):end]) = -1001; %min(-1,S2Thres-1);
            SUM2map{ii}([1:floor(templateSize(1)/SUM2mapBoundaryFraction) end-floor(templateSize(1)/SUM2mapBoundaryFraction):end],:) = -1001; %min(-1,S2Thres-1);
        end

        % compute MAX2, perform surround supression and get activations                  
        tmpActivations = mexc_ComputeMAX2MP( SUM2map, int32(locationPerturbationFraction*partSize/subsampleS2), 0);
        disp(['crop: ' num2str(size(tmpActivations,2)) ' instances in mexc_ComputeMAX2MP']);
     
        activations = [activations,[single(i*ones(1,size(tmpActivations,2)));tmpActivations]];
    end
    activations(2:3,:) = activations(2:3,:) * subsampleS2;
    activatedCluster = ceil( ( activations(5,:) + 1 ) / nTransform );  % index of original activations(5,:) starts from 0
    activatedTransform = activations(5,:) + 1 - (activatedCluster-1) * nTransform;
    activatedImg = activations(1,:);
    initialClusters = activations;

    save( [working_add '/activations0.mat'], 'activations' );
    
    numElementAll = zeros(1, numCluster); 
    for iter = 1:numIter
        %% M step
        syms = cell(numCluster,1);
        syms_singleChain = cell(numCluster,1);
        sketchTemplate = cell(numCluster,1);
        rHatSet=cell(numCluster,1);
        croppeds=cell(numCluster,1);
        toLearn=zeros(numCluster,1);
        
        ticID=tic;
        for cc = 1:numCluster
            
            savingFolder=[output_sparseFrame '/iteration' num2str(iter) '/cluster' num2str(cc) '/'];
            if exist(savingFolder)
               rmdir(savingFolder,'s'); 
               mkdir(savingFolder);
            else
               mkdir(savingFolder);
            end
            
%             for k = 1:buf_length
% 				fprintf(1,'\b');
% 			end
            disp(['run ' num2str(iRS) ': learning iteration ' num2str(iter) ' for cluster ' num2str(cc)]);
% 			fprintf(1,str); drawnow;
% 			buf_length = length(str);
            
            % =====================================
            % crop back and relearn
            % =====================================
            ind = find(activatedCluster == cc);
            mixing(cc) = length(ind);
            aveLogL(cc) = mean(activations(6,ind));
            if isnan(aveLogL(cc))
            	aveLogL(cc) = -1;
            end
            % sample a subset of training postitives, if necessary
            if length(ind) > maxNumClusterMember
                idx = randperm(length(ind));
                ind = ind(idx(1:maxNumClusterMember));
                ind = sort(ind,'ascend');
            end
            nMember = length(ind);
            
            SUM1mapLearn = cell(nMember,numFilter);
            MAX1mapLearn = cell(nMember,numFilter);
            for iFilter = 1:numFilter
                for iMember = 1:nMember
                   MAX1mapLearn{iMember,iFilter} = zeros(templateSize(1),templateSize(2),'single');
                end
            end
            
            cropped = cell(nMember,1);
            currentImg = -1;
            
            for iMember = 1:nMember
                if activatedImg(ind(iMember)) ~= currentImg
                    currentImg = activatedImg(ind(iMember));
                    SUM1MAX1mapName = [working_add '/SUM1MAX1map' 'image' num2str(currentImg)];
                    load(SUM1MAX1mapName, 'SUM1map', 'J' );
                end
                % use mex-C code instead: crop S1 map
				% warning: if the templates are of different sizes, we need to alter outRow and outCol (use only a subset of it).                
                
                
                % for Gabor features (each scale includes sin part and cosin part)  pay attention! Please set tScale=0.   
                numGaborFilter=nScaleGabor*2*numOrient;
                tScale = 0; destHeight = templateSize(1); destWidth = templateSize(2); nScale = nScaleGabor*2; reflection = 1;
                SUM1mapLearn(iMember,1:numGaborFilter) = mexc_CropInstance_FRAME( SUM1map(1+activations(4,ind(iMember)),1:numGaborFilter),...
                    activations(2,ind(iMember))-1,...
                    activations(3,ind(iMember))-1,...
                    rotationRange(activatedTransform(ind(iMember))),tScale,reflection,...
                    outRow{activatedTransform(ind(iMember))},outCol{activatedTransform(ind(iMember))},...
                    numOrient,nScale,destWidth,destHeight );
                
                if nScaleDoG>0
                % for DoG at the end of the feature map (one orientation)
                tScale = 0; destHeight = templateSize(1); destWidth = templateSize(2); nScale = nScaleDoG; reflection = 1; numOrientDoG = 1;
                SUM1mapLearn(iMember,numGaborFilter+1:end) = mexc_CropInstance_FRAME( SUM1map(1+activations(4,ind(iMember)),numGaborFilter+1:end),...
                    activations(2,ind(iMember))-1,...
                    activations(3,ind(iMember))-1,...
                    rotationRange(activatedTransform(ind(iMember))),tScale,reflection,...
                    outRow{activatedTransform(ind(iMember))},outCol{activatedTransform(ind(iMember))},...
                    numOrientDoG,nScale,destWidth,destHeight );
                end
                
        
                % Crop detected image patch for visualization
                srcIm = J{1+activations(4,ind(iMember))};
                tmpNumOrient = 1; nScale=1;
                cropped(iMember) = mexc_CropInstance_FRAME( {single(srcIm)},...
                    activations(2,ind(iMember))-1,...
                    activations(3,ind(iMember))-1,...
                    rotationRange(activatedTransform(ind(iMember))),tScale,reflection,...
                    outRow{activatedTransform(ind(iMember))},outCol{activatedTransform(ind(iMember))},...
                    tmpNumOrient,nScale,destWidth,destHeight );
                
                % local max                              
                for iScale=1:nScaleGabor                    
                    % sin part
                    CgetMAX1(1,double(templateSize(1)),double(templateSize(2)),numOrient,locationShiftLimit,orientShiftLimit,1,...
                        SUM1mapLearn(iMember,1+(iScale-1)*(2*numOrient):numOrient+(iScale-1)*(2*numOrient)),...
                        MAX1mapLearn(iMember, 1+(iScale-1)*(2*numOrient):numOrient+(iScale-1)*(2*numOrient)));
                    
                    % cos part
                    CgetMAX1(1,double(templateSize(1)),double(templateSize(2)),numOrient,locationShiftLimit,orientShiftLimit,1,...
                        SUM1mapLearn(iMember, numOrient+1+(iScale-1)*(2*numOrient): 2*numOrient+(iScale-1)*(2*numOrient)),...
                        MAX1mapLearn(iMember, numOrient+1+(iScale-1)*(2*numOrient): 2*numOrient+(iScale-1)*(2*numOrient)));         
                end                
                % DoG
                for iScale=1:nScaleDoG
                    ind  = 2*numOrient*nScaleGabor+iScale;       
                    CgetMAX1_DoG(1,double(templateSize(1)),double(templateSize(2)),locationShiftLimit,1,SUM1mapLearn(iMember, ind)',MAX1mapLearn(iMember,ind)');   
                end
   
            end
            
            
            for iFilter = 1:numFilter
                rHat{iFilter}=zeros(templateSize(1), templateSize(2),'single');
            end
    
            for iMember = 1:nMember
               for iFilter = 1:numFilter
                rHat{iFilter}= rHat{iFilter}+ MAX1mapLearn{iMember,iFilter}/nMember;
               end
            end
            
            rHatSet{cc}=rHat;
            
            %% changed into double for filter selection
            croppeds{cc}=cell(1,length(cropped));
            for i=1:length(cropped)
               croppeds{cc}{1,i}=double(cropped{i});
            end
                        
            
            im = displayImages(cropped,10,60,60);
            if ~isempty(im)
                imwrite(im,sprintf([output_add '/iter%d_cluter%d.png'],iter,cc));
                toLearn(cc)=true;
            else
                syms{cc} = zeros(templateSize,'single');
                toLearn(cc)=false;
                if iter>1
                	copyfile([working_add sprintf('/learnedmodel%d_iter%d.mat',cc,iter-1)],[working_add sprintf('/learnedmodel%d_iter%d.mat',cc,iter)]);
                else
                	logZ = 1;
		            save([working_add sprintf('/learnedmodel%d_iter%d.mat',cc,iter)], 'logZ');
                end
                                
            end
        end
        disp(['time for cropping: ' num2str(toc(ticID)) ' seconds']); 
               
        
        
       %%% cpu paralleled learning     
        % Set variable to know when to run parallel
       runParallel = 1; 
       % Set number of workers
       noWorkers = 3;   %% depend on how many cores in your cpu     
       %% Set matlabpool as needed
       isOpen = matlabpool('size') > 0;
       isOpenCorr = matlabpool('size') == noWorkers;
       if runParallel && ~isOpenCorr,
          matlabpool close force
          matlabpool(noWorkers)
       end
       if ~runParallel && isOpen,
          matlabpool close force
       end
 
       % assign tasks to the paralleled workers
       numAssignmentList1=ones(1,noWorkers)*floor(numCluster/noWorkers);
       Assignments_left=mod(numCluster,noWorkers);
       numAssignmentList2=[ones(1,Assignments_left), zeros(1,noWorkers-Assignments_left)];
       numAssignmentList=numAssignmentList1 + numAssignmentList2;
       endIndex=cumsum(numAssignmentList);
       startIndex=[1,endIndex(1,1:end-1)+1];
        
        
       ticID=tic;
       storeFiltersBankName = [working_add '/storeFiltersBank.mat'];
       load(storeFiltersBankName, 'filters', 'halfFilterSizes');
        
        temp_result=cell(numCluster,6);  % colum 1 is template, colum 2 is currSample, colum 3 is logZ, colum 4 is templateMP, colum 5 is SketchTemplate, colum 6 is currSample
        spmd          
          
          temp_result=codistributed(temp_result, codistributor1d(1,numAssignmentList));  
            
          for cc = startIndex(labindex):endIndex(labindex)                          
            % now start re-learning
            savingFolder=[output_sparseFrame '/iteration' num2str(iter) '/cluster' num2str(cc) '/'];
            if exist(savingFolder)
               rmdir(savingFolder,'s'); 
               mkdir(savingFolder);
            else
               mkdir(savingFolder);
            end
            
            % filter selection
            if toLearn(cc)
               ticID1=tic;                        
               [template, templateMP, currentSketchTemplate]= filtersSelection_codebook(croppeds{cc}, GaborScaleList, DoGScaleList, numOrient, numSketch,...
                   locationShiftLimit, orientShiftLimit, double(templateSize(1)),double(templateSize(2)), isLocalNormalize, localNormScaleFactor, thresholdFactor, savingFolder);
               template.selectedLambdas=single(zeros(1,template.numSelected));            
               disp(['finished filters selection for cluster ' num2str(cc) ': '  num2str(toc(ticID1)) ' seconds']);
               %sketchTemplate{cc}=currentSketchTemplate;
            
               % sparse FRAME learning
               ticID2=tic;
               isSaved=1;               
               [template,currSample,logZ]=sparseFRAMElearn(template, nIter, filters, rHatSet{cc}, double(templateSize(1)), double(templateSize(2)), halfFilterSizes, ...
                    locationShiftLimit, nTileRow,nTileCol,epsilon,L,lambdaLearningRate,numSample, isSaved,savingFolder);
               disp(['finished filters selection for cluster ' num2str(cc) ': '  num2str(toc(ticID2)) ' seconds']);
               %syms{cc}=currSample; 
               
               temp_result(cc,1)={template};
               temp_result(cc,2)={currSample};
               temp_result(cc,3)={logZ};  
               temp_result(cc,4)={templateMP};    
               temp_result(cc,5)={currentSketchTemplate};  
               temp_result(cc,6)={currSample};   
               
               %save([working_add sprintf( '/learnedmodel%d_iter%d.mat',cc,iter)], 'template', 'logZ', 'currSample','templateMP');
            end
          end 
        end
        temp_result=gather(temp_result);
        
        if matlabpool('size') > 0;
            matlabpool close
        end
        disp(['All Models Learning time: ' num2str(toc(ticID)) ' seconds']); 
             
        % saving template and results of sparse-Frame model
        for cc=1:numCluster
           if toLearn(cc) 
           % the following results are computed by parelleled computing, they cannot be saved properly within the process of parelleled computing  
           template=temp_result{cc,1};
           logZ=temp_result{cc,3};
           currSample=temp_result{cc,2};
           templateMP=temp_result{cc,4};
           save([working_add sprintf( '/learnedmodel%d_iter%d.mat',cc,iter)], 'template', 'logZ', 'currSample', 'templateMP'); 
           
           syms{cc} = temp_result{cc,6};
           sketchTemplate{cc} = temp_result{cc,5};    
           
           % multiple chains of synthesized images
           img=syms{cc};
           gLow = min(img(:));
           gHigh = max(img(:));
           imgSaved = (img-gLow)/(gHigh-gLow);  
           
           % single chain of synthesized image
           crop = img(1:templateSize(1), 1:templateSize(2));
           syms_singleChain{cc} = crop;
           gLow = min(crop(:));
           gHigh = max(crop(:));
           cropSaved = (crop-gLow)/(gHigh-gLow); 
           
           % sketch template
           img=sketchTemplate{cc};
           gLow = min(img(:));
           gHigh = max(img(:));
           sketchSaved = (img-gLow)/(gHigh-gLow); 
           %
           
           %% C = inv(B'B)B'I;  I_reconstruct = B*C = B* inv(B'B)B' I_syn, we compute and store R=B inv(B'B)B' as reconstruction matrix
           B=[];
           for i=1:template.numSelected    
                  Bi = filterTiling(filters{template.selectedFilter(i)}, templateSize(1), templateSize(2), template.selectedRow(i), template.selectedCol(i));
                  B = [B,Bi(:)];
           end
           R=B*(inv(B'*B)*(B'));
           img_projected = reshape(R*crop(:),templateSize(1), templateSize(2));  % R = B*(inv(B'*B)*(B'));
           gLow = min(img_projected(:));
           gHigh = max(img_projected(:));
           img_projected = (img_projected-gLow)/(gHigh-gLow);     
                      
           imwrite(imgSaved,sprintf([output_add '/iter%d_MultipleChainTemplate%d.png'],iter,cc)); 
           imwrite(cropSaved,sprintf([output_add '/iter%d_SingleTemplate%d.png'],iter,cc));
           imwrite(sketchSaved,sprintf([output_add '/iter%d_SketchTemplate%d.png'],iter,cc));
           imwrite(img_projected,sprintf([output_add '/iter%d_projectedTemplate%d.png'],iter,cc));
           end
        end
        
        % ==============================================
        %% E step
        % ==============================================
      flipOrNot=0;
      numTransform=(2*rotateShiftLimit+1)*(1+flipOrNot);
      AllTemplate=cell(numTransform,numCluster);
      
      
      %%%%% test
        for cc = 1:numCluster
              
             disp(['rotating template of cluster' num2str(cc)]);
             load([working_add sprintf( '/learnedmodel%d_iter%d.mat',cc,iter)], 'template', 'logZ');           
             
        
              % shift the origin to the center of the template
             selectedRow = template.selectedRow - round(templateSize(1)/2) +1;
             selectedCol = template.selectedCol - round(templateSize(2)/2) +1;
             selectedOri = mod(template.selectedFilter -1, numOrient ); % Ori index starts from 0;  
             selectedScale = floor((template.selectedFilter -1)/numOrient );  % Scale index starts from 0; 
             S2Templates{cc} = struct('selectedRow',selectedRow, 'selectedCol',selectedCol,'selectedScale',selectedScale,'selectedOri', selectedOri);
          
        end
          
        for cc = 1:numCluster 
            numElement=template.numSelected;
            load([working_add sprintf( '/learnedmodel%d_iter%d.mat',cc,iter)], 'template', 'logZ');   
          
            for iT = 1:nTransform
                templateScaleInd = templateTransform{iT}(1);
                rowScale = templateTransform{iT}(2);
                colScale = templateTransform{iT}(3);
                rotation = templateTransform{iT}(4);
                [tmpSelectedRow tmpSelectedCol tmpSelectedOri tmpSelectedScale] = ...
                    mexc_TemplateAffineTransform( templateScaleInd, rowScale,...
                    colScale, rotation, S2Templates{cc}.selectedRow, S2Templates{cc}.selectedCol,...
                    S2Templates{cc}.selectedOri, S2Templates{cc}.selectedScale, numOrient );
                AllTemplate{iT,cc}.selectedRow = single(tmpSelectedRow)';
                AllTemplate{iT,cc}.selectedCol = single(tmpSelectedCol)';
                AllTemplate{iT,cc}.selectedFilter = single(tmpSelectedOri+tmpSelectedScale*numOrient+1)';               
                AllTemplate{iT,cc}.selectedLambdas = template.selectedLambdas;
                AllTemplate{iT,cc}.logZ = single(logZ);
            end
        end
        
      
        
        if iter == numIter
            codeBook = AllTemplate;
            %save(sprintf('./output/learnedCodeWord.mat'),codeBook);
        end
        
      
      activations = []; % 3 by N matrix, where N is an unknown large number
        
      for i=1:numImage
          disp(['starting to compute SUM2map for Images' num2str(i)]);
          
          ticID=tic;
          
        % compute SUM2
          SUM1MAX1mapName = [working_add '/SUM1MAX1map' 'image' num2str(i)];
          load(SUM1MAX1mapName, 'MAX1map', 'allSizex','allSizey');
          
          SUM2map = cell(nTransform*numCluster,numResolution);
          
          for iCluster=1:numCluster           
            for r=1:numTransform
               [SUM2] = sparseFRAME_SUM2_logZ( single(numResolution), single(allSizex), single(allSizey), single(numFilter), AllTemplate{r,iCluster}, MAX1map);  
               SUM2map(r+(iCluster-1)*numTransform,:)=SUM2;
            end
          end
          
          % random perturbation (to break ties arbitrarily for MAX2)
%           for ii = 1:numel(SUM2map)
%                 SUM2map{ii}(:) = SUM2map{ii}(:) + 1e-3 * ( rand(numel(SUM2map{ii}),1) - .5 );
%           end
          
          % exclude the near-boundary region
          for ii = 1:numel(SUM2map)
            	SUM2map{ii}(:,[1:floor(templateSize(2)/SUM2mapBoundaryFraction) end-floor(templateSize(2)/SUM2mapBoundaryFraction):end]) = -1001; %min(-1,S2Thres-1);
            	SUM2map{ii}([1:floor(templateSize(1)/SUM2mapBoundaryFraction) end-floor(templateSize(1)/SUM2mapBoundaryFraction):end],:) = -1001; %min(-1,S2Thres-1);
          end
                    
          disp(['SUM2map computational time: ' num2str(toc(ticID)) ' seconds']);           
         % save( ['SUM2map' num2str(iter) '.mat'], 'SUM2map' );
                    
          
          disp(['starting to compute MAX2map for Images' num2str(i)]);
          tic
          % compute MAX2
           if iter ==1              
                tmpActivations = mexc_ComputeMAX2MP( SUM2map, int32(locationPerturbationFraction*partSize/subsampleS2),  0);
                disp(['crop: ' num2str(size(tmpActivations,2)) ' instances in mexc_ComputeMAX2MP']);
              
           else
                % discard the activated instances that have a low S2 score
                if iter > floor(numIter/2) % for the later iterations, increase the sparsity
                    locationPerturbationFraction = locationPerturbationFraction_final;
                end
               
                tmpActivations = mexc_ComputeMAX2MP( SUM2map, int32(locationPerturbationFraction*partSize/subsampleS2), 0);
                disp(['crop: ' num2str(size(tmpActivations,2)) ' instances in mexc_ComputeMAX2MP']);               
                
                %tmpActivations=tmpActivations(:,1:5);  
           end
           
           disp(['MAX2 computational time: ' num2str(toc) ' seconds']);
           activations = [activations,[single(i*ones(1,size(tmpActivations,2)));tmpActivations]];   
         
       end
        
        save( [working_add '/activations' num2str(iter) '.mat'], 'activations' );
        
        activations(2:3,:) = activations(2:3,:) * subsampleS2;
        activatedCluster = ceil( ( activations(5,:) + 1 ) / nTransform );
        activatedTransform = activations(5,:) + 1 - (activatedCluster-1) * nTransform;
        activatedImg = activations(1,:);
        % disp(sprintf('on average %.2f activations per image',size(activations,2)/numImage));
  
    end

    %fprintf(1,'\n');

    % compute overall Score
	scores = activations(6,:);
    scores = max(scores,0);
	overallScore = sum(scores); % use mean or sum ?
	
	if overallScore > bestOverallScore  % if we do multiple random startings, this will save the result from the best starting point 
        
        for cc = 1:numCluster
            load([working_add sprintf( '/learnedmodel%d_iter%d.mat',cc,numIter)], 'template', 'logZ', 'templateMP');
            bestModel{cc}=struct('templateFRAME',template, 'logZ',logZ, 'templateMP',templateMP);
        end
        
		bestOverallScore = overallScore;
		bestSyms = syms;
        bestSymsSingle = syms_singleChain;
		bestInitialClusters = initialClusters;
		bestActivations = activations;
        bestMixing = mixing;
        bestAveLogL = aveLogL;
	end

end


fprintf(1,'\n');


% % % % 
% % % % % -- now we have selected the best random starting point --
% % % % %% display the templates and cluster members for the best random starting point
% % % % activations = bestActivations;
% % % % syms = bestSyms;
% % % % syms_singleChain = bestSymsSingle; 
% % % % activatedImg = activations(1,:);
% % % % activatedCluster = ceil( ( activations(5,:) + 1 ) / nTransform ); % starts from 1
% % % % activatedTransform = activations(5,:) + 1 - (activatedCluster-1) * nTransform; % starts from 1
% % % % mixing = bestMixing;
% % % % aveLogL = bestAveLogL;
% % % % cluster_is_nonempty = zeros(numCluster,1);
% % % % for cc = 1:numCluster
% % % % 	ind = find(activatedCluster == cc);
% % % % 	% sample a subset of training postitives, if necessary
% % % % 	if length(ind) > maxNumClusterMember
% % % % 		idx = randperm(length(ind));
% % % % 		ind = ind(idx(1:maxNumClusterMember));
% % % % 		ind = sort(ind,'ascend'); % make sure the imageNo is still in ascending order
% % % % 	end
% % % % 	
% % % % 	nMember = length(ind);
% % % % 	cropped = cell(nMember,1);
% % % % 	currentImg = -1;
% % % % 	for iMember = 1:length(ind)
% % % % 		if activatedImg(ind(iMember)) ~= currentImg
% % % % 			currentImg = activatedImg(ind(iMember));
% % % % 			SUM1MAX1mapName = [working_add '/SUM1MAX1map' 'image' num2str(currentImg)];
% % % %             load(SUM1MAX1mapName, 'J' );
% % % %         end
% % % % 
% % % %         tScale = 0; destHeight = templateSize(1); destWidth = templateSize(2); nScale = 1; reflection = 1; tmpNumOrient = 1;
% % % % 		% Crop detected image patch for visualization
% % % % 		srcIm = J{1+activations(4,ind(iMember))};
% % % %         cropped(iMember) = mexc_CropInstance_FRAME( {single(srcIm)},...
% % % %             activations(2,ind(iMember))-1,...
% % % %             activations(3,ind(iMember))-1,...
% % % %             rotationRange(activatedTransform(ind(iMember))),tScale,reflection,...
% % % %             outRow{activatedTransform(ind(iMember))},outCol{activatedTransform(ind(iMember))},...
% % % %             tmpNumOrient,nScale,destWidth,destHeight );	          
% % % % 	end
% % % % 	
% % % % 	if ~isempty(ind) % cluster is not empty
% % % % 		im = displayImages(cropped,10,templateSize(1),templateSize(2));
% % % % 		imwrite(im,sprintf([output_add '/final_cluster%d.png'],cc));
% % % % 		cluster_is_nonempty(cc) = 1;
% % % % 		%syms{cc} = -single(S2Templates{cc}.commonTemplate);
% % % % 	end
% % % % end
% % % % 
% % % % % calculate the number of templates per image
% % % % for ii=1:numImage
% % % %     ind=find(activations(1,:)==ii);
% % % %     activatedCluster = ceil( ( activations(5,ind) + 1 ) / nTransform );
% % % %     numCI(ii,:)=hist(activatedCluster,1:numCluster);
% % % % end;
% % % % numPOS=numImage;
% % % % totalACT=size(activations,2);
% % % % totalL=sum(bestAveLogL);
% % % % save([working_add '/learning_result.mat'],'bestActivations','bestSyms','bestSymsSingle','bestOverallScore','bestInitialClusters','bestAveLogL','bestMixing');
% % % % 
% % % % % rank the learned templates
% % % % nonempty_clusters = find(cluster_is_nonempty);
% % % % mixing = mixing(nonempty_clusters);
% % % % aveLogL = aveLogL(nonempty_clusters);
% % % % syms_singleChain = syms_singleChain(nonempty_clusters);
% % % % syms_singleChain2 = syms_singleChain;
% % % % numElementAll = numElementAll(nonempty_clusters);
% % % % [sorted idx] = sort( (mixing) .* aveLogL, 'descend' );  % CHANGE XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
% % % % numElementAllsort = numElementAll(idx); 
% % % % 
% % % % for i = 1:numel(sorted)
% % % %     towrite = syms_singleChain{idx(i)};    
% % % %     
% % % %     if range(towrite) < 1
% % % %         towrite(:) = 255;
% % % %     else
% % % %         towrite = uint8(255 * (towrite-min(towrite(:)))/(max(towrite(:))-min(towrite(:))));
% % % %         towrite = double(towrite) - 50;
% % % %     end
% % % %      
% % % %     syms_singleChain2{i} = towrite;
% % % % end
% % % % towrite = displayImages( syms_singleChain2, 10, templateSize(1), templateSize(2), false );
% % % % imwrite(towrite,sprintf([output_add '/template_sorted.png']));
% % % % tempOrder=idx;
% % % % numActiveClusters=length(nonempty_clusters);
% % % % 
% % % % % show the color templates
% % % % colors = colormap(hsv(numCluster));
% % % % non_empty_cluster=find(cluster_is_nonempty);
% % % % for cc=1:length(non_empty_cluster)
% % % %     towrite=zeros(size(syms_singleChain{cc},1),size(syms_singleChain{cc},2),3)+255;
% % % %     [posx,posy]=find(syms_singleChain{cc}(:,:)<-0.1);
% % % %     for jj=1:length(posx)
% % % %         cpos=find(cc==tempOrder);
% % % %         towrite(posx(jj),posy(jj),:)=colors(cpos,:);
% % % %     end;
% % % %     imwrite(towrite,sprintf([output_add '/color_template%d.png'],cc));
% % % % end;
% % % % 
% % % % tr_or_tt='train';
% % % % displayActivations;
% % % % 
% % % % % % % if there are test images, execute on the test images
% % % % % % if exist('./testingImage','dir')  
% % % % % %    TestingProcess;
% % % % % % end;
% % % % % % 
% % % % % % tr_or_tt='testing';
% % % % % % bestActivations=test_activations;
% % % % % % displayActivations;
% % % % % % 
% % % % % % % show html file
% % % % % % GenerateHtml;
end



