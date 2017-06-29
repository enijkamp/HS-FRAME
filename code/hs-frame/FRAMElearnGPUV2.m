function [lambdaF,currSamples, LogZ]=FRAMElearnGPUV2(nIter,filters,rHat,initialSample,initialLambda,initialLogZ,epsilon,L,lambdaLearningRate,numSample,isSaved,savingFolder,isComputelogZ)
% extract auxilary variables
numFilter = numel(filters);
[sx, sy]=size(rHat{1});
halfFilterSizes = zeros(numFilter,1);
for iFilter = 1:numFilter
    halfFilterSizes(iFilter)=(size(filters{iFilter},1)-1)/2;
end
overAllArea = sum((sx-2*halfFilterSizes).*(sy-2*halfFilterSizes));
nTileRow = size(initialSample,1)/sx;
nTileCol = size(initialSample,2)/sy;


% sample 
if isempty(initialSample)
    prevSamples = randn(sx*nTileRow,sy*nTileCol,'single');
else
    prevSamples = initialSample;
end


% lambda and logZ
gradientF = cell(numFilter,1);
lambdaF = cell(numFilter,1);
if isempty(initialLambda)
    for iFilter = 1:numFilter
        lambdaF{iFilter} = zeros(sx,sy,'single');
    end
    initialLogZ =log((2*pi))*(overAllArea/2);
else
   for iFilter = 1:numFilter
        lambdaF{iFilter} = initialLambda{iFilter};
    end
end
for iFilter = 1:numFilter
        % the following lines will be useful when lambda is not initialized as zero;
        %{
        h = halfFilterSizes(iFilter);
        lambdaF{iFilter}([1:h,sx-h+1:sx],:)=0;
        lambdaF{iFilter}(:,[1:h,sy-h+1:sy])=0; 
        %}
        gradientF{iFilter}= zeros(sx,sy);
end

    rHatNorm = 0;
    for iFilter= 1:numFilter
        rHatNorm  = rHatNorm + sum(rHat{iFilter}(:).^2);
    end
    rHatNorm = rHatNorm/overAllArea;
SSD=zeros(nIter,1);
logZRatioSeries = zeros(nIter,1);
for iter = 1:nIter
    disp( [ 'iteration: ' num2str(iter)]);
    %tic 
    [rModel, currSamples]=multiChainHMC_G(numFilter,lambdaF,filters,prevSamples,epsilon,L,numSample,nTileRow,nTileCol);
    if isComputelogZ
    % compute z ratio
     logZRatio = computeLogZRatio(prevSamples,currSamples,filters,gradientF,lambdaLearningRate,100);
     logZRatioSeries(iter)=logZRatio;
    end
    %disp(['one iteration learning time: ' num2str(toc) ' seconds']);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5   
    % compute gradient and do graidnet ascent
    gradientNorm = 0; 
    for iFilter = 1:numFilter
        gradientF{iFilter} = rHat{iFilter}-rModel{iFilter};
        h = halfFilterSizes(iFilter);
        gradientF{iFilter}([1:h,sx-h+1:sx],:)=0;
        gradientF{iFilter}(:,[1:h,sy-h+1:sy])=0;
        aa = gradientF{iFilter}; 
        gradientNorm = gradientNorm + sum(aa(:).^2);
    end
    SSD(iter)=gradientNorm/overAllArea;
    for iFilter = 1:numFilter
        lambdaF{iFilter}=lambdaF{iFilter}+ lambdaLearningRate*gradientF{iFilter};
    end
    prevSamples = currSamples;
    % visualization
    % remove the bounaries and save image
    img = currSamples;
    singleChainImg = currSamples(1:sx,1:sy);
    if isSaved
        %% saving multiple chain synthesized image 
        gLow = min(img(:));
        gHigh = max(img(:));
        disp([ 'min: ' num2str(gLow) ' max: ' num2str(gHigh)]);
        disp([ 'SSD: ' num2str(SSD(iter))]);
        disp([ 'Relative SSD: ' num2str(SSD(iter)/rHatNorm)]);
        img = (img-gLow)/(gHigh-gLow);
        imwrite(img, [savingFolder num2str(iter,'%04d') '.png']);
        %% saving the first chain of synthesized image 
        gLow=min(singleChainImg(:));
        gHigh=max(singleChainImg(:));
        singleChainImg = (singleChainImg-gLow)/(gHigh-gLow);
        imwrite(singleChainImg, [savingFolder num2str(iter,'%04d') '_singleChain.png']);
                
    end
end % iter

if isComputelogZ
% re-estimate logz
 LogZ = initialLogZ + sum(logZRatioSeries);
 disp([' Final LogZ: ' num2str(LogZ)]);
end


if isSaved
    figure;
    plot(1:nIter,SSD);
    saveas(gcf,fullfile(savingFolder,'SSD.png'));
    saveas(gcf,fullfile(savingFolder,'SSD.pdf'));
    figure;
    plot(1:nIter,SSD/rHatNorm);
    saveas(gcf,fullfile(savingFolder,'RelativeSSD.png'));
    saveas(gcf,fullfile(savingFolder,'RelativeSSD.pdf'));
  if isComputelogZ
   bar(1:nIter,logZRatioSeries);
   title(['Initial Log Z: ' num2str(initialLogZ,'%e') ', final logZ: ' num2str(LogZ,'%e')]);
   saveas(gcf,fullfile(savingFolder,'logZRatios.png'));
  end
  saveas(gcf,fullfile(savingFolder,'logZRatios.pdf'));
  disp([' The plots saved!' ]);
  close all;
end

end % end of the whole function
