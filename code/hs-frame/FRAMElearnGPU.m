function [lambdaF,currSamples, LogZ]=FRAMElearnGPU(nIter,filters,rHat,numFilter,sx,sy,halfFilterSizes,nTileRow,nTileCol,epsilon,L,lambdaLearningRate,numSample,isSaved,savingFolder,isComputelogZ)

overAllArea = sum((sx-2*halfFilterSizes).*(sy-2*halfFilterSizes));
lambdaF = cell(numFilter,1);
gradientF = cell(numFilter,1);
for iFilter = 1:numFilter
    lambdaF{iFilter} = zeros(sx,sy);
    % the following lines will be useful when lambda is not initialized as zero;
    %{
    h = halfFilterSizes(iFilter);
    lambdaF{iFilter}([1:h,sx-h+1:sx],:)=0;
    lambdaF{iFilter}(:,[1:h,sy-h+1:sy])=0; 
    %}
    gradientF{iFilter}= zeros(sx,sy);
end

prevSamples = randn(sx*nTileRow,sy*nTileCol);
initialLogZ =log((2*pi))*(overAllArea/2);
SSD=zeros(nIter,1);
logZRatioSeries = zeros(nIter,1);
for iter = 1:nIter
    disp( [ 'iteration: ' num2str(iter)]);
    tic 
    
    %[rModel, currSamples]=multiChainHMC(numFilter,lambdaF,filters,prevSamples,0.01,10,nTileRow,nTileCol);
    [rModel, currSamples]=multiChainHMC_G(numFilter,lambdaF,filters,prevSamples,epsilon,L,numSample,nTileRow,nTileCol);
    if isComputelogZ
    % compute z ratio
     logZRatio = computeLogZRatio(prevSamples,currSamples,filters,gradientF,lambdaLearningRate,100);
     logZRatioSeries(iter)=logZRatio;
    end
    %
    disp(['one iteration learning time: ' num2str(toc) ' seconds']);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5   
    % compute gradient and do graidnet ascent
    gradientNorm = 0; 
    for iFilter = 1:numFilter
        gradientF{iFilter} = rHat{iFilter}-rModel{iFilter};
        h = halfFilterSizes(iFilter);
        gradientF{iFilter}([1:h,sx-h+1:sx],:)=0;
        gradientF{iFilter}(:,[1:h,sy-h+1:sy])=0;
        aa = gradientF{iFilter}; 
        gradientNorm = gradientNorm + sum(abs(aa(:)));
    end
    SSD(iter)=gradientNorm/overAllArea;
    %lambdaLearningRate = 1e-3%gradientNorm;
    for iFilter = 1:numFilter
        lambdaF{iFilter}=lambdaF{iFilter}+ lambdaLearningRate*gradientF{iFilter};
    end
    prevSamples = currSamples;
    % visualization
    % remove the bounaries and save image
    img = currSamples;
    %{
    h = 5;
    img(1:h,:)=0; img(:,1:h)=0;
    for iRow = 1:nTileRow-1
       img(iRow*sx-h:iRow*sx+h+1,:)=0;
    end
    for iCol = 1:nTileCol-1
       img(:,iCol*sy-h:iCol*sy+h+1)=0;
    end
    img(end-h:end,:)=0;
    img(:,end-h:end)=0;
    %}
    if isSaved
    gLow = min(img(:));
    gHigh = max(img(:));
    disp([ 'min: ' num2str(gLow) ' max: ' num2str(gHigh)]);
    disp([ 'SSD: ' num2str(SSD(iter))]);
    %disp([ 'Relative SSD: ' num2str(SSD(iter)/rHatNorm)]);
    img = (img-gLow)/(gHigh-gLow);
    imwrite(img, [savingFolder num2str(iter,'%04d') '.png']);
    end
end

if isComputelogZ
% re-estimate logz
 LogZ = initialLogZ + sum(logZRatioSeries);
 disp([' Final LogZ: ' num2str(LogZ)]);
end


if isSaved
  figure;
  plot(1:nIter,SSD);
  saveas(gcf,fullfile(savingFolder,'SSD.png'));
%  saveas(gcf,fullfile(outPath,'SSD.pdf'));
%   figure;
%   plot(1:nIter,SSD/rHatNorm);
%   saveas(gcf,fullfile(outPath,'RelativeSSD.png'));
%   saveas(gcf,fullfile(outPath,'RelativeSSD.pdf'));
  if isComputelogZ
   bar(1:nIter,logZRatioSeries);
   title(['Initial Log Z: ' num2str(initialLogZ,'%e') ', final logZ: ' num2str(LogZ,'%e')]);
   saveas(gcf,fullfile(savingFolder,'logZRatios.png'));
  end
%  saveas(gcf,fullfile(outPath,'logZRatios.pdf'));
  disp([' The plots saved!' ]);
  close all;
  
end
