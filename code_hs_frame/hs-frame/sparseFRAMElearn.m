function [template,currSample, LogZ]=sparseFRAMElearn(template, nIteration,filters,rHat,sx,sy,halfFilterSize,locationShiftLimit, nTileRow, nTileCol, epsilon,...
    L,lambdaLearningRate,numSample, nPartCol, nPartRow, part_sx, part_sy, isSaved,savingFolder)



nValidFeature= sum((sx-2*halfFilterSize-2*locationShiftLimit).*(sy-2*halfFilterSize-2*locationShiftLimit));
SSD=zeros(nIteration,1);
%initialLogZ = log((2*pi))*(nValidFeature/2);
initialLogZ =0 ;
logZRatioSeries = zeros(nIteration,1);


% the field of enlargeTemplate stores a big tiled template for multiple chains HMC sampling
template.enlargedTemplate = enlargingTemplate(nTileRow, nTileCol, sx, sy, template);
prevSample = single(randn(sx*nTileRow,sy*nTileCol));

rHat_selected = zeros(1,template.numSelected);
for i=1:template.numSelected
    rHat_selected(i)=rHat{template.selectedFilter(i)}(template.selectedRow(i),template.selectedCol(i));
end

gradientF_selected=zeros(1,template.numSelected);
for iter = 1:nIteration
    disp( [ 'iteration: ' num2str(iter)]);
    tic
    
    [rModel_selected, currSample]=multiChainHMC_sparse(template, filters, prevSample, epsilon, L, numSample, nTileRow, nTileCol);
    
    % compute z ratio
    logZRatio = computeLogZRatio_sparse(prevSample,currSample,template,filters,gradientF_selected,lambdaLearningRate, nTileRow, nTileCol);
    logZRatioSeries(iter)=logZRatio;
    %
    
    disp(['one iteration learning time: ' num2str(toc) ' seconds']);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
    % compute gradient and do graidnet ascent
    
    gradientF_selected = rHat_selected-rModel_selected;
    
    
    SSD(iter)=mean(abs(gradientF_selected));
    
    template.selectedLambdas = template.selectedLambdas + lambdaLearningRate * gradientF_selected;
    
    prevSample=currSample;
    
    if isSaved
        % save synthesied image
        gLow = min(currSample(:));
        gHigh = max(currSample(:));
        disp([ 'min: ' num2str(gLow) ' max: ' num2str(gHigh)]);
        disp([ 'SSD: ' num2str(SSD(iter))]);
        img = (currSample-gLow)/(gHigh-gLow);
        imwrite(img, [savingFolder num2str(iter,'%04d') '.png']);
    end
    
end

% re-estimate logz
LogZ = initialLogZ + sum(logZRatioSeries);
disp([' Final LogZ: ' num2str(LogZ)]);

%%%save the multiple chains of (projected) synthesized images separately
%% solution of least squares min |I_syn - B_selected * C|^2
%% C = inv(B'B)B'I;  I_reconstruct = B*C = B* inv(B'B)B' I_syn, we compute and store R=B inv(B'B)B' as reconstruction matrix
B=[];
for i=1:template.numSelected
    Bi = filterTiling(filters{template.selectedFilter(i)}, sx, sy, template.selectedRow(i), template.selectedCol(i));
    B = [B,Bi(:)];
end
R=B*(inv(B'*B)*(B'));
%%%% save each chain of synthesized images respectively in the last iteration
EachChainSynFolder=[savingFolder 'EachChainSyn/'];
if exist(EachChainSynFolder)
    rmdir(EachChainSynFolder,'s');
    mkdir(EachChainSynFolder);
else
    mkdir(EachChainSynFolder);
end
%%%% reconstruct each chain of synthesized images respectively in the last iteration
EachChainReconstructionFolder=[savingFolder 'EachChainReconSyn/'];
if exist(EachChainReconstructionFolder)
    rmdir(EachChainReconstructionFolder,'s');
    mkdir(EachChainReconstructionFolder);
else
    mkdir(EachChainReconstructionFolder);
end
curr_projected=zeros(nTileRow*sx,nTileCol*sy);
currSample2=zeros(nTileRow*sx,nTileCol*sy);
index=0;
for j=1:nTileCol
    for i=1:nTileRow
        index=index+1;
        crop=currSample(1+(i-1)*sx:sx+(i-1)*sx, 1+(j-1)*sy:sy+(j-1)*sy);
        gLow = min(crop(:));
        gHigh = max(crop(:));
        img = (crop-gLow)/(gHigh-gLow);
        imwrite(img,[EachChainSynFolder 'original_chain-' num2str(index,'%04d') '.png']);
        currSample2(1+(i-1)*sx:sx+(i-1)*sx, 1+(j-1)*sy:sy+(j-1)*sy)=img;
        
        %% project the synthesized images onto the domain specified by the selected basis
        img_projected = reshape(R*crop(:),sx,sy);  % R = B*(inv(B'*B)*(B'));
        curr_projected(1+(i-1)*sx:sx+(i-1)*sx, 1+(j-1)*sy:sy+(j-1)*sy)=img_projected;
        gLow = min(img_projected(:));
        gHigh = max(img_projected(:));
        img_projected = (img_projected-gLow)/(gHigh-gLow);
        imwrite(img_projected,[EachChainReconstructionFolder 'reconstruct_chain-' num2str(index,'%04d') '.png']);
        curr_projected(1+(i-1)*sx:sx+(i-1)*sx, 1+(j-1)*sy:sy+(j-1)*sy)=img_projected;
    end
end
imwrite(currSample2,[EachChainSynFolder 'raw.png']);
imwrite(curr_projected,[EachChainReconstructionFolder 'reconstruct.png']);

%%%%%%
% save parts
PartFolder=[savingFolder 'PartSyn/'];
if exist(PartFolder)
    rmdir(PartFolder,'s'); 
    mkdir(PartFolder);
else
    mkdir(PartFolder);
end

index=0;
for i=1:nPartRow
    for j=1:nPartCol
        index=index+1;
        img_parts=img_projected(1+(i-1)*part_sx:part_sx+(i-1)*part_sx, 1+(j-1)*part_sy:part_sy+(j-1)*part_sy);
        gLow = min(img_parts(:));
        gHigh = max(img_parts(:));
        img_parts = (img_parts-gLow)/(gHigh-gLow);  
        imwrite(img_parts,[PartFolder 'part-' num2str(index,'%04d') '.png']);
    end
end



% showing how many lambdas selected.
A=template.selectedLambdas;
numPosLambdas=length(A(A>0));
numNegLambdas=length(A(A<0));
numZeroLambdas=length(A(A==0));
disp([ num2str(numPosLambdas) ' positive lambdas left,' num2str(numNegLambdas) ' negative lambdas left and ' ...
    num2str(numZeroLambdas) ' zero lambdas left.']);



if isSaved
    h = figure;
    plot(1:nIteration,SSD);
    title('Convergence');
    xlabel('Iterations');
    ylabel('SSD');
    saveas(h,[savingFolder 'plot.png'],'png');
    
    disp([' Final LogZ: ' num2str(LogZ)]);
    bar(1:nIteration,logZRatioSeries);
    title(['Initial Log Z: ' num2str(initialLogZ,'%e') ', final logZ: ' num2str(LogZ,'%e')]);
    saveas(gcf,fullfile(savingFolder,'logZRatios.png'));
    disp([' The plots saved!' ]);
    close all;
end



