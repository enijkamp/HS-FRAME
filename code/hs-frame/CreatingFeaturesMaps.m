%% Compute feature maps
close all; 
disp('===============> Creating the SUM1 and MAX1 maps for images');
for img = 1:numImage
    multipleResolutionImageName = [working_add '/multipleResolutionImage' num2str(img)];
    load(multipleResolutionImageName,'J','J_raw','allSizex','allSizey'); 
    
    storeFiltersBankName = [working_add '/storeFiltersBank.mat'];
    load(storeFiltersBankName, 'filters', 'halfFilterSizes', 'numFilter', 'nScaleGabor', 'nScaleDoG');
%     SUM1map = applyfilterBank_MultiResolution(J, filters, halfFilterSizes, numOrient,locationShiftLimit,orientShiftLimit,isLocalNormalize,isSeparate,...
%         localHalfx,localHalfy,localNormScaleFactor,thresholdFactor,nScaleGabor,useDoG);
    
    [SUM1map, MAX1map]= applyfilterBank_MultiResolution_sparseV4(J, filters, halfFilterSizes, numOrient,...
      locationShiftLimit,orientShiftLimit,isLocalNormalize,isSeparate,localNormScaleFactor,thresholdFactor,nScaleGabor,nScaleDoG,1);              
    %% Save the maps
    SUM1MAX1mapName = [working_add '/SUM1MAX1map' 'image' num2str(img)];   
    save(SUM1MAX1mapName, 'SUM1map','MAX1map', 'J','J_raw','allSizex','allSizey');
    
end
