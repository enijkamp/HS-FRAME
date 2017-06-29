function [allFilteredSUM1, allFilteredMAX1] = applyfilterBank_MultiResolution_sparseV4(I, filters, halfFilterSize, nOrient,locationShiftLimit,orientShiftLimit,LocalNormOrNot,isSeparate,scaleFactor,thresholdFactor,nScaleGabor,nScaleDoG, scale)
% filter image by a bank of filters
% I: input images
% allFilter: filter bank
numImage = size(I, 2);    % number of images
numFilter = size(filters, 2);   % number of orientations
allFilteredMAX1 = cell(numImage, numFilter);  % filtered images
allFilteredSUM1 = cell(numImage, numFilter);  % filtered images

for iImg = 1:numImage
    S1 = cell(numFilter,1);
    M1 = cell(numFilter,1);
    
    % SUM1
    
    [sx,sy]=size(I{iImg});
    for iFilter = 1:numFilter
                
        Y = filter2(filters{iFilter},I{iImg});
        S1{iFilter} = abs(single(Y));
        M1{iFilter} = zeros(sx,sy,'single');
    end
       
    
    
    if LocalNormOrNot
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Local Normalization
    
    if isSeparate        
      for iScale=1:nScaleGabor  
        h = halfFilterSize(1+(iScale-1)*(2*nOrient)); 
        [S1(1+(iScale-1)*(2*nOrient) : nOrient+(iScale-1)*(2*nOrient))]= LocalNormalizeV3(S1(1+(iScale-1)*(2*nOrient):nOrient+(iScale-1)*(2*nOrient)),[],h,round(scaleFactor*h),round(scaleFactor*h),thresholdFactor);
        [S1(nOrient+1+(iScale-1)*(2*nOrient): 2*nOrient+(iScale-1)*(2*nOrient))] = LocalNormalizeV3(S1(nOrient+1+(iScale-1)*(2*nOrient): 2*nOrient+(iScale-1)*(2*nOrient)),[],h,round(scaleFactor*h),round(scaleFactor*h),thresholdFactor);
      end
      % for DoG
      for iScale=1:nScaleDoG  
          ind  = 2*nOrient*nScaleGabor+iScale;
          h = halfFilterSize(ind);       
          [S1(ind)]=LocalNormalizeV3(S1(ind),[],h,round(scaleFactor*h),round(scaleFactor*h),thresholdFactor);                
      end
      
    else
        for iScale=1:nScaleGabor
           h = halfFilterSize(1+(iScale-1)*(2*nOrient)); 
           [S1(1+(iScale-1)*(2*nOrient):nOrient+(iScale-1)*(2*nOrient)),S1(nOrient+1+(iScale-1)*(2*nOrient): 2*nOrient+(iScale-1)*(2*nOrient))]= LocalNormalizeV3(S1(1+(iScale-1)*(2*nOrient):nOrient+(iScale-1)*(2*nOrient)),S1(nOrient+1+(iScale-1)*(2*nOrient): 2*nOrient+(iScale-1)*(2*nOrient)),h,round(scaleFactor*h),round(scaleFactor*h),thresholdFactor);
        end
        % for DoG
        for iScale=1:nScaleDoG 
           ind  = 2*nOrient*nScaleGabor+iScale;
           h = halfFilterSize(ind);       
           [S1(ind)]=LocalNormalizeV3(S1(ind),[],h,round(scaleFactor*h),round(scaleFactor*h),thresholdFactor);                
        end        
    end 
    
    %%%%%%%
    for iFilter = 1:numFilter 
        S1{iFilter }=S1{iFilter}*scale;    
    end
     
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    end
        
    
    %MAX1
    for iScale=1:nScaleGabor
      CgetMAX1(1,sx,sy,nOrient,locationShiftLimit,orientShiftLimit,1,S1(1+(iScale-1)*(2*nOrient):nOrient+(iScale-1)*(2*nOrient)),M1(1+(iScale-1)*(2*nOrient):nOrient+(iScale-1)*(2*nOrient)));
      CgetMAX1(1,sx,sy,nOrient,locationShiftLimit,orientShiftLimit,1,S1(nOrient+1+(iScale-1)*(2*nOrient): 2*nOrient+(iScale-1)*(2*nOrient)),M1(nOrient+1+(iScale-1)*(2*nOrient): 2*nOrient+(iScale-1)*(2*nOrient))); 
 
    end
    % DoG
    for iScale=1:nScaleDoG
       ind  = 2*nOrient*nScaleGabor+iScale;       
       CgetMAX1_DoG(1,sx,sy,locationShiftLimit,1,S1(ind),M1(ind)); 
  
    end
       
    
    allFilteredMAX1(iImg,:)=M1;
    allFilteredSUM1(iImg,:)=S1;
end






