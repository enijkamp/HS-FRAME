function M1 = localmax(S1, nOrient, nScaleGabor, nScaleDoG, locationShiftLimit,orientShiftLimit)
 % compute MAX1 map by SUM1 map
 
 numFilter = length(S1);
 M1 = cell(numFilter,1);    
 [sx,sy]=size(S1{1});
 
 for iFilter = 1:numFilter
     M1{iFilter} = zeros(sx,sy,'single');
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