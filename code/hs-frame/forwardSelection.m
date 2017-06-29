function [x_added, y_added, filter_added] = forwardSelection(filters, halfFilterSize, locationShiftLimit, rHat, prevSample, nRow, nCol, sx, sy)

numFilter = length(filters);
rModel = cell(numFilter,1);
gNumChain =gpuArray(nRow*nCol);

gFilters = cell(size(filters));
for iEl = 1:numel(filters)
   gFilters{iEl}=gpuArray(filters{iEl});
end
pgSample = gpuArray(prevSample);

for iFilter = 1:numFilter
      pgRSample=abs(filter2(gFilters{iFilter},pgSample));
      gRSample =gpuArray.zeros(sx,sy);         % If you use a new version of Matlab, please use this statement
      %gRSample = parallel.gpu.GPUArray.zeros(sx,sy);     % If you use an old version of Matlab, please use this statement
      for iRow = 1:nRow
        for iCol = 1:nCol
            gRSample = gRSample + pgRSample((iRow-1)*sx+1:iRow*sx,(iCol-1)*sy+1:iCol*sy);
        end
      end
      gRSample = gRSample/gNumChain;
      rModel{iFilter}=gather(gRSample);
end

gradientF = cell(numFilter,1);
for iFilter = 1:numFilter
    gradientF{iFilter}= zeros(sx, sy);
end

for iFilter = 1:numFilter
    gradientF{iFilter} = rHat{iFilter}-rModel{iFilter};
    %% free boundary
    h = halfFilterSize(iFilter);
    gradientF{iFilter}([1:h+locationShiftLimit,sx-h-locationShiftLimit+1:sx],:)=0;
    gradientF{iFilter}(:,[1:h+locationShiftLimit,sy-h-locationShiftLimit+1:sy])=0;
        
end

previous_max=-10000;
tic
for iFilter = 1:numFilter
    
     [current_max, ind]= max(gradientF{iFilter}(:));
     
     if(current_max > previous_max)
         previous_max = current_max;
         filter_max=iFilter;
         [x_max, y_max] = ind2sub([sx, sy],ind);                  
     end     
end
disp(['xxxxxxxxxxx', num2str(toc)]);

x_added = x_max;
y_added = y_max;
filter_added = filter_max;

