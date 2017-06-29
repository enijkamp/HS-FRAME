function [x_added_list, y_added_list, filter_added_list, gradient_added_list, numAdded] = forwardSelection_multiple(filters, halfFilterSize, locationShiftLimit, rHat, prevSample, nRow, nCol, sx, sy, nPartCol, nPartRow, part_sx, part_sy, gradient_threshold_scale)


x_added_list=[];
y_added_list=[];
filter_added_list=[];
gradient_added_list=[];

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

for iCol=1:nPartCol
    for iRow=1:nPartRow

      previous_max=-10000;    

      for iFilter = 1:numFilter  
          
         gradientF_region=gradientF{iFilter}(1+(iRow-1)*part_sx : part_sx +(iRow-1)*part_sx, 1+(iCol-1)*part_sy : part_sy +(iCol-1)*part_sy);        
          
         [current_max, ind]= max(gradientF_region(:));     
         if(current_max > previous_max)
             previous_max = current_max;
             filter_max=iFilter;
             [x_max_part, y_max_part] = ind2sub([part_sx, part_sy],ind);       %   (x_max_part, y_max_part) are in the coordinate system of each part          
         end     
      end
      
% %       if previous_max>=threshold
% %          x_max = x_max_part +(iRow-1)*part_sx; 
% %          y_max = y_max_part +(iCol-1)*part_sy;
% %       
% %          x_added_list=[x_added_list, x_max];
% %          y_added_list=[y_added_list, y_max];
% %          filter_added_list=[filter_added_list, filter_max];     
% %          gradient_added_list=[gradient_added_list, previous_max]; 
% %       end
      
      
       x_max = x_max_part +(iRow-1)*part_sx; 
       y_max = y_max_part +(iCol-1)*part_sy;
      
       x_added_list=[x_added_list, x_max];
       y_added_list=[y_added_list, y_max];
       filter_added_list=[filter_added_list, filter_max];     
       gradient_added_list=[gradient_added_list, previous_max]; 
       
              
    end
end


addaptive_thrershold = gradient_threshold_scale * max(gradient_added_list);  
keep=find(gradient_added_list>=addaptive_thrershold);
x_added_list=x_added_list(keep);
y_added_list=y_added_list(keep);
filter_added_list=filter_added_list(keep);
gradient_added_list=gradient_added_list(keep);

numAdded=length(filter_added_list);

clear gNumChain;
clear gRSample;
clear pgRSample;
clear pgSample;

