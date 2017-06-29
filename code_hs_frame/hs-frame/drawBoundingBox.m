         
showObjectBoundingBox=1;
showPartBoundingBox=1;


%% set rgb format of grayscale image as background of "matchBoundingBox".    
tmp_img=imageLoaded.imageOriginal;
tmp_img=imresize(tmp_img, size(imageLoaded.ImageMultiResolution{Mind}), 'nearest');
greyImg = rgb2gray(tmp_img);  % images at multiple resolutions    
           
            
rgb_greyImg = cat(3,greyImg,greyImg,greyImg);
matchedBoundingBox = rgb_greyImg;
            


%% prepare the output variable for visualization of matched template
if showObjectBoundingBox
    margin = 2;  % line width of the bounding box
    xx = repmat((1:sx),1, margin*2);
    yy = [];
    for y = [1:margin sy-margin+1:sy]
        yy = [yy,ones(1,sx)*y];
    end
    yy = [yy,repmat((1:sy),1,margin*2)];
    for x = [1:margin sx-margin+1:sx]
        xx = [xx,ones(1,sy)*x];
    end
    inRow = single(xx-floor(sx/2));
    inCol = single(yy-floor(sy/2));
    tScale = 0; rScale = 1; cScale = 1; inO = zeros(numel(inRow),1,'single'); inS = zeros(numel(inRow),1,'single');
    
    [outRow, outCol] = mexc_TemplateAffineTransform(tScale,rScale,cScale,Mrot,inRow,inCol,inO,inS,nOrient);
    
    % directly overwrite the corresponding pixels
    for p = 1:length(outRow)
        x = round(MFx + outRow(p));
        y = round(MFy + outCol(p));
        if x > 0 && x <= size(matchedBoundingBox,1) && y > 0 && y <= size(matchedBoundingBox,2)
            matchedBoundingBox(x,y,:) = [255 0 0];  % color of the bounding box
            % matchedBoundingBox(x,y) = [0.1];
        end
    end
end


%%% draw bounding box for parts
if showPartBoundingBox        
     for iPart = 1:numPart      % draw each bounding box for parts
	      
       r = clusters(c).S3T{bestRotInd}.selectedTransform(iPart)+1; % the transform index now starts from 1 not 0  (rotation) 
       Fx = MFx + round(clusters(c).S3T{bestRotInd}.selectedRow(iPart));
       Fy = MFy + round(clusters(c).S3T{bestRotInd}.selectedCol(iPart));
       
       imagesize = size(MAX2map{r,iPart,Mind});
       
       % set default values of some output variables
       bestPartRes = Mind;
       
       if Fx >= 1 && Fx <= imagesize(1) && Fy >= 1 && Fy <= imagesize(2)
           tmp = MAX2map{r,iPart,Mind};
           
           tmp = MAX2ResolutionTrace{r,iPart,Mind};
           bestPartRes = tmp(Fx,Fy) + 1; % best part resolution
           current_size = size(tmp);
           tmp = MAX2LocTrace{r,iPart,bestPartRes};
           new_size = size(tmp);
           
           Fx = round(Fx*double(new_size)/current_size);
           Fy = round(Fy*double(new_size)/current_size);
           
           if Fx >= 1 && Fx <= size(tmp,1) && Fy >= 1 && Fy <= size(tmp,2);
               translationInd = tmp(Fx,Fy) + 1;
           else
               translationInd = floor(size(M2RowColShift,1)/2);
           end
           
           tmp = MAX2TransformTrace{r,iPart,bestPartRes};
           if Fx >= 1 && Fx <= size(tmp,1) && Fy >= 1 && Fy <= size(tmp,2);
               transformInd = tmp(Fx,Fy) + 1;
           else
               transformInd = floor(numPartRotate/2) + 1;
           end
           
           actualPartRotationInd = transformInd - numPartRotate*(ceil(double(transformInd)/numPartRotate)-1);
           
           Fx = floor( Fx + M2RowColShift(translationInd,1));
           Fy = floor( Fy + M2RowColShift(translationInd,2));
           
       else
           actualPartRotationInd = r;
       end
       actualPartRotation = partRotationRange(actualPartRotationInd);
       
       % find the part location at the higher resolution
       %           Fx = (Fx-1 + .5);
       %           Fy = (Fy-1 + .5);
       
       % show Part Bounding Box
       %tmpMatchedSym = double( ImageMultiResolution{Mind} );
       %tmpMatchedSym = double( imresize(tmpMatchedSym,imageSizeAtBestObjectResolution,'bilinear') );
       %matchedSym = max(matchedSym,tmpMatchedSym);
       
       margin = 2;
       xx = repmat((1:part_sx),1,margin*2);
       yy = [];
       for y = [1:margin part_sy-margin+1:part_sy]
           yy = [yy,ones(1,part_sx)*y];
       end
       yy = [yy,repmat((1:part_sy),1,margin*2)];
       for x = [1:margin part_sx-margin+1:part_sx]
           xx = [xx,ones(1,part_sy)*x];
       end
       inRow = single(xx-floor(part_sx/2)); inCol = single(yy-floor(part_sy/2));
       tScale = 0; rScale = 1; cScale = 1; inO = zeros(numel(inRow),1,'single'); inS = zeros(numel(inRow),1,'single');
       [outRow, outCol] = mexc_TemplateAffineTransform(tScale,rScale,cScale,actualPartRotation,inRow,inCol,inO,inS,nOrient);
       
       % directly overwrite the corresponding pixels
       matchedBoundingBox = imresize(matchedBoundingBox,size(imageLoaded.ImageMultiResolution{bestPartRes}),'nearest');
       for p = 1:length(outRow)
           x = round(outRow(p) + Fx);
           y = round(outCol(p) + Fy);
           if x > 0 && x <= size(matchedBoundingBox,1) && y > 0 && y <= size(matchedBoundingBox,2)
               matchedBoundingBox(x,y,:) = [0 0 255];
           end
       end
       matchedBoundingBox = imresize(matchedBoundingBox,size(imageLoaded.ImageMultiResolution{Mind}),'nearest');
     end
     
end            
            
imwrite(matchedBoundingBox, fullfile(boundingBoxSavingFolder,['boundingBox-cluster-' num2str(c) '-img-' num2str(iImg,'%04d') '.png']));  
% 
            