function destImage = inscribe(rectSize,srcImage,val)
% rectSize : destination lattice
% srcImage: source image
% val : value to fill in the un-attended positions
destImage = ones(rectSize) * val;
sx = size(srcImage,1); sy = size(srcImage,2);
if size(srcImage,3) > 1
    disp('inscribe(): the source image has to be a 2D matrix');
    pause;
end
if sx * rectSize(2) >= sy * rectSize(1)
    newSx = rectSize(1); newSy = floor(.5+newSx*sy/sx );
    srcImage = imresize(srcImage,[newSx,newSy],'nearest');
    destImage(:,(1:newSy)+floor((rectSize(2)-newSy)/2)) = srcImage(:,:);
else
    newSy = rectSize(2); newSx = floor(.5+newSy*sx/sy );
    srcImage = imresize(srcImage,[newSx,newSy],'nearest');
    destImage((1:newSx)+floor((rectSize(1)-newSx)/2),:) = srcImage(:,:);
end
