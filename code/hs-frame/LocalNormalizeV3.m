function [cosinPart, sinePart]=LocalNormalizeV3(cosinPart,sinePart,h,halfSizex,halfSizey,thresholdFactor)
% function to mimic the mexc_localNormalize by Zhangzhang Si.
% suppose for the input image, pixels within h pixels distance from the
% boundar are all 0.
% h is the half of filter size
% sine and cosine parts are actually aboslute value of <I,B>

% if(isempty(sinePart))
%     disp('sinPart is empty, will continue')
% end

sx = size(cosinPart{1},1);
sy = size(cosinPart{1},2);

% compute sum of sum1 over all orientations
S1All = zeros([sx,sy]);
if isempty(sinePart)     % for DoG
    for iOri = 1:length(cosinPart)
        S1All = S1All + cosinPart{iOri}.^2;
    end
else
    for iOri = 1:length(cosinPart)
        S1All =S1All+  0.5.*(cosinPart{iOri}.^2+sinePart{iOri}.^2);
    end
end
S1All = S1All/length(cosinPart);


%integralImage
meanFilter = ones(2*halfSizex+1,2*halfSizey+1);

meanMap = filter2(meanFilter,S1All);
meanMap = sqrt(meanMap/numel(meanFilter));


% % % propagate local average to boundary
% % if ((h+1+halfSizex>=sx-h-halfSizex) || (h+1+halfSizey>=sy-h-halfSizey))
% %    error('Filter size or Halfsize for Local Normalization is too large! Please try smaller filters.'); 
% % end
% % 
% % meanMap = meanMap(h+1+halfSizex:sx-h-halfSizex,h+1+halfSizey:sy-h-halfSizey);
% % meanMap = padarray(meanMap,[h+halfSizex,h+halfSizey],'replicate','both');


% perform normalization
maxAverage = max(meanMap(:));
meanMap = max(meanMap,maxAverage*thresholdFactor);
invMeanMap = single(1./meanMap);
invMeanMap([1:h end-h:end],:)=0;
invMeanMap(:,[1:h end-h:end])=0;

if isempty(sinePart)
    for iOri = 1:length(cosinPart)
        cosinPart{iOri}=cosinPart{iOri}.*invMeanMap;
    end
else
    for iOri = 1:length(sinePart)
       sinePart{iOri} = sinePart{iOri}.*invMeanMap;
       cosinPart{iOri} =cosinPart{iOri}.*invMeanMap;
    end
end
