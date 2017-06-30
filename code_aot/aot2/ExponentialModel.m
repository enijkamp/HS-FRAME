clear; close all; 
mex ClocalNormalizeDouble.c; % local normalization of type double
mex ClocalNormalize.c; % local normalization
%% parameters 
scaleFilter = .7;  % scale of Gabors, length = 17 pixels
numOrient = 16;  % number of orientations
epsilon = .1;
saturation = 6.; % saturation level for sigmoid transformation
binSize = .2;  % binsize for computing histogram of q()
%% parameters for normalization
localHalfx = floor(30*scaleFilter+.5); localHalfy = floor(30*scaleFilter+.5); % the half range for local normalization, has to be quite large
thresholdFactor = .01;  % divide the response by max(average, maxAverage*thresholdFactor)
%% read in background images 
imageName = dir('backgroundImage/*.jpg');  % folder that contains background images as q(I)
numImage = size(imageName, 1); % number of background images
imageSize = zeros(numImage, 2); % size of background images
for i = 1 : numImage
    tmpIm = imread(['backgroundImage' '/' imageName(i).name]);
    if size(tmpIm,3) == 3
        tmpIm = rgb2gray(tmpIm);
    end
    I{i} = imresize(double(tmpIm), 1, 'nearest'); 
    imageSize(i, :) = size(I{i});
end
sizex = min(imageSize(:, 1)); sizey = min(imageSize(:, 2)); 
for i = 1 : numImage
    I{i} = I{i}(1:sizex, 1:sizey); 
end   % make the sizes of the images to be the same
%% filtering background images
disp(['start filtering']);
tic
[allFilter, allSymbol] = MakeFilter(scaleFilter, numOrient);  % generate Gabor filters 
halfFilterSize = (size(allFilter{1}, 1)-1)/2;  % half size of Gabor
allFilteredImage = ApplyFilterfftSame(I, allFilter, localHalfx, localHalfy, 1, thresholdFactor);   % filter training images
disp(['filtering time: ' num2str(toc) ' seconds']);
%% compute histogram of q()
numBin = floor(saturation/binSize)+1;  % binnumbers
histog = zeros(numBin, 1);  % store F  
mex Chistogram.c;   % compile C code
disp(['start histogramming']);
tic
Chistogram(numImage, numOrient, allFilteredImage, halfFilterSize, sizex, sizey, binSize, numBin, histog, saturation);
disp(['histogramming time: ' num2str(toc) ' seconds']);
%% compute stored lambda, expectation, logZ
r = (0:(numBin-1))*binSize;
figure; plot(r, histog); 
title('background density'); xlabel('sigmoid(response)'); ylabel('density');
numStoredPoint = 50; 
storedExpectation = zeros(numStoredPoint, 1); 
Z = zeros(numStoredPoint, 1); 
for (k=1:numStoredPoint)
   lambda = (k-1.)/10.; 
   p = exp(lambda*r).*(histog'); 
   Z(k) = sum(p*binSize); p = p/Z(k); 
   storedExpectation(k) = sum(r.*p*binSize);
end
storedlambda = (0:(numStoredPoint-1))/10.; 
storedLogZ = log(Z); 
figure; plot(storedlambda, storedExpectation);
title('mean versus lambda'); xlabel('lambda'); ylabel('mean');
figure; plot(storedlambda, storedLogZ);
title('logZ versus lambda'); xlabel('lambda'); ylabel('logZ');

Correlation = CorrFilter(allFilter, epsilon);  % correlation between filters 

save 'storedExponentialModel' numStoredPoint storedlambda storedExpectation storedLogZ ...
      saturation scaleFilter numOrient halfFilterSize allFilter allSymbol ...
      localHalfx localHalfy thresholdFactor Correlation; 
  




