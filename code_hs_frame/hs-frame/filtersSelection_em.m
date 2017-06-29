% clear all;
% close all;

% filtersSelectionV3.m is for outputing the deformed templates for each
% training images


function [template, deformedTemplate] = filtersSelection_em(I, GarborScaleList, DoGScaleList, numOrient, numSketch, numTopFeatureToShow, locationShiftLimit, orientShiftLimit, sx, sy,savingFolder)

%% set parameters 

numScale = length(GarborScaleList); %number of scales for Gabor wavelets
numDoG = length(DoGScaleList);

epsilon = .1; 
Upperbound = 6.; 
SHUTUP = 0.; 

numImg = length(I);

%% read training images
% imgFold = 'Image';
% imgName = dir([imgFold '/*.jpg']);
% numImg = size(imgName, 1);
% clear I;
% for (i = 1 : numImg)
%     img = imread([imgFold '/' imgName(i).name]);
%     img = imresize(img, [sx,sy],'nearest');
%     img = im2double(rgb2gray(img));
%     img=img-mean(img(:));
%     img=img/std(img(:));
%     I{i} = img;
% end


%% filter images
tic;
[allFilterR, allFilterI, allSymbol] = MakeFilterBank(GarborScaleList, DoGScaleList, numOrient);
disp(['finished generating filter bank: '  num2str(toc) ' seconds']);

tic;
[Crr, Cri, Cir, Cii] = CorrFilter(allFilterR, allFilterI); %calculate correlation among filters
disp(['finished computing filters correlation: ' num2str(toc) ' seconds']);

tic;
[allFilteredR, allFilteredI, allFiltered] = ApplyFilter(I, allFilterR, allFilterI);
disp(['finished images filtering: ' num2str(toc) ' seconds']);

half = zeros(1, numOrient*numScale+numDoG); 
for (i=1:numOrient*numScale+numDoG)
    half(i) = (size(allFilterR{1, i}, 1)-1)/2;
end

%% prepare for storing experimental results
% multiple-scaled Gobar template
syma = zeros(sx, sy); 
syma_top = zeros(sx, sy); 
if numDoG>0
    numDoGTemplate=1;
else
    numDoGTemplate=0;
end
  
sym = cell(1, numScale + numDoGTemplate); 
for (i = 1 : numScale + numDoGTemplate)
   sym{i} = zeros(sx, sy); 
end
    
% recontructed images
Asym = cell(1, numImg); 
for (i = 1 : numImg)
    Asym{i} = zeros(sx, sy); 
end

% deformed Gabor templates
Asyma = cell(1, numImg); 
for (i = 1 : numImg)
    Asyma{i} = zeros(sx, sy); 
end

% deformed Gabor templates visualized by only a few of top features
Asyma_top = cell(1, numImg); 
for (i = 1 : numImg)
    Asyma_top{i} = zeros(sx, sy); 
end

% deformed DoG templates
Asym0 = cell(1, numImg*(numScale+1));
for (i = 1 : numImg*(numScale+1))
    Asym0{i} = zeros(sx, sy);  
end

filterSelected = zeros(1, numSketch);
xSelected = zeros(1, numSketch);
ySelected = zeros(1, numSketch);
%numSketchSelected = 0;

deformed_filterSelected=zeros(numImg,numSketch);
deformed_xSelected=zeros(numImg,numSketch);
deformed_ySelected=zeros(numImg,numSketch);

%% do matching pursuit
tic;
MatchingPursuitV3(numImg, numOrient*numScale+numDoG, numOrient, numDoG, allFilteredR,...
                allFilteredI, allFiltered, Crr, Cri, Cir, Cii, half, sx, sy,...
                allFilterR(1, :), allFilterI(1, :), allSymbol(1, :), syma, syma_top, sym,...
                Asym, Asym0, Asyma, Asyma_top, locationShiftLimit, orientShiftLimit, numSketch,...
                epsilon, Upperbound, SHUTUP, filterSelected, xSelected, ySelected, ...
                deformed_filterSelected, deformed_xSelected, deformed_ySelected, numTopFeatureToShow);
 
disp(['finished Matching Pursuit: ' num2str(toc) ' seconds']);
  
 numSketchSelected=length(filterSelected);
 %% deformed template for every training image          
 filterSelected = filterSelected(1:numSketchSelected); 
 xSelected = xSelected(1:numSketchSelected);
 ySelected = ySelected(1:numSketchSelected);
 
 
 deformed_filterSelected = deformed_filterSelected(:,1:numSketchSelected);
 deformed_xSelected = deformed_xSelected(:,1:numSketchSelected);
 deformed_ySelected = deformed_ySelected(:,1:numSketchSelected);
 
 deformedTemplate=cell(1,numImg);
 for i=1:numImg
     [deformed_filter, deformed_x, deformed_y, num]= ...
         indexTtransformation(deformed_filterSelected(i,:), deformed_xSelected(i,:), deformed_ySelected(i,:), sx, sy, numOrient, numScale, numDoG);
     
     deformedTemplate{i}=struct('selectedFilter',single(deformed_filter),'selectedRow',single(deformed_x),'selectedCol',single(deformed_y),'numSelected',single(num));
 end
 %%
 

 %[indEliminated, numSelectedFeatures]= createMask(filterSelected,xSelected,ySelected, sx, sy, numOrient, numScale, numDoG);
 [filterSelected, xSelected,ySelected, numSelectedFeatures]= indexTtransformation(filterSelected, xSelected, ySelected, sx, sy, numOrient, numScale, numDoG);

 template=struct('selectedFilter',single(filterSelected),'selectedRow',single(xSelected),'selectedCol',single(ySelected), 'numSelected',single(numSelectedFeatures));
% show and save experimental results
ShowAndSaveResults;

%% generate html to show filter selection experimental results
GenerateHtmlForFilterSelection;




