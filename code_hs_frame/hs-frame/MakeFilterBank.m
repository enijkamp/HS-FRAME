function [ allFilterR, allFilterI, allSymbol] = MakeFilterBank(GarborScaleList, DoGScaleList, numOfOrient)

numOfGarborScale = length(GarborScaleList);
numOfDoG = length(DoGScaleList);

allOrient  = (0:numOfOrient-1) * 180/numOfOrient;
allFilterR = cell(1, numOfOrient*numOfGarborScale + numOfDoG);
allFilterI = cell(1, numOfOrient*numOfGarborScale + numOfDoG);
allSymbol  = cell(1, numOfOrient*numOfGarborScale + numOfDoG);

% Gabor filter
for s = 1:numOfGarborScale
    
    for o = 1:numOfOrient        
        n = o + (s-1)*numOfOrient;    
        [allFilterR{1, n}, allFilterI{1, n}, allSymbol{1, n}] = GaborFilter2(GarborScaleList(s), allOrient(o), s);    
       
    end;    
end;

% DoG filter
for s = 1:numOfDoG    
    n = numOfOrient*numOfGarborScale+s;
    [allFilterR{1, n}, allFilterI{1, n}, allSymbol{1, n}] = dog(DoGScaleList(s), 1.-s*.4);
    
end;
