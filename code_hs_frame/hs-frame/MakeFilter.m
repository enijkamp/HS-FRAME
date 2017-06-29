function [allFilter, allSymbol] = MakeFilter(scaleFilter, numOrient) 
% make Gabor filters at fixed scale with numOrient orientations
% "allFilter" contains all the Gabor sine and cosine pairs
% "allSymbol" contains bars for display purpose
allOrient  = (0:numOrient-1) * 180/numOrient; 
allFilter = cell(1, numOrient); 
allSymbol = cell(1, numOrient); 
for o = 1 : numOrient
     [allFilter{1, o}, allSymbol{1, o}] = GaborFilter(scaleFilter, allOrient(o));
end
    



    

