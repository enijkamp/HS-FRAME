function [newfilterSelected, newxSelected,newySelected, numSelectedFeatures]= indexTtransformation(filterSelected,xSelected,ySelected, sx, sy, numOrient, numScale, numDoG)

numGabor = numOrient*numScale;

newfilterSelected=[];
newxSelected=[];
newySelected=[];


numFilter = numGabor + numDoG;
count = cell(numFilter,1);
for iFilter = 1:numFilter
    count{iFilter}=zeros(sx,sy);
end


for i=1:length(filterSelected)
    
 if (count{filterSelected(i)}(xSelected(i),ySelected(i))==0)
     
   if filterSelected(i)<=numGabor  % Gabor
        iScale=floor((filterSelected(i)-1)/numOrient)+1;
        newfilter=(iScale-1)*numOrient+ filterSelected(i);
        newfilterSelected = [newfilterSelected, newfilter, newfilter+numOrient];
        newxSelected = [newxSelected, xSelected(i), xSelected(i)];
        newySelected = [newySelected, ySelected(i), ySelected(i)];
   else   % DoG
        newfilter=numScale*numOrient+filterSelected(i);
        newfilterSelected = [newfilterSelected, newfilter];
        newxSelected = [newxSelected, xSelected(i)];
        newySelected = [newySelected, ySelected(i)];
   end
   
   count{filterSelected(i)}(xSelected(i),ySelected(i))=1;
 else
   disp(['filter:' num2str(filterSelected(i)) ' x:' num2str(xSelected(i)) ' y:' num2str(ySelected(i)) ' is selected more than once! We skip it']);  
 end
 
end


numSelectedFeatures=length(newfilterSelected);