function [allSelectedRow, allSelectedCol, allSelectedFilter]=rotateS2Template(selectedRow,selectedCol,selectedFilter, nScaleGabor, numOrient, partRotationRange, numPart)

numPartRotate = length(partRotationRange);

allSelectedRow = cell(numPart, numPartRotate); 
allSelectedCol = cell(numPart, numPartRotate);
allSelectedFilter = cell(numPart, numPartRotate);

for iPart = 1:numPart

 numElement=length(selectedRow{iPart});
   
                        
 for (r=1:numPartRotate)  
     
    rot = partRotationRange(r);
    theta = -pi*rot/numOrient;
    sintheta = sin(theta);
    costheta = cos(theta);    
     
    allSelectedRow{iPart,r}=zeros(1,numElement);
    
    if (rot == 0)
       for (i=1: numElement)
           allSelectedRow{iPart,r}(i) = selectedRow{iPart}(i);
           allSelectedCol{iPart,r}(i) = selectedCol{iPart}(i);
           allSelectedFilter{iPart,r}(i) = selectedFilter{iPart}(i); 
       end
                      
    else
       for (i=1: numElement)
           allSelectedRow{iPart,r}(i) = selectedRow{iPart}(i)*costheta + selectedCol{iPart}(i)*sintheta;
           allSelectedCol{iPart,r}(i) = -selectedRow{iPart}(i)*sintheta + selectedCol{iPart}(i)*costheta;
                
           if  selectedFilter{iPart}(i)> nScaleGabor * numOrient *2  % DoG
                   allSelectedFilter{iPart,r}(i) = selectedFilter{iPart}(i);
                
           else   % Gabor
                   iScale=floor((selectedFilter{iPart}(i)-1)/numOrient); % index of iScale starts from 0;
                   iOri= mod(selectedFilter{iPart}(i)-1,numOrient);  % index of iOri starts from 0;
                   newOri = iOri + rot;
                   if (newOri>=numOrient)
                       newOri = newOri - numOrient;
                   end
                   if (newOri<0)
                       newOri = newOri + numOrient;
                   end
                   allSelectedFilter{iPart,r}(i) = iScale*numOrient + newOri+1; % index of filter starts from 1;  
          end
       end
    end
 end
 
end 

 
        