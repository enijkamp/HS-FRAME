function [resultR, resultI, resultE] = ApplyFilter(I, allFilterR, allFilterI)

nimage = size(I, 2); 
[sx, sy] = size(I{1}); 
no = size(allFilterR, 2);
resultR = cell(nimage, no); 
resultI = cell(nimage, no); 
resultE = cell(nimage, no); 
adjustedI = cell(nimage, 1); 

for i = 1 : nimage
   [sx, sy] = size(I{i});  % size of images  
   J = I{i}; v = var(J(:)); ave = sqrt(v);  
   for o = 1 : no
       re = filter2(allFilterR{1, o}, I{i}, 'same');
       im = filter2(allFilterI{1, o}, I{i}, 'same');
       energy = re.*re + im.*im; 
   
       
%        resultR{i, o} = re/ave; 
%        resultI{i, o} = im/ave; 
%        resultE{i, o} = energy/v/2.; 
       
       resultR{i, o} = re; 
       resultI{i, o} = im; 
       resultE{i, o} = energy/2; 
       
   end
   %adjustedI{i} = (J-mean(J(:)))/ave; 
end