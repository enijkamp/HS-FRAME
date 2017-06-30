function [allFiltered] = applyfilterfft(I, allFilter, localHalfx, localHalfy, thresholdFactor);
% filter image by a bank of filters
% I: input images
% allFilter: filter bank
numImage = size(I, 2);    % number of images
numOrient = size(allFilter, 2);   % number of orientations
h = (size(allFilter{1}, 1)-1)/2;  % half size of filters
allFiltered = cell(numImage, numOrient);  % filtered images
for (i = 1 : numImage)
   [sx, sy] = size(I{i});  % size of images
   tot = 0.; 
   fftI = fft2(I{i}, sx+h+h, sy+h+h);
   for (o = 1 : numOrient)
      out = ifft2(fftI.*fft2(allFilter{1, o}, sx+h+h, sy+h+h)); 
      filtered = out(h+1:h+sx, h+1:h+sy); 
      re = real(filtered); im = imag(filtered); 
      energy = re.*re + im.*im;  % compute the local energy
      energy([1:h sx-h:sx],:) = 0.; energy(:,[1:h sy-h:sy]) = 0.; 
      allFiltered{i, o} = energy; 
      tot = tot + sum(sum(energy((h+1):(sx-h-1), (h+1):(sy-h-1))))/(sx-2*h-1)/(sy-2*h-1); 
   end
   ave = tot/numOrient; 
   for o = 1 : numOrient
       allFiltered{i, o} = allFiltered{i, o}/ave;  % normalizing by whitening
   end
   ClocalNormalize(sx, sy, numOrient, h, localHalfx, localHalfy, allFiltered(i, :), thresholdFactor); 
end
    
 
