%% show experimental results
output = fullfile(savingFolder, 'output');
if ~exist(output,'dir')
   mkdir(output);
else
   rmdir(output,'s');
   mkdir(output);
end
% figure; 
% col = 6;%each row has 6 images 
% row = 1 + floor((numImg+1.)/2); 
% 
% for (i = 1:numScale+numDoGTemplate) %show templates, the last one is for DoG
%    subplot(row, col, i); imshow(-sym{i}, []);  % sym include all scales of templates and one DoG template
% end
% subplot(row, col, numScale+numDoGTemplate+1); imshow(-syma, []);  % syma is one comprehensive template 
% showIndex = 6; 
% for (i = 1 : numImg)
%     showIndex = showIndex + 1;  subplot(row, col,  showIndex); 
%     imshow(I{i}, []);       % all images
%     showIndex = showIndex + 1;  subplot(row, col,  showIndex); 
%     imshow(-Asyma{i}, []);  % sketch templates for all images
%     showIndex = showIndex + 1;  subplot(row, col,  showIndex); 
%     imshow(Asym{i}, []);    % reconstructed images
% end

%% save experimental results
for (i = 1:numScale+numDoGTemplate) %save learned templates
   towrite =  -sym{i}; 
   if max(towrite(:))-min(towrite(:))~=0
       towrite = (towrite-min(towrite(:)))/(max(towrite(:))-min(towrite(:)));
   end;
   imwrite(towrite, [output filesep sprintf('LearnedTemplate%d.png',i)],'PNG'); 
end

towrite =  -syma;
if max(towrite(:))-min(towrite(:))~=0
       towrite = (towrite-min(towrite(:)))/(max(towrite(:))-min(towrite(:)));
end;
imwrite(towrite, [output filesep sprintf('LearnedTemplate%d.png',numScale+numDoGTemplate+1)],'PNG'); 

towrite =  -syma_top;
if max(towrite(:))-min(towrite(:))~=0
       towrite = (towrite-min(towrite(:)))/(max(towrite(:))-min(towrite(:)));
end;
imwrite(towrite, [output filesep sprintf('LearnedTemplate_top.png')],'PNG');    


space = 4;
for (i = 1 : numImg)%save images, filters and recovered images
    h = size(I{i},1);
    w = size(I{i},2);
    towrite = ones(h, 4*w+3*space); % the ones make the background white, zeros make it black
    towrite(1:h,1:w) = (I{i}-min(I{i}(:)))/(max(I{i}(:))-min(I{i}(:)));%original image
    towrite(1:h,w+space+1:2*w+space) = (-Asyma{i}-min(-Asyma{i}(:)))/(max(-Asyma{i}(:))-min(-Asyma{i}(:)));%filters
    towrite(1:h,2*w+2*space+1:3*w+2*space) = (Asym{i}-min(Asym{i}(:)))/(max(Asym{i}(:))-min(Asym{i}(:)));%filters   
    residual = I{i} - Asym{i};
    towrite(1:h,3*w+3*space+1:4*w+3*space) = (residual-min(residual(:)))/(max(residual(:))-min(residual(:)));%filters
    imwrite(towrite, [output filesep sprintf('ReconstructedImage%d.png',i)],'PNG');
end

for (i=1:numImg)
  towrite =  -Asyma_top{i};
  if max(towrite(:))-min(towrite(:))~=0
       towrite = (towrite-min(towrite(:)))/(max(towrite(:))-min(towrite(:)));
  end;
  imwrite(towrite, [output filesep sprintf('topFeatureTemplate%d.png',i)],'PNG');   
end


space = 4;
for (i = 1 : numImg)%save images, filters and recovered images
    
    h = size(I{i},1);
    w = size(I{i},2);
    towrite = ones(h, (numScale+numDoGTemplate)*w+(numScale+numDoGTemplate-1)*space); % the ones make the background white, zeros make it black
    
    for(j=1:numScale+numDoGTemplate)
        
       p=-Asym0{(i-1)*(numScale+numDoGTemplate)+j};
       towrite(1:h,(j-1)*(w+space)+1:(j-1)*(w+space)+w) = (p-min(p(:)))/(max(p(:))-min(p(:)));    
    
    end
      
    imwrite(towrite, [output filesep sprintf('DeformedTemplate%d.png',i)],'PNG');
end


