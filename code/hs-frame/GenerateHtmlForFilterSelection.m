% GenerateHtml - Generates html documentation for experimental results.
%close all;
html_dir = savingFolder;
imgFold = 'images';
if ~exist(html_dir,'dir')
    mkdir(html_dir);
end
html_path = sprintf('%s%s.html',html_dir,imgFold);
html = fopen(html_path,'w');

html_img_dir = [html_dir imgFold filesep];
if ~exist(html_img_dir,'dir')
    mkdir(html_img_dir);
end
numOfTemplate = numScale + numDoGTemplate + 1;

%% html header
% note: modified on Jan 13, 2010. Some url's are now absolute. Only 
%   copied files are linked with relative url.
tmp = '<!DOCTYPE HTML PUBLIC "-//W3C//DTD HTML 4.01 Transitional//EN">\n';
tmp = [tmp '<html>\n'];
tmp = [tmp '<head>\n'];
tmp = [tmp '<title>Image synthesis by multi-scale Gabors and DoGs </title>\n'];
tmp = [tmp '</head>\n'];
fprintf(html,tmp);
fprintf(html, '<body>\n');
fprintf( html, '<div id="header">\n');
fprintf( html, '<h1>Image synthesis by multi-scale Gabors and DoGs</h1></div>\n' );
fprintf( html, '\n<div id="content">\n');


%% explain the parameters
fprintf(html, '<div style="border-top:1 solid #dddddd; margin-top:0.3em;"></div> ');
fprintf(html, '\n<p><b>(I) Parameters</b>: <br>');
fprintf( html, sprintf('The number of training images is %d. ',numImg));
fprintf( html, sprintf('Each image is resized to  %d(height) and %d(width). ',sx,sy));
fprintf( html, sprintf('The number of sketches used to reconstruct each image is fixed to %d.',numSketch));
fprintf( html, sprintf('Two kinds of filters are used: %d scales Gabor and %d scales DoG', numScale, numDoG));
fprintf( html, sprintf('The sizes of the Gabor elements are %d, %d, and %d respectively', 35, 25, 27));
fprintf( html, sprintf('The sizes of DoG elements are are %d, and %d respectively, ', 77, 55));
fprintf( html, sprintf('The DoG elements are illustrated by circles. The radius of a circle is about half of that of the blob represented')); 
fprintf( html, sprintf('by the corresponding DoG elements. Larger circles are darker than smaller ones.'));
fprintf( html, sprintf('Local perturbation algong X or Y axis is fixed to %d, while local pertubation along orientation is %d', locationShiftLimit, orientShiftLimit));
fprintf( html, '\n<br><br>');
fprintf( html, '</p>\n ');


%% show learned templates
fprintf(html, '<div style="border-top:1 solid #dddddd; margin-top:0.3em;"></div> ');
fprintf(html, '\n<p><b>(II) Learned templates (Layout of Gabor and DoG): </b>. ');
fprintf( html, '\n</p>' );
fprintf( html, '\n<p>' );
for cc = 1:numOfTemplate
    new_img_name = sprintf('LearnedTemplate%d.png',cc);
    src = fullfile(output, new_img_name);
    dst = fullfile(html_dir,imgFold,new_img_name);
    copyfile(src,dst);
    fprintf( html, sprintf('\n(%d)', cc));
    fprintf( html, '\n<img src="%s" alt="" width="70"/>', fullfile(imgFold,new_img_name));
end;
fprintf( html, '\n<br><br><br>');
fprintf( html, '\n</p>\n' );


%% training images and their reconstructed images
fprintf(html, '<div style="border-top:1 solid #dddddd; margin-top:0.3em;"></div> ');
fprintf(html, '\n<p><b>(III) Traing images (shown in column 1) , Filters (shown in column 2), their reconstructed images (shown in column 3), and residual images (shown in column 4):</b></p>');
% read the training examples
for (i = 1:numImg)
    new_img_name = sprintf('ReconstructedImage%d.png', i);
    src = fullfile(output,new_img_name);   
    dst = fullfile(html_dir,imgFold,new_img_name);
    copyfile(src,dst);
end

for (i = 1:numImg)
    new_img_name = sprintf('DeformedTemplate%d.png', i);
    src = fullfile(output,new_img_name); 
    dst = fullfile(html_dir,imgFold,new_img_name);
    copyfile(src,dst);
end

% generate corresponding html
fprintf( html, '\n<p>' );
for (i = 1:numImg)
    new_img_name = sprintf('ReconstructedImage%d.png', i);
    new_img_name1 = sprintf('DeformedTemplate%d.png', i);
    fprintf( html, '\n<p>' );
    fprintf( html, '\n<img src="%s" alt="" height=80/>', fullfile(imgFold, new_img_name) );
    fprintf( html, '\n<img src="%s" alt="" height=80/>', fullfile(imgFold, new_img_name1) );
    fprintf( html, '\n</p>' );
end
fprintf( html, '\n</p>\n' );
fprintf( html, '\n<br><br><br><br><br><p>' );
fprintf( html, '\n</p>' );


%% finishing off
fprintf( html, '\n\n\n</div> \n');
fprintf( html, '<div id="last" class="footer"></div>' );
fprintf(html, '</body> </html> \n');
fclose(html);
disp('finished generating Html... go to document folder and click the html');
close all

if exist(output,'dir')
   rmdir(output,'s');
end
