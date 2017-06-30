function GenerateHtml(numIteration)
% GenerateHtml - Generates html documentation for the learning results.
%

 end
% % whether show bounding boxes for selected partial templates
% boundSelectedParts = true;
% 
% % load the parameters
% load Config
% 
% %numIteration = 25;
% 
% % delete the previous version
% html_dir = 'document/html';
% if ~exist(html_dir,'dir')
%     mkdir(html_dir);
% end
% imFolder = sprintf('HABMorph_%s_files',category);
% html_path = sprintf('%s/%s.html',html_dir,imFolder(1:end-6));
% html = fopen(html_path,'w');
% 
% html_img_dir = sprintf('document/html/%s/',[imFolder]);
% if ~exist(html_img_dir,'dir')
%     mkdir(html_img_dir);
% end
% 
% %% html header
% tmp = '<!DOCTYPE HTML PUBLIC "-//W3C//DTD HTML 4.01 Transitional//EN">\n';
% tmp = [tmp '<html>\n'];
% tmp = [tmp '<head>\n'];
% tmp = [tmp '<title>Learning hierarchical active basis template</title>\n'];
% tmp = [tmp '<link rel="stylesheet" href="http://www.stat.ucla.edu/~zzsi/plain_and_simple.css" type="text/css" media="screen" />\n'];
% tmp = [tmp '<script type="text/javascript" src="http://www.stat.ucla.edu/~zzsi/1.js"></script>\n'];
% tmp = [tmp '</head>\n'];
% fprintf(html,tmp);
% fprintf(html, '<body>\n');
% fprintf( html, '<div id="header">\n');
% fprintf( html, '<h1>Weakly supervised learning of hierarhical active basis</h1></div>\n' );
% fprintf( html, '\n<div id="content">\n');
% 
% %% link to project page
% fprintf( html, '\n<p><a href="http://www.stat.ucla.edu/~zzsi/hab.html">Project home</a></p>\n' );
% 
% %% table of content
% fprintf( html, '<div id="content">\n' );
% fprintf( html, '<div id="TableOfContents">\n' );
% fprintf( html, '<p>Contents</p>\n' );
% fprintf( html, '<ul>\n' );
% fprintf( html, '<li>\n' );
% fprintf( html, '<A href="#download">Download</A>\n' );
% fprintf( html, '</li>\n' );
% fprintf( html, '<li>\n');
% fprintf( html, '<a href="#traindata">Training examples</a>\n' );
% fprintf( html, '</li>\n' );
% fprintf( html, '<li>\n' );
% fprintf( html, '<a href="#parts">Learned partial templates</a>\n' );
% fprintf( html, '</li>\n' );
% fprintf( html, '<li>\n' );
% fprintf( html, '<a href="#matched">Shape matching by hierarchical active basis</a>\n' );
% fprintf( html, '</li>\n' );
% fprintf( html, '</ul>\n' );
% fprintf( html, '</div>\n' );
% 
% %% for download
% fprintf( html, '<div style="border-bottom:1 solid #dddddd; margin-top:0.3em;"></div>\n' );
% fprintf( html, '<a name="download"></a> <table cellspacing="10" cellpadding="10" class="center" width="60%%">\n' );
% fprintf( html, '\n<tr><td>\n' );
% fprintf( html, sprintf('\n<b>Code and data: (<a href="%s">ZIP</a>).</b>\n',zipfilename) );
% fprintf( html, '\n</td>\n' );
% fprintf( html, '\n<td>\n' );
% fprintf( html, '\n<a href="http://www.stat.ucla.edu/~zzsi/HAB/hab_changelog.html">Change Log</a>\n' );
% fprintf( html, '\n</td>\n' );
% fprintf( html, '\n</tr>\n' );
% fprintf( html, '\n<tr>\n' );
% fprintf( html, '\n<td colspan=2 align=left>\n' );
% fprintf( html, '\nRun StartFromHere.m in matlab to learn the hierarchical active basis template. Run GenerateHtml.m to generate this webpage.\n' );
% fprintf( html, '\n</td>\n' );
% fprintf( html, '\n</tr>\n' );
% fprintf( html, '\n</table>\n' );
% 
% %% training examples
% fprintf(html, '<div style="border-top:1 solid #dddddd; margin-top:0.3em;"></div><h2> ');
% fprintf(html, '<a name="traindata"></a>Training examples');
% fprintf(html, ' </h2>');
% fprintf( html, ['\n<p>We learn a 2-level hierarchical active basis template from images with articulated objects with unknown scale, rotation and translation. The object '...
%     sprintf('template size is %d * %d pixels; the part template size is %d * %d pixels.', templateSize(1), templateSize(2), partSizeX, partSizeY)...
%     ' Each part may translate and rotate locally relative to the object center. And each edge element can translate and rotate locally relative to the part center. We use EM-type algorithm to obtain the part and object templates as well as the hierarchical deformation on each image.'...
%     sprintf(' Number of EM iterations = %d. ',numIteration) ...
%     sprintf(' Each part is an active basis model of Gabor elements with local shift <= %d pixels, and orientation shift <= %d * PI/%d. ', locationShiftLimit, orientShiftLimit, numOrient) ...
%     sprintf(' In total, there are %d Gabor elements in the object template. ', numElement ) ...
%     ' </p>\n']);
% % read the training examples
% Iname = dir(['positiveImage/' category '/*.jpg']);
% % move the images to img/ folder
% for i = 1:length(Iname)
%     new_img_name = Iname(i).name;
%     src = sprintf(['positiveImage/' category '/%s'],new_img_name);
%     dst = sprintf('%s/%s',html_img_dir, new_img_name);
%     copyfile(src,dst);
% end
% % generate corresponding html
% fprintf( html, '\n<p>' );
% for i = 1:length(Iname)
%     fprintf( html, '<img src="%s" alt="" height=80/>', sprintf('%s/%s',imFolder,Iname(i).name) );
% end
% fprintf( html, '\n</p>\n' );
% % show the starting image
% fprintf( html, sprintf('\n<h3> Learning from example %d</h3>\n', starting ) );
% fprintf( html, '\n<p> The template is initialized from the following example image. </p>' );
% fprintf( html, '\n<p>' );
% fprintf( html, '<img src="%s" alt="" height = 100 />', sprintf('%s/%s',imFolder,Iname(starting).name) );
% fprintf( html, '\n</p> <p></p>\n' );
% 
% 
% %% show the whole template
% fprintf( html, '<div style="border-bottom:1 solid #dddddd; margin-top:0.3em;"></div>\n <h2>Learned object template</h2> <a name="parts"></a>' );
% new_img_name = sprintf( 'template_iter%d.png', numIteration-1 );
% src = sprintf('output/%s',new_img_name);
% dst = sprintf('%s%s', html_img_dir,new_img_name);
% copyfile(src,dst);
% fprintf( html, sprintf('\n<img src="%s" alt="" height=%d/>\n<br />',sprintf('%s/%s',imFolder,new_img_name),sizeTemplatex) );
% 
% %{
% %% show learned candidate part templates, together with information gain
% fprintf( html, ['<div style="border-bottom:1 solid #dddddd; margin-top:0.3em;"></div>\n <h2>Learned partial templates</h2> <a name="parts"></a><p>Learned final active basis templates' ...
%     ' at non-overlapping subwindows. The number to the lower-left of each template' ...
%     ' is the average MAX2 score of the plate (partial template) evaluated on image examples.' ...
%     ' The color of the partial template indicates its information gain as the "hotness" measure. '] );
% if boundSelectedParts
%     fprintf( html, ...
%         ['The Hierarhical Active Basis template is composed of selected partial templates with high information gains.' ...
%         ' They are indicated by bounding boxes. ' ]);
% end
% fprintf( html, '</p>\n\n<p>' );
% 
% % find the range of average MAX2 score
% load(sprintf('working/aveMAX2_iter%d.mat',numIteration));
% MAX2 = aveMAX2;
% max_MAX2 = max(MAX2(:));
% lowerbound_proportion = .6;
% lowerbound_MAX2 = max_MAX2 * lowerbound_proportion;
% 
% % dipslay the color legend indicating "hotness" of partial templates
% im = zeros(20,500,3);
% for j = 1:size(im,2)
% 	s = j / size(im,2);
% 	[r g b] = hotness(s);
% 	im(:,j,1) = r; im(:,j,2) = g; im(:,j,3) = b;
% end
% legendFileName = 'ColorLegend.png';
% imwrite(uint8(im),sprintf('%s/%s',html_img_dir,legendFileName));
% fprintf( html, '\n <p> low <img src="%s" alt="color legend" width = 500/> high </p>', ...
%     sprintf('%s/%s',imFolder,legendFileName) );
% 
% 
% count = 0;
% 
% for ix = 1:length(PartLocX0)
%     for iy = 1:length(PartLocY0)
%         count = count + 1;
%     	aveMAX2score = MAX2(iy+(ix-1)*length(PartLocY0));
%         fprintf( html, '%.1f', aveMAX2score );
%         % move the result image to /img
%         new_img_name = sprintf('template_iter%d_part%d_%d.png',numIteration, PartLocX0(ix),PartLocY0(iy));
%         src = sprintf('output/%s',new_img_name);
%         dst = sprintf('%s/%s', html_img_dir,new_img_name);
%         im = imread(src);
%         mask = im(:,:,1) < 255;
% 		if size(im,3) == 1
% 			im = repmat(im,[1,1,3]);
% 		end
% 		im = double(im);
% 		
% 		% apply a custom color to it
% 		ind = ( max(aveMAX2score,lowerbound_MAX2) - lowerbound_MAX2 ) / (max_MAX2 - lowerbound_MAX2);
% 		tmp = im(:,:,1);
% 		[r g b] = hotness(ind);
% 		tmp(mask) = r;
% 		im(:,:,1) = tmp; % red
% 		tmp(mask) = g;
% 		im(:,:,2) = tmp; % green
% 		tmp(mask) = b;
% 		im(:,:,3) = tmp; % blue
%         
%         % add a bounding box if it is selected
%         if boundSelectedParts && PartOnOff(count) > 0
%             im([1:3,end-2:end],:,:) = 100; im(:,[1:3,end-2:end],:) = 100;
%         end
% 		
% 		imwrite(uint8(im),dst);
%         fprintf( html, '\n<img src="%s" alt="" height=80/>', sprintf('%s/%s',imFolder,new_img_name) );
%     end
%     fprintf( html, '\n<br />\n\n' );
% end
% fprintf( html, '\n</p>\n' );
% %}
% 
% 
% %% matched HAB template (SUM3 templates)
% fprintf( html, ['<div style="border-top:1 solid #dddddd; margin-top:0.3em;"></div>\n<a name="matched"></a>' ...
%     '\n <h2>Shape matching by Hierarchical Active Basis</h2> \n <p> A hierarchical active basis (HAB) template is ' ...
%     'composed of a small number of partial templates which can translate and rotate locally relative to other partial templates.' ... 
%     'The following pairs of images show the matched HAB templates' ...
%     ' on example images. </p>\n']);
% fprintf( html, '\n<p>' );
% % move the images to img/ folder
% for i = 1:length(Iname)
%     % matched template
%     new_img_name = sprintf('matchedS3T%d_iter%d.png',i,numIteration);
%     src = sprintf('output/%s',new_img_name);
%     dst = sprintf('%s/%s',html_img_dir, new_img_name);
%     % copyfile(src,dst);
%     im = imread(src);
%     if size(im,3) == 1
%         im = repmat(im,[1,1,3]);
%     end
%     im([1:3,end-2:end],:,:) = 100; im(:,[1:3,end-2:end],:) = 100;
%     imwrite(im,dst);
% 
%     % matched image patch
%     new_img_name = sprintf('matchedS3image%d.png',i);
%     src = sprintf('output/%s',new_img_name);
%     dst = sprintf('%s/%s',html_img_dir, new_img_name);
%     % copyfile(src,dst);
%     im = imread(src);
%     if size(im,3) == 1
%         im = repmat(im,[1,1,3]);
%     end
%     im([1:3,end-2:end],:,:) = 100; im(:,[1:3,end-2:end],:) = 100;
%     imwrite(im,dst);
% end
% % generate corresponding html
% fprintf( html, '\n<p>' );
% numCol = 4;
% 
% for i = 1:length(Iname)
%     fprintf( html, '<img src="%s" alt="" height=120/>', sprintf('%s/%s',imFolder,sprintf('matchedS3image%d.png',i)) );
%     fprintf( html, '<img src="%s" alt="" height=120/>\n', sprintf('%s/%s',imFolder,sprintf('matchedS3T%d_iter%d.png',i,numIteration)) );
%     if rem(i,numCol) == 0
%     	fprintf( html, '<br />\n' );
%     end
% end
% fprintf( html, '\n</p>\n' );
% 
% 
% %% finishing off
% fprintf( html, '\n\n\n</div> \n');
% fprintf( html, '<div id="last" class="footer"></div>' );
% fprintf(html, '</body> </html> \n');
% fclose(html);
% disp('finished generating Html.. check here:');
% disp([pwd '/' html_path]);
% 
% %
% %% Tool: color map
% %
% function [r g b] = hotness(s)
% r = 200 + 50 * s;
% g = 100 + 150 * (1-s);
% b = 100 + 150 * (1-s);
