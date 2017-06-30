
%% Prepare output variables for learning
selectedOrient = cell(numCandPart, 1);  % orientation and location of selected Gabors
selectedx = cell(numCandPart, 1);
selectedy = cell(numCandPart, 1);
selectedlambda = cell(numCandPart, 1); % weighting parameter for scoring template matching
selectedLogZ = cell(numCandPart, 1); % normalizing constant

largerSelectedOrient = cell(numCandPart, 1);  % orientation and location of selected Gabors
largerSelectedx = cell(numCandPart, 1);
largerSelectedy = cell(numCandPart, 1);
largerSelectedlambda = cell(numCandPart, 1); % weighting parameter for scoring template matching
largerSelectedLogZ = cell(numCandPart, 1); % normalizing constant

commonTemplate = cell( numCandPart, 1 );
for i = 1:length(commonTemplate)
	commonTemplate{i} = single(zeros(partSizeX, partSizeY)); % template of active basis 
end

deformedTemplate = cell(1, numImage); % templates for training images 
for img = 1 : numImage
    deformedTemplate{img} = single(zeros(partSizeX, partSizeY));  
end

allSelectedx = cell(numCandPart, numRotate); 
allSelectedy = cell(numCandPart, numRotate);
allSelectedOrient = cell(numCandPart, numRotate);
largerAllSelectedx = cell(numCandPart, numRotate);
largerAllSelectedy = cell(numCandPart, numRotate);
largerAllSelectedOrient = cell(numCandPart, numRotate);

%% Initialize learning from the starting image
numImage0 = 1; locationShiftLimit0 = 0; orientShiftLimit0 = 0; 
SUM1mapLearn0 = cell(1, numOrient);
load( sprintf('working/ImageAndFeature_%d.mat',starting) );
for orient = 1 : numOrient
    SUM1mapLearn0{1, orient} = single(zeros(sizeTemplatex, sizeTemplatey));
    sizex = size(SUM1mapFind{originalResolution,orient},1);
    sizey = size(SUM1mapFind{originalResolution,orient},2);
    Ccopy(SUM1mapLearn0{1, orient}, SUM1mapFind{originalResolution,orient}, startx-1, starty-1, 0, 0, sizeTemplatex, sizeTemplatey, sizex, sizey, 0); 
end

deformedTemplate0{1} = single(zeros(sizeTemplatex, sizeTemplatey));

disp('start from single image learning'); 
tic

tmpSelectedx = zeros(numElement,1);
tmpSelectedy = zeros(numElement,1);
tmpSelectedo = zeros(numElement,1);
tmpSelectedlambda = zeros(numElement,1);
tmpSelectedlogz = zeros(numElement,1);
tmpCommonTemplate = zeros(sizeTemplatex,sizeTemplatey,'single');
CsharedSketch(numOrient, locationShiftLimit0, orientShiftLimit0, subsample, ... % about active basis  
	   numElement, numImage0, sizeTemplatex, sizeTemplatey, SUM1mapLearn0, ... % about training images 
	   halfFilterSize, Correlation, allSymbol(1, :), ... % about filters
	   numStoredPoint, storedlambda, storedExpectation, storedLogZ, ... % about exponential model 
	   tmpSelectedo, tmpSelectedx, tmpSelectedy, tmpSelectedlambda, tmpSelectedlogz, ... % learned parameters
	   tmpCommonTemplate, deformedTemplate0); % learned templates


% split the object template into non-overlapping partial templates
for iPart = 1:numCandPart
	ind = find( tmpSelectedx >= PartLocX(iPart) & tmpSelectedx < PartLocX(iPart) + partSizeX & ...
			tmpSelectedy >= PartLocY(iPart) & tmpSelectedy < PartLocY(iPart) + partSizeY );
    selectedOrient{iPart} = single( tmpSelectedo(ind) );
    selectedx{iPart} = single( floor( tmpSelectedx(ind) - PartLocX(iPart) ) );
    selectedy{iPart} = single( floor( tmpSelectedy(ind) - PartLocY(iPart) ) );
    selectedlambda{iPart} = single( tmpSelectedlambda(ind) );
    selectedLogZ{iPart} = single( tmpSelectedlogz(ind) );
    filename = sprintf('output/template_iter%d_part%d_%d.png',0,PartLocX(iPart),PartLocY(iPart));
    % add a small margin to make sure all Gabor elements are displayed fully
    im = displayMatchedTemplate([partSizeX+2*halfFilterSize partSizeY+2*halfFilterSize],selectedx{iPart}+halfFilterSize,...
            selectedy{iPart}+halfFilterSize,selectedOrient{iPart},zeros(length(ind),1,'single'),selectedlambda{iPart},allSymbol,numOrient);
    towrite = -double(im);
    towrite = uint8(255 * (towrite-min(towrite(:)))/(1e-10+max(towrite(:))-min(towrite(:))));
    if max(towrite(:)) == 0
        towrite(:) = 255;
    end
    imwrite(towrite,filename);
end

% split the object template into larger overlapping partial templates (context sensitive)
for iPart = 1:numCandPart
	ind = find( tmpSelectedx >= PartLocX(iPart)-partMarginX & tmpSelectedx < PartLocX(iPart) + partSizeX + partMarginX & ...
			tmpSelectedy >= PartLocY(iPart) - partMarginX & tmpSelectedy < PartLocY(iPart) + partSizeY + partMarginY );
    largerSelectedOrient{iPart} = single( tmpSelectedo(ind) );
    largerSelectedx{iPart} = single( tmpSelectedx(ind) - PartLocX(iPart) + partMarginX );
    largerSelectedy{iPart} = single( tmpSelectedy(ind) - PartLocY(iPart) + partMarginY );
    largerSelectedlambda{iPart} = single( tmpSelectedlambda(ind) );
    largerSelectedLogZ{iPart} = single( tmpSelectedlogz(ind) );
    filename = sprintf('output/largertemplate_iter%d_part%d_%d.png',0,PartLocX(iPart),PartLocY(iPart));
    im = displayMatchedTemplate([partSizeX+2*partMarginX+2*halfFilterSize partSizeY+2*partMarginY+halfFilterSize*2],largerSelectedx{iPart}+halfFilterSize,...
            largerSelectedy{iPart}+halfFilterSize,largerSelectedOrient{iPart},zeros(length(ind),1,'single'),largerSelectedlambda{iPart},allSymbol,numOrient);
    towrite = -double(im);
    towrite = uint8(255 * (towrite-min(towrite(:)))/(1e-10+max(towrite(:))-min(towrite(:))));
    if max(towrite(:)) == 0
        towrite(:) = 255;
    end
    imwrite(towrite,filename);
end

RotateTemplate;

disp(['mex-C learning time: ' num2str(toc) ' seconds']);


PartOnOff = ones(numCandPart,1); % all parts are selected initially
S3SelectedRow = zeros(1,numCandPart,'single');
S3SelectedCol = zeros(1,numCandPart,'single');
S3SelectedOri = zeros(1,numCandPart,'single');
for iPart = 1:numCandPart
	S3SelectedRow(iPart) = PartLocX(iPart) - 1 + floor(partSizeX/2);
	S3SelectedCol(iPart) = PartLocY(iPart) - 1 + floor(partSizeY/2);
end
allS3SelectedRow = zeros(numRotate,numCandPart,'single');
allS3SelectedCol = zeros(numRotate,numCandPart,'single');
allS3SelectedOri = zeros(numRotate,numCandPart,'single');
RotateS3Template;

