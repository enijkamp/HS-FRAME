function [S2T, S3T] = hierachicalTemplate(numPart, part_sx, part_sy, sx, sy, rotateShiftLimit, nOrient, numRotate, template, nScaleGabor, partRotationRange,PartLocX, PartLocY)



%%%%%%%% handeling the large template
S3SelectedRow = zeros(1,numPart);
S3SelectedCol = zeros(1,numPart);
S3SelectedOri = zeros(1,numPart);

for iPart = 1:numPart
    S3SelectedRow(iPart) = PartLocX(iPart) - 1 + floor(part_sx/2);  % index starts from 0?????
    S3SelectedCol(iPart) = PartLocY(iPart) - 1 + floor(part_sy/2);
end

for i = 1 : numPart   % move the center of the template to the origin 
	S3SelectedRow(i) = S3SelectedRow(i) - floor(sx/2);
	S3SelectedCol(i) = S3SelectedCol(i) - floor(sy/2);
end

[allS3SelectedRow, allS3SelectedCol, allS3SelectedOri]=rotateS3Template(S3SelectedRow,S3SelectedCol,S3SelectedOri,rotateShiftLimit,nOrient,numRotate,numPart);

%%%%%%%%% split the object template into non-overlapping partial templates

selectedS2Filter = cell(numPart, 1);%orientation and location of selected Gabors
selectedS2Row = cell(numPart, 1);
selectedS2Col = cell(numPart, 1);
selectedS2lambda = cell(numPart, 1);%weighting parameter for scoring template matching
 
PartOnOff = ones(numPart,1); % all parts are selected initially
selectedPart = find(PartOnOff);


for iPart = 1:numPart
   ind = find( template.selectedRow >= PartLocX(iPart) & template.selectedRow < PartLocX(iPart) + part_sx & ...
          template.selectedCol >= PartLocY(iPart) & template.selectedCol < PartLocY(iPart) + part_sy );
   selectedS2Filter{iPart} = template.selectedFilter(ind);
   selectedS2Row{iPart} = floor( template.selectedRow(ind) - PartLocX(iPart) +1 );  % index starts from 1
   selectedS2Col{iPart} = floor( template.selectedCol(ind) - PartLocY(iPart) +1 );  % index starts from 1
   selectedS2lambda{iPart} = template.selectedLambdas(ind);  
   
end

for iPart = 1:numPart     % move the center of the part template to the origin 
    selectedS2Row{iPart} = selectedS2Row{iPart}-floor(part_sx/2);  % index starts from 1
    selectedS2Col{iPart} = selectedS2Col{iPart}-floor(part_sy/2);  % index starts from 1
end

[allS2SelectedRow, allS2SelectedCol, allS2SelectedFilter]=rotateS2Template(selectedS2Row,selectedS2Col,selectedS2Filter, nScaleGabor, nOrient, partRotationRange, numPart);


numPartRotate = length(partRotationRange);

% organize the S2 part template
S2T = cell( numPartRotate, numPart );
for iPart = 1:numPart
   for r = 1:numPartRotate
      S2T{r,iPart} = struct( 'selectedRow',single( round(allS2SelectedRow{iPart,r}(:))),...
            'selectedCol', single( round(allS2SelectedCol{iPart,r}(:))),...
            'selectedFilter', single(allS2SelectedFilter{iPart,r}(:)),...           
            'selectedLambdas', single(selectedS2lambda{iPart}(:)),...
            'numSelected', single(length(allS2SelectedRow{iPart,r})) );
   end
end

% organize the S3 large template
S3T = cell(numRotate,1);
for r = 1:numRotate % this is the rotation of the S3 template
     %rot = rotationRange(r);
     selectedTransform = zeros(length(selectedPart),1,'single');
     for j = 1:length(selectedPart)
           selectedTransform(j) = find( allS3SelectedOri(r,j) == partRotationRange );  % change the actual orientation into the index of transformation
     end
       S3T{r} = struct( 'selectedRow',single( round(allS3SelectedRow(r,selectedPart))),...
             'selectedCol', single( round(allS3SelectedCol(r,selectedPart))),...
             'selectedInd', single(selectedPart) - 1,...        % index starts from 0;
             'selectedTransform', selectedTransform - 1,...
             'selectedLambdas', ones(length(selectedPart),1,'single'),...
             'selectedLogZ', single( 0*ones(length(selectedPart),1) ) );
end

