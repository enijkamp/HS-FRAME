function feature = featureSPM(maps,nLevel,t)
% compute the max feature response in cells of different sizes
factor = 2;

nMap = numel(maps);

%% compute the features in finest scale
nRow = factor^(nLevel-1);
nCol = nRow;

pyramid = cell(nLevel,1);

hist = zeros(nMap,nRow,nCol);
for iMap = 1:nMap
    [sx sy]=size(maps{iMap});
    x_step_width = round(sx/nRow);
    y_step_width = round(sy/nCol);
    for iRow = 1:nRow
        for iCol = 1:nCol
            subMap = maps{iMap}((iRow-1)*x_step_width+1:min(sx,iRow*x_step_width),(iCol-1)*y_step_width+1:min(sy,iCol*y_step_width));
            maxVal = max(0,max(subMap(:))-t); % find maximum in the cell, subject to thresholding
            hist(iMap,iRow,iCol)=maxVal;
        end
    end
end
pyramid{nLevel}=hist;


%% combine finest scale features to more coarse scales
for iLevel = nLevel-1:-1:1
    jLevel = iLevel +1;
    nRow = size(pyramid{jLevel},2);
    nCol = size(pyramid{jLevel},3);

    nRow = nRow/factor;
    nCol = nCol/factor;
    hist = zeros(nMap,nRow,nCol);
    for iRow = 1:nRow
        for iCol = 1:nCol
            maxVal = zeros(nMap,1)-Inf;

            for jRow = (iRow-1)*factor+1:iRow*factor
                for jCol = (iCol-1)*factor+1:iCol*factor
                    maxVal = max(maxVal,pyramid{jLevel}(:,jRow,jCol));
                end
            end
            pyramid{iLevel}(:,iRow,iCol)=maxVal;
        end
    end
end

%% lump the featurescells into a vector
feature=[];
for iLevel = nLevel:-1:1
    feature = [feature pyramid{iLevel}(:)'];
end
