function sym = displayTemplate(selectedRow,selectedCol,selectedO, selectedS,selectedMean,allsymbol,nGaborOri,subsample)
% displayTemplate - Display the symbolic template of an active basis template. 

%% determine the rectangular size of the template
% Note: The row/column indices are relative to the origin (0,0).
maxScale = max(selectedS); % starts from 0
margin = ceil( size(allsymbol{maxScale*nGaborOri+1},1) / 2 ); % to allow enough space for Gabor symbol
selectedRow = selectedRow * subsample;
selectedCol = selectedCol * subsample;
selectedS = selectedS + 1; % now starts from 1
top = min(selectedRow);
down = max(selectedRow);
left = min(selectedCol);
right = max(selectedCol);
height = down - top + 2 * margin;
width = right - left + 2 * margin;
% make it square:
if height > width
    width = height;
else
    height = width;
end

%% display Gabor symbols one by one
nGaborScale = numel(allsymbol) / nGaborOri;
sym = zeros(height,width);
for k = 1:length(selectedRow)
    scale = selectedS(k);
    ori = selectedO(k) + 1;
    col = selectedCol(k) + floor(width/2);
    row = selectedRow(k) + floor(height/2);
    if scale < 1 || scale > nGaborScale
        continue;
    end
    sym = drawGaborSymbol( sym, allsymbol, row, col, ori, nGaborOri, scale, sqrt(selectedMean(k)) );
end
sym = uint8(255 * (sym-min(sym(:)))/(max(sym(:))-min(sym(:))));


