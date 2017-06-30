function sym = displayMatchedTemplate(latticeSize, selectedRow,selectedCol,selectedO, selectedS,selectedMean,allsymbol,nGaborOri)
% displayTemplate - Display the matched symbolic template of an active basis template. 

%% display Gabor symbols one by one
nGaborScale = numel(allsymbol) / nGaborOri;
sym = zeros(latticeSize);
for k = 1:length(selectedRow)
    scale = selectedS(k) + 1;
    ori = selectedO(k) + 1;
    col = selectedCol(k);
    row = selectedRow(k);
    if scale < 1 || scale > nGaborScale
        continue;
    end
    sym = drawGaborSymbol( sym, allsymbol, row, col, ori, nGaborOri, scale, sqrt(selectedMean(k)) );
end
sym = uint8(255 * (sym-min(sym(:)))/(max(sym(:))-min(sym(:))));


