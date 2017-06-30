denseX = -floor(partSizeX/2) + (1:partSizeX);
denseY = -floor(partSizeY/2) + (1:partSizeY);
count = 0;
inRow = zeros(length(denseX)*length(denseY),1,'single');
inCol = zeros(length(denseX)*length(denseY),1,'single');
for y = denseY
    for x = denseX
        count = count+1;
        inRow(count) = x;
        inCol(count) = y;
    end
end

% non-overlapping partial templates
for iPart = 1:numCandPart
	centerx = floor(partSizeX/2);
    centery = floor(partSizeY/2);
    clusters{c}.selectedx{iPart} = clusters{c}.selectedx{iPart}-centerx;
    clusters{c}.selectedy{iPart} = clusters{c}.selectedy{iPart}-centery;
	for r = 1:length(partRotationRange)
		rot = partRotationRange(r);
        tScale = 0; rScale = 1; cScale = 1;
        inS = zeros(numel(inRow),1,'single');
        [clusters{c}.allSelectedx{iPart,r} clusters{c}.allSelectedy{iPart,r} tmpO]...
            = mexc_TemplateAffineTransform(tScale,rScale,cScale,rot,...
            single(clusters{c}.selectedx{iPart}),...
            single(clusters{c}.selectedy{iPart}),single(clusters{c}.selectedOrient{iPart}),inS,numOrient);
        tmpO(tmpO<0) = tmpO(tmpO<0) + numOrient;
        tmpO(tmpO>=numOrient) = tmpO(tmpO>=numOrient) - numOrient;
        clusters{c}.allSelectedOrient{iPart,r} = tmpO;
	end
end

% overlapping partial templates
for iPart = 1:numCandPart
	centerx = floor(partSizeX/2 + partMarginX);
    centery = floor(partSizeY/2 + partMarginY);
    clusters{c}.largerSelectedx{iPart} = clusters{c}.largerSelectedx{iPart}-centerx;
    clusters{c}.largerSelectedy{iPart} = clusters{c}.largerSelectedy{iPart}-centery;
	for r = 1:length(partRotationRange)
		rot = partRotationRange(r);
        tScale = 0; rScale = 1; cScale = 1;
        inS = zeros(numel(inRow),1,'single');
        [clusters{c}.largerAllSelectedx{iPart,r} clusters{c}.largerAllSelectedy{iPart,r} tmpO]...
            = mexc_TemplateAffineTransform(tScale,rScale,cScale,rot,...
            single(clusters{c}.largerSelectedx{iPart}),...
            single(clusters{c}.largerSelectedy{iPart}),single(clusters{c}.largerSelectedOrient{iPart}),inS,numOrient);
        tmpO(tmpO<0) = tmpO(tmpO<0) + numOrient;
        tmpO(tmpO>=numOrient) = tmpO(tmpO>=numOrient) - numOrient;
        clusters{c}.largerAllSelectedOrient{iPart,r} = tmpO;
	end
end



