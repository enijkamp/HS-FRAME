

%% Generate multi-resolution images and compute SUM1, MAX1 maps
disp('start filtering training images at all resolutions'); tic
for img = 1:numImage
    ImageMultiResolution = cell(1,numResolution);
    for j=1:numResolution
        resolution = .6+(j-1)*.2; % so that .8+(originalResolution-1)*.2 = 1
        ImageMultiResolution{j} = imresize(I{img}, resolution, 'nearest');  % images at multiple resolutions
    end
    save(sprintf('working/ImageAndFeature_%d.mat',img),'ImageMultiResolution');

    SUM1mapFind = ApplyFilterfft(ImageMultiResolution, allFilter,...
        localHalfx, localHalfy, thresholdFactor); % filtering images at multiple resolutions
    mexc_Sigmoid(saturation, SUM1mapFind);
    save(sprintf('working/ImageAndFeature_%d.mat',img),'SUM1mapFind','-append');

    MAX1map = cell(size(SUM1mapFind));
    M1Trace = cell(size(SUM1mapFind));
    for iRes = 1:numResolution
        [MAX1map(iRes,:) M1Trace(iRes,:) M1RowShift M1ColShift M1OriShifted] = ...
             mexc_ComputeMAX1( numOrient, SUM1mapFind(iRes,:), locationShiftLimit,...
             orientShiftLimit, 1 );
    end
    save(sprintf('working/ImageAndFeature_%d.mat',img),'MAX1map','M1Trace','M1RowShift','M1ColShift','M1OriShifted','-append');
end

disp(['filtering time: ' num2str(toc) ' seconds']);
