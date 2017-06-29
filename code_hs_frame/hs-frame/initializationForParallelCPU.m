function [nGPU]=initializationForParallelCPU(runParallel) 

if runParallel && (gpuDeviceCount<2)
    runParallel=0;
    disp(['Parallel computing will not be used because there are no mutiple GPUs setting in this machine, it has been automatically switched to the normal computing mode.']);
end

isOpen = matlabpool('size') > 0;
isOpenCorrect = matlabpool('size') == gpuDeviceCount;

nGPU=gpuDeviceCount;

if runParallel 
    nGPU=gpuDeviceCount;
    disp(['Parallel computing with ' num2str(nGPU) ' GPUs starts' ]);
    if ~isOpenCorrect,
      matlabpool close force
      matlabpool(nGPU)
    end
else
    nGPU=1;
    if isOpen,
      matlabpool close force
    end
end

end