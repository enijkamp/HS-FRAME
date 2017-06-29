for task_id=1:16
    for seed=1:10
       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

         para = getParametersForExp(task_id);
         disp(['generating initialization matrix for task ' num2str(task_id) ' with seed ' num2str(seed)]);


         % set seed for random numbers generation
          rng(seed);

          inPath = ['./positiveImages/' para.foldName];
          files = dir(fullfile(inPath,'*.jpg'));
          numImage=length(files);
          numCluster=para.numCluster;
          
          MAX3scoreAll = rand(numImage, numCluster);   % randomly assign members to different cluster

          if ~exist([inPath '/initialization'],'dir')
               mkdir([inPath '/initialization']);
          end
          save([inPath '/initialization/weightingMatrix' num2str(seed) '.mat'],'MAX3scoreAll');
    end
end
