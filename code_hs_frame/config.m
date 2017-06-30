function para = config()

para.name = 'cluster_5_toy'; % Choose configuration name 'cluster_5_toy', 'cluster_5_full', 'cluster_5_base', 'cluster_11_full', 'cluster_11_base'
para.noWorkers = 64;  % Set number of workers, which depends on how many cores in your cpu  


para.dataPath = 'dataset/AnimalFace/';
para.categoryNames = {'BearHead','CatHead','ChickenHead','CowHead', ...
                  'DeerHead','DogHead','DuckHead','EagleHead','ElephantHead', ...
                  'HumanHead','LionHead','MonkeyHead',...
                  'MouseHead','PandaHead','PigeonHead','PigHead','RabbitHead','SheepHead', ...
                  'TigerHead','WolfHead'};
        

%% 5 clusters (toy)

if strcmp(para.name, 'cluster_5_toy')
    
    para.task_ids = 1:2;

    para.method = 'two_stage';  % 'two_stage' : matching pursuit learning,  'one_stage': generative boosting learning

    para.nPartCol = 2;
    para.nPartRow = 2;
    para.part_sx = 50;
    para.part_sy = 50;

    para.GaborScaleList = [0.7];
    para.DoGScaleList = [];

    para.numCluster = 5;
    para.numWavelet = 300;   % default 300, 370
    para.numEMIteration = 10;
    para.isLocalNormalize = true;

    para.relativePartRotationRange=1*(-1:1);
    para.relativePartLocationRange=1;

    para.resolutionShiftLimit = 1;
    para.rotateShiftLimit = 3;   % template rotation  from -rotateShiftLimit to rotateShiftLimit, eg. (1)-2:2 if rotateShiftLimit=2 (2)0 is without rotation
    para.ratioDisplacementSUM3=0;   % default=0. Compute all values in SUM3 map

    para.locationShiftLimit=3;  % arg-max
    para.orientShiftLimit=1;    % arg-max

    para.numResolution = 3;     
    
    para.useSUM3 = 1;

end

%% 5 clusters

if strcmp(para.name, 'cluster_5_full')
    
    para.task_ids = 1:20;

    para.method = 'two_stage';  % 'two_stage' : matching pursuit learning,  'one_stage': generative boosting learning

    para.nPartCol = 2;
    para.nPartRow = 2;
    para.part_sx = 50;
    para.part_sy = 50;

    para.GaborScaleList = [0.7];
    para.DoGScaleList = [];

    para.numCluster = 5;
    para.numWavelet = 300;   % default 300, 370
    para.numEMIteration = 10;
    para.isLocalNormalize = true;

    para.relativePartRotationRange=1*(-1:1);
    para.relativePartLocationRange=1;

    para.resolutionShiftLimit = 1;
    para.rotateShiftLimit = 3;   % template rotation  from -rotateShiftLimit to rotateShiftLimit, eg. (1)-2:2 if rotateShiftLimit=2 (2)0 is without rotation
    para.ratioDisplacementSUM3=0;   % default=0. Compute all values in SUM3 map

    para.locationShiftLimit=3;  % arg-max
    para.orientShiftLimit=1;    % arg-max

    para.numResolution = 3;   
    
    para.useSUM3 = 1;
    
end


%% 5 cluster (without parts / baseline)

if strcmp(para.name, 'cluster_5_base')
    
    para.task_ids = 1:20;

    para.method = 'two_stage';  % 'two_stage' : matching pursuit learning,  'one_stage': generative boosting learning

    para.nPartCol = 2;
    para.nPartRow = 2;
    para.part_sx = 50;
    para.part_sy = 50;

    para.GaborScaleList = [0.7];
    para.DoGScaleList = [];

    para.numCluster = 5;
    para.numWavelet = 300;   % default 300, 370
    para.numEMIteration = 10;
    para.isLocalNormalize = true;

    para.relativePartRotationRange=0*(-1:1);
    para.relativePartLocationRange=0;

    para.resolutionShiftLimit = 0;
    para.rotateShiftLimit = 3;   % template rotation  from -rotateShiftLimit to rotateShiftLimit, eg. (1)-2:2 if rotateShiftLimit=2 (2)0 is without rotation
    para.ratioDisplacementSUM3=0;   % default=0. Compute all values in SUM3 map

    para.locationShiftLimit=3;  % arg-max
    para.orientShiftLimit=1;    % arg-max

    para.numResolution = 3;    
    
    para.useSUM3 = 1;

end

%% 11 cluster

if strcmp(para.name, 'cluster_11_full')
    
    para.task_ids = 1:20;

    para.method = 'two_stage';  % 'two_stage' : matching pursuit learning,  'one_stage': generative boosting learning

    para.nPartCol = 2;
    para.nPartRow = 2;
    para.part_sx = 50;
    para.part_sy = 50;

    para.GaborScaleList = [0.7];
    para.DoGScaleList = [];

    para.numCluster = 11;
    para.numWavelet = 300;   % default 300, 370
    para.numEMIteration = 10;
    para.isLocalNormalize = true;

    para.relativePartRotationRange=1*(-1:1);
    para.relativePartLocationRange=1;

    para.resolutionShiftLimit = 1;
    para.rotateShiftLimit = 3;   % template rotation  from -rotateShiftLimit to rotateShiftLimit, eg. (1)-2:2 if rotateShiftLimit=2 (2)0 is without rotation
    para.ratioDisplacementSUM3=0;   % default=0. Compute all values in SUM3 map

    para.locationShiftLimit=3;  % arg-max
    para.orientShiftLimit=1;    % arg-max

    para.numResolution = 3; 
    
    para.useSUM3 = 1;

end

%% 11 cluster (without parts / baseline)

if strcmp(para.name, 'cluster_11_base')
    
    para.task_ids = 1:20;

    para.method = 'two_stage';  % 'two_stage' : matching pursuit learning,  'one_stage': generative boosting learning

    para.nPartCol = 2;
    para.nPartRow = 2;
    para.part_sx = 50;
    para.part_sy = 50;

    para.GaborScaleList = [0.7];
    para.DoGScaleList = [];

    para.numCluster = 11;
    para.numWavelet = 300;   % default 300, 370
    para.numEMIteration = 10;
    para.isLocalNormalize = true;

    para.relativePartRotationRange=0*(-1:1);
    para.relativePartLocationRange=0;

    para.resolutionShiftLimit = 0;
    para.rotateShiftLimit = 3;   % template rotation  from -rotateShiftLimit to rotateShiftLimit, eg. (1)-2:2 if rotateShiftLimit=2 (2)0 is without rotation
    para.ratioDisplacementSUM3=0;   % default=0. Compute all values in SUM3 map

    para.locationShiftLimit=3;  % arg-max
    para.orientShiftLimit=1;    % arg-max

    para.numResolution = 3; 
    
    para.useSUM3 = 1;

end
