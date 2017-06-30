function para = config()

para.name = 'cluster_5'; % Choose configuration name 'cluster_5_toy', 'cluster_5', 'cluster_11'
para.noWorkers = 8;  % Set number of workers, which depends on how many cores in your cpu  


para.dataPath = 'dataset/AnimalFace/';
para.categoryNames = {'BearHead','CatHead','ChickenHead','CowHead', ...
                  'DeerHead','DogHead','DuckHead','EagleHead','ElephantHead', ...
                  'HumanHead','LionHead','MonkeyHead',...
                  'MouseHead','PandaHead','PigeonHead','PigHead','RabbitHead','SheepHead', ...
                  'TigerHead','WolfHead'};
        

%% 5 clusters (toy)

if strcmp(para.name, 'cluster_5_toy')
    
    para.categoryNames = {'BearHead','CatHead'};
    
    para.task_ids = 1:2;
    para.numCluster = 5;
    para.numResolution = 3;   
    
    % to be frequently adjusted:
    para.numIteration = 10;  % number of iterations
    para.partRotationRange = 2*(-2:2); % absolute part rotation (rotation of partial templates)
    para.maxPartRelativeRotation = 2;
    para.resolutionShiftLimit = 1;
    para.rotationRange = 2*(-1:1); % whole object rotation


    % to be occationally adjusted
    para.numElement = 300; % number of Gabors in active basis
    para.locationPerturbFraction = .2; % part perturbation
    para.locationShiftLimit = 2; % shift in normal direction = locationShiftLimit*subsample pixels
    para.orientShiftLimit = 1; % shift in orientation
    para.subsampleS2 = 3;
    para.subsampleM2 = 1;
    
end

%% 5 clusters

if strcmp(para.name, 'cluster_5')
    
    para.task_ids = 1:20;
    para.numCluster = 5;
    para.numResolution = 3;   
    
    % to be frequently adjusted:
    para.numIteration = 10;  % number of iterations
    para.partRotationRange = 2*(-2:2); % absolute part rotation (rotation of partial templates)
    para.numPartRotate = length(para.partRotationRange);
    para.maxPartRelativeRotation = 2;
    para.resolutionShiftLimit = 1;
    para.minRotationDif = (sin(para.maxPartRelativeRotation*pi/para.numOrient)-sin(0))^2 + (cos(para.maxPartRelativeRotation*pi/para.numOrient)-cos(0))^2 + 1e-10;
    para.rotationRange = 2*(-1:1); % whole object rotation
    para.numRotate = length(para.rotationRange);

    % to be occationally adjusted
    para.numElement = 300; % number of Gabors in active basis
    para.locationPerturbFraction = .2; % part perturbation
    para.locationShiftLimit = 2; % shift in normal direction = locationShiftLimit*subsample pixels
    para.orientShiftLimit = 1; % shift in orientation
    para.subsampleS2 = 3;
    para.subsampleM2 = 1;
    
end

%% 11 cluster

if strcmp(para.name, 'cluster_11')
    
    para.task_ids = 1:20;
    para.numCluster = 11;
    para.numResolution = 3;   
    
    % to be frequently adjusted:
    para.numIteration = 10;  % number of iterations
    para.partRotationRange = 2*(-2:2); % absolute part rotation (rotation of partial templates)
    para.numPartRotate = length(para.partRotationRange);
    para.maxPartRelativeRotation = 2;
    para.resolutionShiftLimit = 1;
    para.minRotationDif = (sin(para.maxPartRelativeRotation*pi/para.numOrient)-sin(0))^2 + (cos(para.maxPartRelativeRotation*pi/para.numOrient)-cos(0))^2 + 1e-10;
    para.rotationRange = 2*(-1:1); % whole object rotation
    para.numRotate = length(para.rotationRange);

    % to be occationally adjusted
    para.numElement = 300; % number of Gabors in active basis
    para.locationPerturbFraction = .2; % part perturbation
    para.locationShiftLimit = 2; % shift in normal direction = locationShiftLimit*subsample pixels
    para.orientShiftLimit = 1; % shift in orientation
    para.subsampleS2 = 3;
    para.subsampleM2 = 1;

end
