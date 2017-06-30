function para = config()

para.name = 'cluster_5_toy'; % Choose configuration name 'cluster_5_toy', 'cluster_5_full', 'cluster_5_base', 'cluster_11_full', 'cluster_11_base'
para.noWorkers = 8;  % Set number of workers, which depends on how many cores in your cpu  


para.dataPath = 'dataset/AnimalFace/';
para.categoryNames = {'BearHead','CatHead','ChickenHead','CowHead', ...
                  'DeerHead','DogHead','DuckHead','EagleHead','ElephantHead', ...
                  'HumanHead','LionHead','MonkeyHead',...
                  'MouseHead','PandaHead','PigeonHead','PigHead','RabbitHead','SheepHead', ...
                  'TigerHead','WolfHead'};
        

%% 5 clusters (toy)

if strcmp(para.name, 'cluster_5_toy')
    
    para.categoryNames = {'CatHead'};
    
    para.task_ids = 1;
    para.numCluster = 5;
    para.numEMIteration = 10;
    para.numResolution = 3;   
    
end

%% 5 clusters

if strcmp(para.name, 'cluster_5_full')
    
    para.task_ids = 1:20;
    para.numCluster = 5;
    para.numEMIteration = 10;
    para.numResolution = 3;   
    
end

%% 11 cluster

if strcmp(para.name, 'cluster_11_full')
    
    para.task_ids = 1:20;
    para.numCluster = 11;
    para.numEMIteration = 10;
    para.numResolution = 3;  

end
