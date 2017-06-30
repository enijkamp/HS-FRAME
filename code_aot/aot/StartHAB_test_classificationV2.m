
clear
imgPath = '../AnimalFace/';
% codebook_path = './output/codebook_old.mat';
categoryNames = {'BearHead','CatHead','ChickenHead','CowHead', ...
    'DeerHead','DogHead','DuckHead','EagleHead','ElephantHead', ...
    'HumanHead','LionHead','MonkeyHead',...
    'MouseHead','PandaHead','PigeonHead','PigHead','RabbitHead','SheepHead', ...
    'TigerHead','WolfHead'};


%%%%%%%%%%%%%%%%%%%
numIteration = 2;  % number of iterations
codebook = [];
for iClass = 1:1
    template_file = sprintf('./template/hab_task%d_iter%d.mat', iClass, numIteration);
    load(template_file);
    codebook = [codebook ; clusters];
end



for iClass=1:length(categoryNames)
    
    categoryName = categoryNames{iClass};
    imageFolder = ['../AnimalFace/' categoryName]; % folder of training images
    imageName = dir([imageFolder '/*.jpg']);
    imgList = dir(fullfile(imgPath, categoryNames{iClass},'*.jpg'));
    
    for img= 1:length(imgList)
        
        disp([' start processing image ' num2str(img) ' in class ' categoryName]); tic       
        tmpIm = imread(fullfile(imgPath, categoryNames{iClass},imgList(img).name));       
        if size(tmpIm,3) == 3
            tmpIm = rgb2gray(tmpIm);
        end
        imagesBatch{img}=tmpIm;
    end  
     
    
    
    Features = BatchExtractFeaturesV2(codebook, imagesBatch);        
        
    
end  
