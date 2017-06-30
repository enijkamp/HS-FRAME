imgPath = '../AnimalFace/';
categoryNames = {'BearHead','CatHead','ChickenHead','CowHead', ...
    'DeerHead','DogHead','DuckHead','EagleHead','ElephantHead', ...
    'HumanHead','LionHead','MonkeyHead',...
    'MouseHead','PandaHead','PigeonHead','PigHead','RabbitHead','SheepHead', ...
    'TigerHead','WolfHead'};
numClass = length(categoryNames);
for iClass = 1:numClass
	imgList = [dir(fullfile(imgPath, categoryNames{iClass},'*.jpg')); dir(fullfile(imgPath, categoryNames{iClass},'*.JPG'))];
	disp( ['category ',  num2str(iClass), ': ' categoryNames{iClass}, ' has ' num2str(length(imgList)) ' images']);
	imgList = [dir(fullfile(imgPath, [categoryNames{iClass} '_test'],'*.jpg')); dir(fullfile(imgPath, [categoryNames{iClass} '_test'],'*.JPG'))];
	disp( ['category ',  num2str(iClass), ': ' categoryNames{iClass}, '_test has ' num2str(length(imgList)) ' images']);
end