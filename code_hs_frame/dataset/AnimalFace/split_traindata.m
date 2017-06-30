dictionarys = dir;
dictionarys = dictionarys(3:end);
rng(1);
for ff = 1:length(dictionarys)
    folder = dictionarys(ff);
    mkdir([folder.name '_test']);
    imgList = dir(fullfile(folder.name, '*.jpg'));
    num_train=round(length(imgList)/2);
    idx= randperm(length(imgList));
    
    for iImg = num_train +1 :length(imgList)
        movefile([folder.name '/' imgList(idx(iImg)).name], [folder.name '_test/' imgList(idx(iImg)).name])
    end
end