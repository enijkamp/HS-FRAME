clear;
close all;

noWorkers = 4;
seed = 1;
rng(seed); % set the seed so that experiments are reproducible
read_feature_from_data = false;
run_svm = true;
% if false, use SUM2 instead

imgPath = '../AnimalFace/';
codebook_path = './output/codebook_old.mat';
categoryNames = {'BearHead','CatHead','ChickenHead','CowHead', ...
                 'DeerHead','DogHead','DuckHead','EagleHead','ElephantHead', ...
                  'HumanHead','LionHead','MonkeyHead',...
                  'MouseHead','PandaHead','PigeonHead','PigHead','RabbitHead','SheepHead', ...
                  'TigerHead','WolfHead'};
             % categoryNames = {'toy'};
img_size = 90;
numClass = 1% length(categoryNames);
numIteration = 2;

trainImgs = cell(0);
testImgs = cell(0);
trainLabels = [];
testLabels = [];

for iClass = 1:numClass  
    
	imgList = dir(fullfile(imgPath, categoryNames{iClass},'*.jpg'));
	disp( ['category ',  num2str(iClass), ': ' categoryNames{iClass}, ' has ' num2str(length(imgList)) ' images']);
	for iImg= 1:length(imgList)
		img = imread(fullfile(imgPath, categoryNames{iClass},imgList(iImg).name));
		trainImgs{end+1}=img;
		trainLabels = [trainLabels; iClass];
	end
	imgList = dir(fullfile(imgPath, [categoryNames{iClass} '_test'],'*.jpg'));
	disp( ['category ',  num2str(iClass), ': ' categoryNames{iClass}, '_test has ' num2str(length(imgList)) ' images']);
	for iImg = 1:length(imgList)
		img = imread(fullfile(imgPath,[categoryNames{iClass} '_test'],imgList(iImg).name));
		testImgs{end+1}=img;
		testLabels = [testLabels; iClass];
	end
end

codeBook = [];
for iClass = 1:numClass
    template_file = sprintf('./template/hab_task%d_iter%d.mat', iClass, numIteration);
    load(template_file);
    codeBook = [codeBook ; clusters];
end

%save(sprintf('./output/AllcodeBook.mat'),'codeBook');
% extract features for training and testing images

if ~isempty(gcp('nocreate')),delete(gcp('nocreate'));end
%parpool(noWorkers);

if read_feature_from_data && exist(['./output/trainFeatures.mat'],'file')
    disp(['====> Found cache data for train feature, skip extract train image.']);
else
    disp(['Extracting features for training images: noWorkers: ' num2str(noWorkers)]);
    batch_spm = cell(noWorkers, 1);
    batch_img = cell(noWorkers, 1);
    numImage = length(trainImgs);
    idx = floor(numImage / noWorkers * (0:noWorkers));
    idx(noWorkers+1) = numImage;
    for batch = 1:noWorkers
        batch_img{batch} = trainImgs(idx(batch) + 1 : idx(batch+1));
        batch_spm{batch} = BatchExtractFeatures(batch_img{batch},codeBook, idx(batch) + 1 : idx(batch+1));
    end
    trainFeatures = [];
    for batch = 1:noWorkers
        trainFeatures = [trainFeatures ; batch_spm{batch}];
    end
    save('./output/trainFeatures.mat', '-v7.3', 'trainFeatures');
    clear('trainFeatures');
end
if read_feature_from_data && exist(['./output/testFeatures.mat'],'file')
    load('./output/testFeatures.mat');
    disp(['====> Found cache data for train feature, using it.']);
else
    disp(['Extracting features for testing images: noWorkers: ' num2str(noWorkers)]);
    batch_spm = cell(noWorkers, 1);
    batch_img = cell(noWorkers, 1);
    numImage = length(testImgs);
    idx = floor(numImage / noWorkers * (0:noWorkers));
    idx(noWorkers+1) = numImage;
    for batch = 1:noWorkers
        batch_img{batch} = testImgs(idx(batch) + 1 : idx(batch+1));
        batch_spm{batch} = BatchExtractFeatures(batch_img{batch},codeBook, idx(batch) + 10001 : idx(batch+1) + 10000);
    end

    testFeatures = [];
    for batch = 1:noWorkers
        testFeatures = [testFeatures ; batch_spm{batch}];
    end
    save('./output/testFeatures.mat', '-v7.3', 'testFeatures');
end

if run_svm
    
    load('./output/trainFeatures.mat');
    % train the model using muli-class svm
    addpath('./liblinear/matlab');

    libLinearOptions = ['-s 4 -c 10 -B 0'];
    model = train(trainLabels,sparse(trainFeatures),libLinearOptions);

    % test the model, and get confusion matrix
    [testLabelsHat, acc]= predict(testLabels,sparse(testFeatures),model);
    save('svm_result.mat', 'model', 'acc', 'testLabelsHat');

end
%{ 
Details of liblinar options
    -s type : set type of solver (default 1)
      for multi-class classification
             0 -- L2-regularized logistic regression (primal)
             1 -- L2-regularized L2-loss support vector classification (dual)
             2 -- L2-regularized L2-loss support vector classification (primal)
             3 -- L2-regularized L1-loss support vector classification (dual)
             4 -- support vector classification by Crammer and Singer
             5 -- L1-regularized L2-loss support vector classification
             6 -- L1-regularized logistic regression
             7 -- L2-regularized logistic regression (dual)
      for regression
            11 -- L2-regularized L2-loss support vector regression (primal)
            12 -- L2-regularized L2-loss support vector regression (dual)
            13 -- L2-regularized L1-loss support vector regression (dual)
    -c cost : set the parameter C (default 1)
    -p epsilon : set the epsilon in loss function of epsilon-SVR (default 0.1)
    -e epsilon : set tolerance of termination criterion
            -s 0 and 2
                    |f'(w)|_2 <= eps*min(pos,neg)/l*|f'(w0)|_2,
                    where f is the primal function and pos/neg are # of
                    positive/negative data (default 0.01)
            -s 11
                    |f'(w)|_2 <= eps*|f'(w0)|_2 (default 0.001) 
            -s 1, 3, 4 and 7
                    Dual maximal violation <= eps; similar to libsvm (default 0.1)
            -s 5 and 6
                    |f'(w)|_inf <= eps*min(pos,neg)/l*|f'(w0)|_inf,
                    where f is the primal function (default 0.01)
            -s 12 and 13\n"
                    |f'(alpha)|_1 <= eps |f'(alpha0)|,
                    where f is the dual function (default 0.1)
    -B bias : if bias >= 0, instance x becomes [x; bias]; if < 0, no bias term added (default -1)
    -wi weight: weights adjust the parameter C of different classes (see README for details)
    -v n: n-fold cross validation mode
    -q : quiet mode (no outputs)
%}

