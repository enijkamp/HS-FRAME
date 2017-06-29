% clear
clear
close all;

% code
addpath('./hs-frame');

% load
read_codebook = true;   
read_feature_from_data = true;
run_svm_params = false;
run_svm = true;

% config
para = config();
noWorkers = para.noWorkers;
useSUM3 = para.useSUM3;

% run
seed = 1;
rng(seed);

imgPath = para.dataPath;
codebook_path = ['output/' para.name '/model/codebook.mat'];
categoryNames = para.categoryNames;
img_size = 90;
numClass = length(para.task_ids);
numIteration = 10;
spm_numLayers = 4;
spm_threshold = 0;
para.numCategory = length(para.task_ids);

trainImgs = cell(0);
testImgs = cell(0);
trainLabels = [];
testLabels = [];
for iClass = para.task_ids
	imgList = dir(fullfile(imgPath, categoryNames{iClass},'*.jpg'));
	disp( ['category ',  num2str(iClass), ': ' categoryNames{iClass}, ' has ' num2str(length(imgList)) ' images']);
	for iImg= 1:length(imgList)
		img = imread(fullfile(imgPath, categoryNames{iClass},imgList(iImg).name));
		img = imresize(img,[img_size img_size]);
		trainImgs{end+1}=img;
		trainLabels = [trainLabels; iClass];
	end
	imgList = dir(fullfile(imgPath, [categoryNames{iClass} '_test'],'*.jpg'));
	disp( ['category ',  num2str(iClass), ': ' categoryNames{iClass}, '_test has ' num2str(length(imgList)) ' images']);
	for iImg = 1:length(imgList)
		img = imread(fullfile(imgPath,[categoryNames{iClass} '_test'],imgList(iImg).name));
		img = imresize(img,[img_size img_size]);
		testImgs{end+1}=img;
		testLabels = [testLabels; iClass];
	end
end

codeBook = [];
if read_codebook && exist(codebook_path,'file')
    load(codebook_path);
else
    for iClass = para.task_ids
        min_max3 = 1e8;
        select_seed = 0;
        for seed = 1:5
            template_file = sprintf(['output/' para.name '/template/template_task%d_seed%d_iter%d.mat'], iClass, seed, numIteration);
            if exist(template_file, 'file')
                load(template_file);
                max3 = sum(min(MAX3scoreAll, [], 2));
                if max3<min_max3
                    select_seed = seed;
                    min_max3 = max3;
                end
            end
        end
        if select_seed == 0, error('no template found!'),end

        template_file = sprintf(['output/' para.name '/template/template_task%d_seed%d_iter%d.mat'], iClass, select_seed, numIteration);
        load(template_file);
        codeBook = [codeBook ; clusters];

    end
    if ~exist(['output/' para.name '/model'],'dir'),mkdir(['output/' para.name '/model']),end
    save(codebook_path,'codeBook');
end

% extract features for training and testing images

if ~isempty(gcp('nocreate')),delete(gcp('nocreate'));end
parpool(noWorkers);

if read_feature_from_data && exist(['./output/' para.name '/trainFeatures.mat'],'file')
    disp(['====> Found cache data for train feature, skip extract train image.']);
else
    disp(['Extracting features for training images: noWorkers: ' num2str(noWorkers)]);
    batch_spm2 = cell(noWorkers, 1);
    batch_spm3 = cell(noWorkers, 1);
    batch_img = cell(noWorkers, 1);
    numImage = length(trainImgs);
    idx = floor(numImage / noWorkers * (0:noWorkers));
    idx(noWorkers+1) = numImage;
    parfor batch = 1:noWorkers
        batch_img{batch} = trainImgs(idx(batch) + 1 : idx(batch+1));
        [batch_spm2{batch}, batch_spm3{batch}] = BatchExtractFeatures(batch_img{batch},codeBook,para, idx(batch) + 1 : idx(batch+1));
    end

    trainFeatures = [];
    trainFeatures2 = [];
    trainFeatures3 = [];
    for batch = 1:noWorkers
        trainFeatures2 = [trainFeatures2 ; batch_spm2{batch}];
        trainFeatures3 = [trainFeatures3 ; batch_spm3{batch}];
    end
    save(['./output/' para.name '/trainFeatures.mat'], '-v7.3', 'trainFeatures2', 'trainFeatures3');
    clear('trainFeatures2');
    clear('trainFeatures3');
end
if read_feature_from_data && exist(['./output/' para.name '/testFeatures.mat'],'file')
    load(['./output/' para.name '/testFeatures.mat']);
    disp(['====> Found cache data for train feature, using it.']);
else
    disp(['Extracting features for testing images: noWorkers: ' num2str(noWorkers)]);
    numImage = length(testImgs);
    idx = floor(numImage / noWorkers * (0:noWorkers));
    idx(noWorkers+1) = numImage;
    parfor batch = 1:noWorkers
        batch_img{batch} = testImgs(idx(batch) + 1 : idx(batch+1));
        [batch_spm2{batch}, batch_spm3{batch}] = BatchExtractFeatures(batch_img{batch},codeBook,para, idx(batch) + 10001 : idx(batch+1) + 10000);
    end

    testFeatures = [];
    testFeatures2 = [];
    testFeatures3 = [];
    for batch = 1:noWorkers
        testFeatures2 = [testFeatures2 ; batch_spm2{batch}];
        testFeatures3 = [testFeatures3 ; batch_spm3{batch}];
    end
    save(['./output/' para.name '/testFeatures.mat'], '-v7.3', 'testFeatures2', 'testFeatures3');
end

if run_svm_params || run_svm
    
    load(['./output/' para.name '/trainFeatures.mat']);

    if useSUM3 
        trainFeatures = trainFeatures3;
        testFeatures = testFeatures3;
    else
        trainFeatures = trainFeatures2;
        testFeatures = testFeatures2;
    end

    addpath('./liblinear-multicore-2.11-1/matlab');
end

if run_svm_params
    
    svmResults = cell(1,6);
    i = 1;
    for s = [0,1,3]
        for c = 10:5:20
            for B = 0:10:50
                disp([]);
                libLinearOptions = ['-s ' num2str(s) ' -c ' num2str(c) ' -B ' num2str(B) ' -n ' num2str(noWorkers) ' -q'];
                model = train(trainLabels,sparse(trainFeatures),libLinearOptions);
                [testLabelsHat, acc, decision_values] = predict(testLabels,sparse(testFeatures),model);
                svmResults{i,1} = i;
                svmResults{i,2} = libLinearOptions;
                svmResults{i,3} = acc;
                svmResults{i,4} = model;
                svmResults{i,5} = testLabelsHat;
                svmResults{i,6} = decision_values;
                disp([num2str(i) ' acc ' num2str(acc(1)) ' options ' libLinearOptions]);
                i = i + 1;
            end
        end
    end
    
    save(['output/' para.name '/svmResults.mat'], 'svmResults');
    
    for i = 1:size(svmResults)
        display([num2str(i) ' ' num2str(svmResults{i,3}(1)) ' ' svmResults{i,2}]);
    end
    
end

if run_svm

    libLinearOpt = ['-s 1 -c 5 -B 0 -n ' num2str(noWorkers) ' -q'];
    model = train(trainLabels,sparse(trainFeatures),libLinearOpt);
    [testLabelsHat, acc, decision_values] = predict(testLabels,sparse(testFeatures),model);
    disp(['acc ' num2str(acc(1)) ' {options ' libLinearOpt '}']);
    
    save(['output/' para.name '/svm_accuracy.mat'], 'acc');
    save(['output/' para.name '/svm_model.mat'], 'model');

    addpath('./piotrs-toolbox-3.50/classify');
    addpath('./piotrs-toolbox-3.50/matlab');
    T = zeros(max(testLabels),length(testLabels));
    Y = T;
    for iSample = 1:length(testLabels)
            T(testLabels(iSample),iSample)=1;
            Y(testLabelsHat(iSample),iSample)=1;
    end
    [~,confMat,~,~] = confusion(T,Y);
    cm=confMat;
    for i=1:size(cm,1)
        cm(i,:)=cm(i,:)/sum(cm(i,:));
    end
    h=figure;
    confMatrixShow(cm, categoryNames,{'FontSize',10},[],0);
    print(h, '-dpdf', ['output/' para.name '/confusion.pdf']);
    saveas(h, ['output/' para.name '/confusion.png'], 'png');
    
    disp('done.');
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

