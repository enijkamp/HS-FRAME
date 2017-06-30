% clear
clear;
close all;

% code
addpath('./aot');
addpath('./liblinear-multicore-2.11-1/matlab');

% options
read_feature_from_data = true;
run_svm_params = false;
run_svm = true;
run_confusion = true;

% config
para = config();

% fix seed
rng(1);

% load
imgPath = para.dataPath;
categoryNames = para.categoryNames;
numClass = length(para.task_ids);

if read_feature_from_data && exist(['./output/' para.name '/HAB_All_Features.mat'],'file')
    load(['output/' para.name '/HAB_All_Features.mat'], 'trainFeatures', 'testFeatures', 'trainLabels', 'testLabels');
else
    trainLabels = [];
    testLabels = [];
    for iClass = para.task_ids
        imgList = [dir(fullfile(imgPath, categoryNames{iClass},'*.jpg')); dir(fullfile(imgPath, categoryNames{iClass},'*.JPG'))];
        disp( ['category ',  num2str(iClass), ': ' categoryNames{iClass}, ' has ' num2str(length(imgList)) ' images']);
        for iImg= 1:length(imgList)
            trainLabels = [trainLabels; iClass];
        end
        imgList = [dir(fullfile(imgPath, [categoryNames{iClass} '_test'],'*.jpg')); dir(fullfile(imgPath, [categoryNames{iClass} '_test'],'*.JPG'))];
        disp( ['category ',  num2str(iClass), ': ' categoryNames{iClass}, '_test has ' num2str(length(imgList)) ' images']);
        for iImg = 1:length(imgList)
            testLabels = [testLabels; iClass];
        end
    end

    big_train = [];
    big_test = [];

    for iClass = para.task_ids
        load(['./output/' para.name '/trainFeatures_' categoryNames{iClass} '.mat']);
        load(['./output/' para.name '/testFeatures_' categoryNames{iClass} '.mat']);
        big_train = [big_train ; trainFeatures];
        big_test = [big_test ; testFeatures];
    end
    trainFeatures = big_train;
    testFeatures = big_test;

    if length(trainLabels)~=size(trainFeatures, 1), error('train data error, please check!'), end
    if length(testLabels)~=size(testFeatures, 1), error('test data error, please check!'), end

    save(['./output/' para.name '/HAB_All_Features.mat'], '-v7.3', 'trainFeatures', 'testFeatures', 'trainLabels', 'testLabels');
end

if run_svm_params || run_svm  
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
end

if run_confusion

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

